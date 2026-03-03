## Analysis: ComputeSourceWindowStartingFromSource in GDAL

### What it does

The function is the **forward probe** — it transforms a grid of source pixels
into destination pixel space and checks which ones land inside the current
destination window. It then updates the source bounding box to include only
the source pixels that actually contribute.

### Algorithm (lines 2620-2748)

1. **One-time setup** (cached in `privateData`, reused across calls):
   - Create a (nStepCount+2) x (nStepCount+2) grid across the **entire source raster**
   - The +2 adds points at the very edges (0.5 pixels in from each border)
   - Default nStepCount = 21, so 23x23 = 529 points
   - Transform ALL points forward: source pixel -> destination pixel
   - Cache the results (transform is expensive, do it once)

2. **Per-call filtering** (lines 2720-2747):
   - For each cached point, check if its destination coordinates fall inside
     the current destination window [nDstXOff, nDstXOff+nDstXSize] x
     [nDstYOff, nDstYOff+nDstYSize]
   - If yes, update the source min/max with that point's SOURCE coordinates
   - This tightens the source bounding box to only include source pixels
     whose forward transform actually lands in the destination tile

### Key insight: it MERGES, it does not trim

At line 3205, `ComputeSourceWindowStartingFromSource` is called with pointers
to dfMinXOut/dfMaxXOut that already contain the backward pass result. The
function updates these with `std::min`/`std::max` — so it can only EXPAND
the source window, never shrink it.

Lines 2740-2743:
```cpp
*padfSrcMinX = std::min(*padfSrcMinX, dfSrcX);
*padfSrcMinY = std::min(*padfSrcMinY, dfSrcY);
*padfSrcMaxX = std::max(*padfSrcMaxX, dfSrcX);
*padfSrcMaxY = std::max(*padfSrcMaxY, dfSrcY);
```

This EXPANDS the source window. The forward probe ADDS source pixels that the
backward probe might have missed (because the backward transform failed for
some destination points). It does NOT trim the source window.

So for the antimeridian case, this function would not help trim the full-width
read — it would only potentially expand it further.

### When it's called

Line 3203: only called when `bUseGrid == true`. The grid mode is triggered by:
- Failed corner transforms (any of the 4 dst corners fail backward transform)
- Special points in destination extent (poles)
- Failed edge sampling points
- `SAMPLE_GRID=YES` warp option

### The actual antimeridian "fix" in GDAL

The antimeridian is not fixed by `ComputeSourceWindowStartingFromSource`.
It is fixed by **CollectChunkListInternal** (lines 1456-1624):

1. `ComputeSourceWindow` returns a source window and a `dfSrcFillRatio`
2. If `dfSrcFillRatio < 0.5` AND destination is large (>100 pixels),
   **subdivide the destination** (line 1531-1536)
3. Split along the longest dimension (halve it), recurse on each half
4. Each half gets its own `ComputeSourceWindow` call
5. Eventually the halves are small enough that they are entirely on one side
   of the antimeridian, producing compact source windows

This is the recursive bifurcation approach. It does not detect the antimeridian
explicitly — it detects that the source fill ratio is low (a full-width source
window for a small destination means most source pixels are wasted) and
subdivides until the ratio improves.

### What this means for cogcache

1. `ComputeSourceWindowStartingFromSource` is a **safety net for missing
   source pixels** when the backward transform fails. It is not a source
   window optimisation. Implementing it adds robustness (catches pixels
   the backward probe misses) but does not solve the antimeridian read
   efficiency problem.

2. The antimeridian efficiency fix is `CollectChunkListInternal`'s recursive
   subdivision driven by `fill_ratio`. We already return `fill_ratio` from
   our `compute_source_window`. The R planning layer can use it:

   ```r
   sw <- rust_compute_source_window(...)
   if (sw$fill_ratio < 0.5 && tile_size > 16) {
     # subdivide destination and recurse
   }
   ```

3. For a pure Rust core that "does an efficient job emulating GDAL," we want:
   - `ComputeSourceWindowStartingFromSource` for correctness (captures
     source pixels missed by backward probe)
   - Recursive destination subdivision driven by fill_ratio for efficiency
   - Both are implementable in Rust as part of the source_window module

### Implementation plan for Rust

**Phase 1: Forward probe (ComputeSourceWindowStartingFromSource)**

New function in source_window.rs:

```rust
pub fn refine_source_window_from_source(
    transformer: &impl Transformer,
    src_raster_size: [i32; 2],
    dst_off: [i32; 2],
    dst_size: [i32; 2],
    current_sw: &mut SourceWindow,  // expanded in-place
    step_count: usize,
)
```

- Create a (step_count+2)^2 grid across the source raster
- Transform forward (dst_to_src = false)
- Filter: keep only points where dst coords fall in [dst_off, dst_off+dst_size]
- Expand current_sw to include those source points
- Cache the forward-transformed grid for reuse across tiles (this is the
  privateData pattern from GDAL)

Note: the forward transform goes source pixel -> destination pixel. Our
Transformer trait already supports dst_to_src = false.

**Phase 2: Recursive destination subdivision**

New function:

```rust
pub fn collect_chunk_list(
    transformer: &impl Transformer,
    src_raster_size: [i32; 2],
    dst_off: [i32; 2],
    dst_size: [i32; 2],
    resample_padding: i32,
    min_fill_ratio: f64,       // threshold, GDAL uses 0.5
    min_dst_size: i32,         // stop recursion, GDAL uses 100
) -> Vec<ChunkPlan>
```

Where:
```rust
pub struct ChunkPlan {
    pub src_window: SourceWindow,
    pub dst_off: [i32; 2],
    pub dst_size: [i32; 2],
}
```

- Compute source window for the full destination
- If fill_ratio < min_fill_ratio AND dst is large enough, split destination
  along longest dimension and recurse
- Return a list of ChunkPlan — each can be warped independently
- The R layer calls this once per tile and gets the optimal set of source reads

This gives the R layer a single call that returns the optimal set of
source reads for a destination tile, handling antimeridian splits and
memory budgets automatically.

### Special points (pole detection)

GDAL pre-computes where the poles (lat +/-89.9999) land in destination pixel
space (lines 695-712). If a pole falls within the destination raster extent,
it forces grid sampling mode (which triggers the forward probe).

The mechanism: transform (lon=0, lat=+/-89.9999) through
`GDALTransformLonLatToDestGenImgProjTransformer` to get destination pixel
coordinates. Store as `aDstXYSpecialPoints`. At ComputeSourceWindow time
(line 3094-3116), if any special point is inside the destination raster,
switch to grid mode.

For vwarp, we could generalise this: transform a set of known singularities
(poles, antimeridian points at various latitudes) to destination coordinates
and flag any destination tile that contains one. This would be a one-time
setup per destination grid, not per tile.

### Connection to the orthogonal grid lines discussion

The forward probe IS the "source grid lines in destination space" view.
The backward probe IS the "destination grid lines in source space" view.
GDAL combines both — backward first (cheap, usually sufficient), forward
as refinement (more expensive, catches edge cases).

The recursive subdivision is driven by the relationship between the two:
when the backward-derived source window is much larger than what the
forward probe says contributes, the fill ratio is low and subdivision
helps. This is exactly the "two-way insight" from the earlier discussion.
