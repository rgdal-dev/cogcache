# Design Notes & Known Edge Cases

**cogcache warp pipeline — living document**

This file records design decisions, hard calls, and known limitations in the
Rust warp implementation. It's intended for anyone reading or depending on
this code. If something here surprises you, that's the point — better to
find out here than in a silent 1-pixel shift.

Last updated: 2026-03-02

GDAL source reference: commit `a3b7b01d3e` (3.13.0dev, 2026-02-23).
TODO: record full hash from `git -C /path/to/gdal rev-parse HEAD`.

---

## Things that map cleanly

**Grid arithmetic** (`grid.rs`). Affine math — identical semantics in
R/vaster, GDAL, and Rust. Cell centres, extent↔geotransform conversion,
pixel↔geo coordinate mapping. No ambiguity here.

**The three-step GenImgProjTransformer** (`transform.rs`). Composition of
dst pixel → dst geo → src geo → src pixel. Each step is well-defined.
The geotransform steps are trivial affine multiplies. The CRS step is PROJ.
Maps directly from `GDALGenImgProjTransform()` in gdaltransformer.cpp.

**The nearest-neighbour warp kernel** (`warp.rs`). Per-scanline loop: build
coords, transform, truncate, sample. Maybe 50 lines of Rust replacing a
function that's ~270 lines in GDAL (because GDAL inlines mask/density
handling we don't need yet).

---

## Hard calls we've made

### 1. Truncation vs round for nearest neighbour

**Decision**: truncation (`as i32`), matching GDAL's `static_cast<int>()`.

**Context**: Our earlier prototype used `.round()` and matched gdalraster's
`pixel_extract()` at 65536/65536 pixels. GDAL's warp kernel uses truncation,
which picks a different pixel when the fractional coordinate is exactly X.5.

**Risk**: Switching to truncation may break our perfect match with
`pixel_extract` while getting closer to GDAL `warp`. The test suite checks
both comparisons, so we'll know immediately if this matters.

**What to watch for**: If truncation disagrees with `pixel_extract`, we need
to decide which is our ground truth. For science (exact coordinates), round
is arguably more correct. For GDAL compatibility, truncation is right. We may
want to support both via the resampling trait.

### 2. The +0.5 cell-centre convention

**Decision**: Follow GDAL's convention — the caller passes `col + 0.5` as
the pixel coordinate, and the transformer applies the geotransform directly.

**Context**: In vaster/R, `x_from_col()` bakes in the +0.5 offset. In GDAL's
warp kernel, the caller explicitly constructs `iDstX + 0.5` and the
transformer treats it as a raw coordinate. We follow GDAL here because the
transformer is a general-purpose function that shouldn't assume cell centres.

**Risk**: If someone feeds 0-based integer column indices to the transformer
without adding 0.5, they'll get coordinates for the top-left corner of the
pixel, not the centre. This is a half-pixel shift — it's silent, it's wrong,
and it's exactly the kind of bug that takes days to find.

**Mitigation**: `transform_scanline()` in transform.rs handles the +0.5
addition, so callers using the high-level API don't need to think about it.
The raw `Transformer::transform()` is for cases where you know what you're
doing.

### 3. Two PROJ objects instead of one

**Decision**: Create separate `Proj` instances for dst→src and src→dst.

**Context**: GDAL uses a single `OGRCoordinateTransformation` and calls it
with a direction flag. The `proj` Rust crate's `Proj::convert()` doesn't
expose a direction argument — it always transforms in the "forward" direction
of the pipeline it was created with.

**Impact**: Doubles PROJ initialisation cost. Irrelevant for tile serving
(create once, reuse across tiles). Could matter if creating thousands of
short-lived transformers.

**Future**: The `proj` crate may add direction support. Or we could use the
raw `proj-sys` FFI to call `proj_trans()` with `PJ_FWD`/`PJ_INV` directly.
Not worth the complexity now.

### 4. No rotated geotransform support

**Decision**: `inv_geotransform()` returns `None` for rotated grids
(gt[2] != 0 or gt[4] != 0).

**Context**: GDAL's `GDALInvGeoTransform()` handles the full 2×3 affine
matrix inversion. We only handle the non-rotated case because COGs and
standard projected data are always north-up.

**Risk**: Silently rejecting rotated grids is the right failure mode (better
than producing wrong pixel coordinates). But someone working with airborne
imagery or certain HDF products might hit this. The error message should be
clear about why.

**Future**: Full 2×2 matrix inversion is trivial to add — it's just a
2×2 determinant and adjugate. Deferred until we have a test case.

---

## Resolved: bilinear/cubic/lanczos diffs vs GDAL when downsampling

### 9. GDAL kernel scaling for downsampled warps

**Status**: Understood, not a bug. Documented 2026-03-03.

**Symptom**: At 4.5:1 source/destination ratio, bilinear max_diff=1672
against GDAL. At 1.3:1 ratio, 65536/65536 exact match (max_diff=0).
Cubic and Lanczos show the same pattern with slightly larger diffs
(1766, 1870).

**Root cause**: When the source-to-destination pixel ratio exceeds ~1,
GDAL's warp kernel scales the filter support by the local downsampling
factor (`dfXScale`/`dfYScale` in gdalwarpkernel.cpp ~L3800). For
bilinear at 4.5:1, the kernel effectively covers ~9×9 source pixels
instead of 2×2. This is antialiased resampling — it prevents aliasing
artefacts by averaging over the source pixels that contribute to each
destination pixel. Our kernels always use the textbook filter width
(2×2 for bilinear, 4×4 for cubic, 6×6 for lanczos).

**Key insight**: This is NOT the same as overview selection. `-ovr NONE`
controls which overview band to read from (always full-res with NONE).
The kernel scaling is a separate mechanism that operates on full-res
source pixels with a wider convolution. There is no gdalwarp flag to
disable it — it's intrinsic to the warp kernel when downsampling.

When overviews exist, `-ovr AUTO` selects an overview close to the
output resolution so the kernel scaling stays near 1:1. When overviews
don't exist (or are disabled), the kernel scaling does the equivalent
work. Same output quality, different path.

**Verification** (GEBCO 2024, Fiji LCC tile [4,3]):

| Ratio | Kernel   | Match     | Max diff | Notes                   |
|-------|----------|-----------|----------|-------------------------|
| 1.3:1 | bilinear | 65536/65536 | 0      | bit-identical to GDAL   |
| 4.5:1 | nearest  | 65536/65536 | 0      | NN unaffected by scaling|
| 4.5:1 | bilinear | 5543/65536  | 1672   | kernel scaling diff     |
| 4.5:1 | cubic    | 12654/65536 | 1766   | wider kernel, larger diff|
| 4.5:1 | lanczos  | 11185/65536 | 1870   | widest kernel, largest  |

Additional verification: coordinates confirmed identical across Rust
`proj` crate, sf, and reproj (zero diff at 100 sampled pixels).
Bilinear manually recomputed in R matches Rust output exactly (-5699).
Synthetic small-tif without overviews produces same GDAL output as
COG. GDAL's answer of -4030 at the max-diff pixel requires source
coordinates ~2 pixels from our landing point, consistent with the
wider kernel averaging over the surrounding bathymetric shelf.

**Options for matching GDAL's downsampled behaviour**:

1. **Select overviews first** (preferred for vwarp): the planning layer
   already knows the source/dst ratio from `compute_source_window`.
   Select the appropriate overview level before the kernel runs, keeping
   the kernel simple at ~1:1 ratio. Matches GDAL's `-ovr AUTO` path.

2. **Implement kernel scaling**: add the `dfXScale`/`dfYScale` logic
   from gdalwarpkernel.cpp. This widens filter support and renormalises
   weights based on local downsampling ratio. More complex kernel, but
   handles cases where overviews don't exist.

3. **Document and accept**: the kernel is correct for 1:1 and
   upsampling. For downsampling without overviews, results differ from
   GDAL. Users who need GDAL-identical output at reduced resolution
   should pre-select overviews.

---

## Known limitations (not yet addressed)

### 5. PROJ axis order

We rely on `Proj::new_known_crs()` to normalise axis order to
easting/northing (lon/lat). This is correct for EPSG codes. But:

- Raw PROJ strings may not get normalised
- WKT2 with explicit lat/lon axis order may surprise us
- GDAL uses `OGR_CT_FORCE_TRADITIONAL_GIS_ORDER` as a safety net

We haven't tested edge cases. For the common case (EPSG codes for UTM
and WebMercator), this works. For exotic CRS definitions passed as raw
WKT, there may be silent x/y swaps.

**Mitigation**: Always prefer EPSG codes over raw WKT/PROJ strings when
constructing transformers. Document this for users.

### 6. Truncation semantics for negative coordinates

Rust's `as i32` on a float truncates toward zero (same as C's
`static_cast<int>()`). For negative fractional source coords (destination
pixels that project outside the source image):

```
(-0.3f64) as i32  →  0    (same as C)
(-1.7f64) as i32  → -1    (same as C)
```

The bounds check in the warp kernel catches all negative indices, so this
is safe. But it means pixels at the very edge of the source image might
be included or excluded differently than with round. For nearest neighbour
at the margins, this produces at most a 1-pixel difference.

### 7. Antimeridian and poles (the big one)

GDAL has ~200 lines of special-case handling in `ComputeSourceWindow()`
for:

- **Antimeridian wrapping**: destination tiles that cross ±180° longitude
- **Pole singularities**: destination tiles that include a pole, where
  the inverse projection maps to a ring of source pixels
- **"Inside-out" projections**: where edge-only sampling misses interior
  extrema (e.g., WGS84 → polar stereographic)

Our code currently does nothing about any of these. For the Hobart/Antarctic
use case, this will matter when destination tiles approach the pole.

**Symptoms**: Source window computation requests the entire source image,
or misses needed source pixels entirely. Warped output will have gaps or
garbage near the pole.

**Plan**: Implement grid sampling (GDAL's `SAMPLE_GRID` option) as a
fallback for non-linear transforms. Detect the need for it by checking
whether corner transforms fail or produce unreasonably large source
windows. This is Phase 5 work.

### 8. proj.db location

The `proj` crate with `bundled_proj` builds its own libproj but doesn't
know where the system proj.db lives. This produces warnings like
`proj_create: Cannot find proj.db`. The transforms still work for common
CRS pairs (EPSG codes are handled without the database for basic cases).

**Fix**: Set `PROJ_DATA` environment variable before calling any Rust
transform functions:

```r
Sys.setenv(PROJ_DATA = "/usr/share/proj")
## or: Sys.setenv(PROJ_DATA = sf::sf_proj_search_paths()[1])
```

For production use, this should be set in the package `.onLoad()`.

---

## Test oracle

Our ground truth hierarchy:

1. **Mathematical correctness**: Does the coordinate transform produce the
   right geographic location? Verified against reproj/PROJ in R.

2. **gdalraster::pixel_extract**: Per-pixel exact PROJ transform + gdalraster's
   own geotransform arithmetic. Our "gold standard" for pixel identity.

3. **GDAL warp**: The production warper. Uses ApproxTransformer (interpolated
   coordinates), which produces different pixel assignments than exact
   per-pixel PROJ. We expect and accept these differences.

Current status (2026-03-02, Sentinel-2 B04 over Macquarie Island):
- Rust GenImgProjTransformer vs pixel_extract: **65536/65536 identical** ✅
- Rust per-scanline warp vs legacy warp_map: **65536/65536 identical** ✅
- Rust warp vs GDAL warp: 5753/65536 identical, max diff 65 (expected, ApproxTransformer)
- Rust scanline warp: 0.068s vs pixel_extract 0.413s (6× faster, exact per-pixel PROJ)

Note: `rust_gen_img_proj_transform` output differs from manual R pixel coord
computation by exactly 0.5 everywhere — this is the cell-centre convention.
The transformer expects `col + 0.5` input (GDAL convention), which
`transform_scanline()` provides. The R test computed fractional pixel coords
without the +0.5 offset. Confirmed as a test artefact, not a code bug.

---

## GDAL source cross-references

For anyone reading the GDAL source alongside this code:

| Our code | GDAL equivalent | File:line (commit a3b7b01d3e) |
|---|---|---|
| `grid::x_from_col` | `GDALApplyGeoTransform` | gdal_alg.h |
| `grid::inv_geotransform` | `GDALInvGeoTransform` | gdaltransformer.cpp |
| `GenImgProjTransformer` | `GDALGenImgProjTransform` | gdaltransformer.cpp ~L2800 |
| `GenImgProjTransformer::new` | `GDALCreateGenImgProjTransformer2` | gdaltransformer.cpp ~L1600 |
| `Transformer` trait | `GDALTransformerFunc` typedef | gdal_alg.h |
| `transform_scanline` | top of `GWKNearestThread` | gdalwarpkernel.cpp ~L5520 |
| `warp_nearest` | `GWKNearestThread` | gdalwarpkernel.cpp ~L5510 |
| `ApproxTransformer` | `GDALApproxTransform` | gdaltransformer.cpp ~L3500 |
| `compute_source_window` | `ComputeSourceWindow` | gdalwarpoperation.cpp ~L2751-3367 |
| `collect_chunk_list` | `CollectChunkListInternal` | gdalwarpoperation.cpp ~L1456-1624 |
| `find_discontinuity_range` | *(no GDAL equivalent — GDAL uses memory budget)* | — |
| `build_source_grid` | `ComputeSourceWindowStartingFromSource` setup | gdalwarpoperation.cpp ~L2656-2703 |
| `refine_from_source` | `ComputeSourceWindowStartingFromSource` per-call | gdalwarpoperation.cpp ~L2720-2747 |
