# Rust Warper: Design from GDAL Internals

**Michael Sumner — March 2026**

## Provenance

GDAL source analysis based on the system build reported as:
`GDAL 3.13.0dev-a3b7b01d3e-dirty (released 2026-02-23)`

The commit `a3b7b01d3e` is the reference point. The warper files
(`gdalwarpkernel.cpp`, `gdaltransformer.cpp`, `gdalwarpoperation.cpp`) are
stable — last significant changes to the core warper logic predate this by
months/years. But we should pin the exact commit for reproducibility.

TODO: record the full commit hash from `git -C /path/to/gdal rev-parse HEAD`
on the machine where the analysis was done.

## Project framing

**cogcache** is a minimal alpha implementation. It works, it proves the
pipeline end-to-end, but it is not the foundation we build on. The goal is:

1. **A set of tight Rust crates** that reimplement the GDAL warp API as a
   lighter body of work that reflects and credits its origins. Community-adoptable,
   with clear lineage to the GDAL C++ source.

2. **DuckDB/Parquet for caching and config** — not now, but the architecture
   should not preclude it. Keep an eye on where GDAL's block cache and
   string-soup config would map to structured storage.

3. **Exposed in R first** — via extendr, because that's where we're comfortable
   and can iterate fastest. cogcache as a single R package is fine for now,
   as long as the Rust core is structured to be extractable as standalone
   crates later.

The Rust crate(s) should be usable without R. The extendr wrappers are a
thin skin over the real API.

## What we've proven

The full pipeline works end-to-end in Rust with pixel-identical results to
gdalraster's `pixel_extract`:

```
rust_warp_map()          → source pixel indices     ✅ 65536/65536 match
rust_fetch_decode_tile() → raw pixel values          ✅ pixel-identical to ds$read
rust_apply_warp()        → warped output             ✅ pixel-identical to pixel_extract
```

What remains: replacing the naive per-pixel-exact-PROJ approach with GDAL's
actual warper pipeline, which is faster and (for the general case) handles
edge cases we currently don't.

---

## GDAL's warp pipeline, decomposed

When you call `gdalwarp -r near src.tif dst.tif`, this is the actual call stack:

```
gdalwarp_lib.cpp:  GDALWarp()
  → GDALWarpOperation::ChunkAndWarpImage()
    → for each chunk:
        ComputeSourceWindow()         [Layer 3]
        WarpRegionToBuffer()
          → load source pixels
          → GDALWarpKernel::PerformWarp()
            → PerformWarp() dispatch
              → GWKNearestThread()    [Layer 1+2]
                → for each scanline:
                    pfnTransformer()  [Layer 5]
                    for each pixel:
                      sample source pixel
```

Six layers. For our tile-serving case (single 256×256 output tile from a COG),
three of these collapse:

- **Chunking** (Layer 4): unnecessary, output IS one chunk
- **Masking** (Layer 6): start with simple nodata value check
- **Config** (Layer 7): Rust struct, not string-soup

That leaves three layers that matter:

### Layer 1: Warp Kernel — `GWKNearestThread` / `GWKGeneralCaseThread`

The core loop. For nearest neighbour, the GDAL code (simplified) is:

```cpp
// gdalwarpkernel.cpp, GWKNearestThread (simplified)
for (int iDstY = 0; iDstY < nDstYSize; iDstY++) {
    // Fill arrays with dest pixel centres
    for (int iDstX = 0; iDstX < nDstXSize; iDstX++) {
        padfX[iDstX] = iDstX + 0.5 + nDstXOff;
        padfY[iDstX] = iDstY + 0.5 + nDstYOff;
        padfZ[iDstX] = 0;
    }

    // Transform whole scanline: dst pixel coords → src pixel coords
    pfnTransformer(pTransformerArg, TRUE, nDstXSize,
                   padfX, padfY, padfZ, pabSuccess);

    for (int iDstX = 0; iDstX < nDstXSize; iDstX++) {
        if (!pabSuccess[iDstX]) continue;

        // Nearest neighbour: floor
        int iSrcX = static_cast<int>(padfX[iDstX]);
        int iSrcY = static_cast<int>(padfY[iDstX]);

        // Bounds check against source window
        if (iSrcX < nSrcXOff || iSrcX >= nSrcXOff + nSrcXSize ||
            iSrcY < nSrcYOff || iSrcY >= nSrcYOff + nSrcYSize)
            continue;

        // Read source pixel, write to dest
        dst[iDstY * nDstXSize + iDstX] =
            src[(iSrcY - nSrcYOff) * nSrcXSize + (iSrcX - nSrcXOff)];
    }
}
```

**Rust equivalent**: this is exactly what `rust_apply_warp` does, minus the
per-scanline transformer call. Our current approach pre-computes all transforms
then applies. GDAL interleaves transform and sample per scanline.

**Design choice**: GDAL's per-scanline approach is better for the
ApproxTransformer (see Layer 5) because it can check error per scanline.
Our Rust version should also work per-scanline.

### Layer 2: Resampling Kernels

For nearest neighbour, it's just `floor()` (actually `static_cast<int>()` which
truncates toward zero — same as floor for positive values).

For bilinear/cubic/lanczos, the kernel needs fractional source coordinates
and reads a neighbourhood of source pixels. The key functions:

```
Bilinear: 2×2 neighbourhood, weights = (1 - dx) * (1 - dy) etc.
Cubic:    4×4 neighbourhood, Mitchell-Netravali weights
Lanczos:  6×6 neighbourhood, sinc windowed weights
```

**Rust design**: trait-based, generic over data type:

```rust
trait Resample {
    fn radius() -> usize;  // 0 for NN, 1 for bilinear, 2 for cubic, 3 for lanczos
    fn weight(dx: f64) -> f64;
    fn sample<T: Pixel>(src: &RasterBuffer<T>, x: f64, y: f64) -> T;
}

struct NearestNeighbour;
struct Bilinear;
struct Cubic;
struct Lanczos;
```

Monomorphisation gives us the specialisation GDAL achieves with 30 separate
functions, but from one generic implementation.

### Layer 3: Source Window Computation — `ComputeSourceWindow`

Given a destination region, compute which source pixels are needed. GDAL's
algorithm:

1. Sample points along the **edges** of dest region (21 per edge by default)
2. Transform each to source pixel coordinates via the transformer
3. Take bounding box of successful transforms
4. Add padding for kernel radius
5. Clamp to source raster bounds

**Edge cases** (where GDAL's 20 years of fixes live):

- **Grid sampling**: when edge-only misses interior extrema (polar projections).
  Enabled when corner transforms fail.
- **Forward transform fallback**: when inverse is unreliable, also samples
  source→dest to find overlap.
- **Antimeridian/pole handling**: special points in dest space that force
  grid sampling.

**Rust design**: for tile serving, the dest region is always small (256×256)
and the transform is smooth (Mercator↔UTM). Edge sampling with 4 corners +
edge midpoints is sufficient. Grid sampling as opt-in for pathological cases.

```rust
fn compute_source_window(
    dst_gt: &GeoTransform,
    dst_region: &[usize; 4],  // xoff, yoff, xsize, ysize
    transformer: &impl Transformer,
    kernel_radius: usize,
) -> Option<SourceWindow> {
    // Transform edge sample points
    // Take bounding box
    // Add kernel padding
    // Clamp to source bounds
}
```

### Layer 5: Coordinate Transformer — `GDALGenImgProjTransform`

This is the composition that maps dst pixels → src pixels:

```
dst pixel (col, row)
  → dst geo (x, y)     via dst geotransform
  → src geo (x, y)     via PROJ (CRS transform)
  → src pixel (col, row) via inverse src geotransform
```

Three steps, all well-defined. The geotransform steps are trivial affine math.
The CRS step is PROJ.

**ApproxTransformer** wraps this and adds interpolation:

1. Given N points on a scanline, transform only the endpoints + midpoint
2. Check: is the true midpoint close to the linearly interpolated midpoint?
3. If error < threshold (default 0.125 pixels): interpolate all interior points
4. If error > threshold: recursively subdivide and retry

This is a classic adaptive refinement. For smooth transforms (UTM↔Mercator),
almost all scanlines pass on first check → massive speedup (transform 3 points
instead of 256).

**Rust design**: this is the piece that makes a real warper vs our proof of concept.

```rust
trait Transformer {
    fn transform(&self, dst_to_src: bool, x: &mut [f64], y: &mut [f64]) -> Vec<bool>;
}

/// Exact per-pixel transform via PROJ
struct ProjTransformer {
    proj: proj::Proj,
    src_gt: GeoTransform,
    dst_gt: GeoTransform,
}

/// Adaptive interpolating wrapper (GDAL's ApproxTransformer)
struct ApproxTransformer<T: Transformer> {
    inner: T,
    max_error: f64,  // pixels
}
```

The `ApproxTransformer` is where the smoothing comes from that we see in
GDAL's output. It's not bilinear resampling of pixels — it's bilinear
interpolation of the *coordinate transform itself*. At 155m dest pixels
mapping to 10m source pixels, a 0.125 pixel error in source coordinates
means the sampled source pixel can shift by ~2 source pixels. This explains
our observed differences.

---

## Proposed Rust module structure

```
cogcache-core/
├── src/
│   ├── lib.rs
│   ├── grid.rs           // GridSpec, GeoTransform, extent_dim_to_gt, etc.
│   │                     // Mirrors vaster semantics
│   ├── transform.rs      // Transformer trait, ProjTransformer, ApproxTransformer
│   ├── kernel.rs         // Resample trait, NN/Bilinear/Cubic/Lanczos
│   ├── source_window.rs  // compute_source_window (from GDAL Layer 3)
│   ├── warp.rs           // warp_tile() — the main entry point
│   ├── decode.rs         // DEFLATE/JPEG/WebP tile decoders (existing)
│   └── fetch.rs          // HTTP range fetch (existing, → async later)
```

### `grid.rs` — the foundation

```rust
pub struct GridSpec {
    pub crs: String,          // WKT or EPSG
    pub extent: [f64; 4],     // xmin, xmax, ymin, ymax (vaster convention)
    pub dim: [usize; 2],      // ncol, nrow
}

impl GridSpec {
    pub fn geotransform(&self) -> GeoTransform { ... }
    pub fn x_from_col(&self, col: f64) -> f64 { ... }
    pub fn y_from_row(&self, row: f64) -> f64 { ... }
    pub fn col_from_x(&self, x: f64) -> f64 { ... }
    pub fn row_from_y(&self, y: f64) -> f64 { ... }
    pub fn xy_from_cell(&self, cell: usize) -> (f64, f64) { ... }
    pub fn cell_from_xy(&self, x: f64, y: f64) -> Option<usize> { ... }
}
```

### `warp.rs` — the main entry point

```rust
pub fn warp_tile<T: Pixel, R: Resample>(
    src: &GridSpec,
    dst: &GridSpec,
    src_pixels: &[T],       // source pixel buffer (decoded tiles)
    src_window: &SourceWindow,
    transformer: &impl Transformer,
    nodata: Option<T>,
) -> Vec<T> {
    let mut output = vec![nodata.unwrap_or_default(); dst.dim[0] * dst.dim[1]];

    for row in 0..dst.dim[1] {
        // Build scanline of dest pixel centres
        let mut xs: Vec<f64> = (0..dst.dim[0])
            .map(|col| dst.x_from_col(col as f64))
            .collect();
        let mut ys: Vec<f64> = vec![dst.y_from_row(row as f64); dst.dim[0]];

        // Transform: dst geo → src geo (via PROJ or ApproxTransformer)
        let success = transformer.transform(true, &mut xs, &mut ys);

        for col in 0..dst.dim[0] {
            if !success[col] { continue; }

            // Convert src geo → src pixel
            let src_col = src.col_from_x(xs[col]);
            let src_row = src.row_from_y(ys[col]);

            // Resample
            let dst_idx = row * dst.dim[0] + col;
            output[dst_idx] = R::sample(src_pixels, src_window, src_col, src_row);
        }
    }

    output
}
```

This is the whole warp kernel. ~30 lines of Rust. It replaces GDAL's
9,283-line gdalwarpkernel.cpp because generics handle the type/resampling
specialisation, and we've separated I/O from computation.

---

## Implementation order

### Phase 1 (now): GenImgProjTransformer

This is the closest thing to what we already do — the back-transform-coords-to-source
approach. GDAL's `GDALGenImgProjTransform` composes three steps:

1. **dst pixel → dst geo**: apply dst geotransform (affine, trivial)
2. **dst geo → src geo**: PROJ coordinate transform (the expensive bit)
3. **src geo → src pixel**: apply inverse src geotransform (affine, trivial)

Our `rust_warp_map` already does this, but all at once for every pixel. The
GenImgProjTransformer does it **per scanline** and takes arrays of points.
This is the right interface because:

- It matches GDAL's transformer callback signature
- It's the natural unit for the ApproxTransformer to wrap
- Per-scanline means we can interleave transform and sample

Implement as:

```rust
trait Transformer {
    fn transform(&self, dst_to_src: bool,
                 x: &mut [f64], y: &mut [f64]) -> Vec<bool>;
}

struct GenImgProjTransformer {
    proj: proj::Proj,       // dst CRS → src CRS
    src_gt: [f64; 6],       // source geotransform
    src_inv_gt: [f64; 6],   // inverse source geotransform
    dst_gt: [f64; 6],       // dest geotransform
    dst_inv_gt: [f64; 6],   // inverse dest geotransform
}
```

The `transform(dst_to_src=true, ...)` call:
- input: dst pixel coords (col + 0.5, row + 0.5)
- applies dst_gt to get dst geo coords
- calls PROJ to get src geo coords
- applies src_inv_gt to get src pixel coords
- output: src pixel coords (fractional)

This is exactly what GDAL does in `GDALGenImgProjTransform()` at
gdaltransformer.cpp ~line 2800.

### Phase 2: ApproxTransformer

Wraps GenImgProjTransformer with adaptive scanline interpolation:

1. Transform endpoints + midpoint of scanline exactly
2. Check interpolation error at midpoint
3. If error < threshold: linearly interpolate all interior points
4. If error > threshold: recursively subdivide

This is the piece that makes GDAL fast for tile serving and also explains
the pixel-level differences we see.

### Phase 3: Per-scanline warp kernel

Restructure the warp loop to work per-scanline instead of pre-computing
all transforms. This naturally combines with the transformer:

```rust
for row in 0..dst_nrow {
    let (xs, ys, ok) = transformer.transform_scanline(row);
    for col in 0..dst_ncol {
        if ok[col] { output[row * ncol + col] = resample(src, xs[col], ys[col]); }
    }
}
```

### Phase 4: Resampling kernels

Bilinear, cubic, lanczos — pure math, generic over pixel type.

### Phase 5: Source window computation

Edge sampling + kernel padding. Needed for efficient tile fetching.

### Phase 6: Async I/O + tile assembly in Rust

Move the tile fetch/decode/assemble loop from R into Rust.
Replace ureq (sync) with reqwest (async) for parallel tile fetches.

### Phase 7: Extract standalone crates

- `warp-core`: grid, transform, kernel, source window (no I/O)
- `warp-cog`: COG-specific fetch + decode + warp
- `cogcache`: R package wrapping the above via extendr

---

## The key insight from our comparison

Our exact per-pixel PROJ transform gives **mathematically more correct** results
than GDAL's ApproxTransformer — we match `pixel_extract` perfectly while GDAL's
`warp` disagrees with both. But GDAL's approach is:

1. **Faster** — transforms ~3 points per scanline instead of 256
2. **Standard** — everyone expects GDAL's output, warts and all
3. **Configurable** — `-et` controls accuracy/speed tradeoff

For a production warper, we should support both modes:
- `TransformMode::Exact` — current behaviour, pixel-identical to `pixel_extract`
- `TransformMode::Approximate(tolerance)` — GDAL-like adaptive interpolation

The user picks based on their needs: exact for science, approximate for tile serving.
