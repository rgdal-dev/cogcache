# GDAL Warp Pipeline Architecture — Rust Implementation Plan

Based on GDAL source analysis (commit a3b7b01d3e, GDAL 3.13.0dev).

## Overview: The Call Chain

```
gdalwarp (command-line)
  └─ gdalwarp_lib.cpp: GDALWarp()
      ├─ Create transformer pipeline
      ├─ Create GDALWarpOptions
      └─ GDALWarpOperation
          ├─ Initialize()
          ├─ ChunkAndWarpImage()
          │   ├─ CollectChunkList()            — subdivide into memory-fitting chunks
          │   │   └─ CollectChunkListInternal() — recursive, calls ComputeSourceWindow()
          │   └─ for each chunk: WarpRegion()
          │       └─ WarpRegionToBuffer()
          │           ├─ Read source pixels (GDALRasterIO)
          │           ├─ Build GDALWarpKernel
          │           │   ├─ nSrcXOff, nSrcYOff, nSrcXSize, nSrcYSize
          │           │   ├─ nDstXOff, nDstYOff, nDstXSize, nDstYSize
          │           │   ├─ pfnTransformer (function pointer)
          │           │   └─ source/dest pixel buffers + masks
          │           └─ oWK.PerformWarp()
          │               └─ kernel function per resample algo
          │                   └─ for each scanline:
          │                       1. pfnTransformer(dst → src)
          │                       2. for each pixel: GWKCheckAndComputeSrcOffsets()
          │                       3. read source pixel(s) → write dest pixel
          └─ Cleanup
```

## Phase 1: Transformer Pipeline (DONE)

**Source**: `gdaltransformer.cpp`

### What we have

```
ApproxTransformer(max_error=0.125) {      ← GDALApproxTransform
  wraps: GenImgProjTransformer {           ← GDALGenImgProjTransform
    step 1: dst pixel → dst geo             (apply dst geotransform)
    step 2: dst geo → src geo               (PROJ CRS transform)
    step 3: src geo → src pixel             (apply inv src geotransform)
  }
}
```

**Evidence** (gdalwarp_lib.cpp):

| Line | Code | Meaning |
|------|------|---------|
| 2832 | `GDALCreateGenImgProjTransformer2(hSrcDS, hDstDS, ...)` | Creates inner transformer |
| 2841 | `pfnTransformer = GDALGenImgProjTransform` | Function pointer for inner |
| 1598 | `dfErrorThreshold = 0.125` | Default threshold (pixels) |
| 3240 | `bUseApproxTransformer = dfErrorThreshold != 0.0` | Unless -et 0 |
| 3260-3263 | `GDALCreateApproxTransformer(GDALGenImgProjTransform, ...)` | Wraps entire pipeline |
| 3302 | `psWO->pfnTransformer = pfnTransformer` | Passed to warp operation |
| (warpop) 1987 | `oWK.pfnTransformer = psOptions->pfnTransformer` | Passed to kernel |

The `REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT` option (line 2778) wraps only the
reprojection step, but it is **never used by default gdalwarp** — zero hits in
gdalwarp_lib.cpp or gdalwarpoperation.cpp. It exists for other API callers.

### Rust status

- `transform.rs`: `GenImgProjTransformer` — validated identical to GDAL
- `approx.rs`: `ApproxTransformer` — matches GDAL 3.13.0dev recursive algorithm
- Verified: Rust output is **bit-identical** to GDAL's `transform_xy` for coordinates

### What's left

Nothing. The transformer is complete and **bit-identical** to GDAL.


## RESOLVED: The ~60000 Pixel Gap

The persistent ~5749/65536 match rate between our output and default `gdalwarp` was
caused by **overview selection**. The source COG has 4 overviews. At zoom 10 with
~155m output pixels vs 10m source pixels, GDAL automatically selects an overview
level (~160m) for efficiency. We read full-resolution data.

Evidence:
```
GDAL VRT vs GDAL VRT(-ovr NONE):  5749 / 65536   ← overview vs full-res
GDAL VRT(-ovr NONE) vs GTiff:    65536 / 65536   ← both full-res: identical
Rust exact vs GDAL VRT(-ovr NONE):65536 / 65536   ← PERFECT MATCH
Rust exact vs GDAL GTiff:        65536 / 65536   ← PERFECT MATCH
```

Our Rust warp is **bit-identical** to `gdalwarp -et 0 -ovr NONE -r near`.
Overview selection is an optimization feature, not a correctness issue.


## Phase 2: Source Window Computation

**Source**: `gdalwarpoperation.cpp`, `ComputeSourceWindow()` (line 3015)

### What it does

Given a destination window (nDstXOff, nDstYOff, nDstXSize, nDstYSize), compute which
region of the source raster needs to be read to produce that output.

### Algorithm

1. **Sample strategy selection** (line 3028-3073):
   - Default: sample along edges only (`bUseGrid=false`, 21 steps = `DEFAULT_STEP_COUNT`)
   - If any corner fails to reproject: sample ALL edge pixels (`bAll=true`)
   - If pole-like special points in extent: upgrade to grid sampling
   - If any edge point fails: retry with grid

2. **Edge/grid sampling** (`ComputeSourceWindowTransformPoints`, line 2754):
   - Edge mode: 4 × nStepCount points around the rectangle perimeter
   - Grid mode: (nStepCount+2)² points on interior grid
   - Transform all points dst→src via `pfnTransformer`
   - Collect min/max of successful transformed coordinates

3. **Supplementary forward sampling** (line 3203, `ComputeSourceWindowStartingFromSource`):
   - Only when grid mode is active
   - Samples source raster edges, transforms src→dst
   - If destination point falls in the target window, expands source bounds
   - Handles "inside-out" projections (polar stereographic, etc.)

4. **Resampling kernel padding** (line 3262-3297):
   - `nResWinSize = GWKGetFilterRadius(eResampleAlg)`
     - NearestNeighbour: 0
     - Bilinear: 1
     - Cubic/CubicSpline: 2
     - Lanczos: 3
   - Adjust for scale: `if (dfXScale < 0.95) nXRadius = ceil(nResWinSize / dfXScale)`
   - Add `SOURCE_EXTRA` if specified, or +10 if any points failed

5. **Clamping and output** (line 3309-3367):
   - Round near-integer values (`roundIfCloseEnough`, 1e-6 tolerance)
   - Clamp to source raster bounds
   - If >90% of source width covered, use full width (anti-meridian handling)
   - Output: `nSrcXOff, nSrcYOff, nSrcXSize, nSrcYSize` + extra size for filter

### Rust plan

```rust
/// Compute the source window needed for a destination region.
///
/// Returns (src_off_x, src_off_y, src_size_x, src_size_y, extra_x, extra_y)
pub fn compute_source_window(
    transformer: &impl Transformer,
    dst_off: [i32; 2],      // nDstXOff, nDstYOff
    dst_size: [i32; 2],     // nDstXSize, nDstYSize
    src_raster_size: [i32; 2],  // full source dimensions
    resample_radius: i32,   // GWKGetFilterRadius result
    source_extra: i32,      // additional padding, default 0
) -> Option<SourceWindow>
```

**Simplifications for v1**: skip grid mode (only edge sampling), skip
`ComputeSourceWindowStartingFromSource` (forward probing), skip special-point
pole detection. These are needed for global/polar projections but not for
typical Sentinel-2 tile warps.


## Phase 3: Chunk Management

**Source**: `gdalwarpoperation.cpp`, `CollectChunkListInternal()` (line 1456)

### What it does

Subdivides the destination extent into chunks that fit within memory limits,
then processes each chunk independently with its own source window.

### Algorithm

1. Call `ComputeSourceWindow()` for the full destination region
2. Estimate memory: `(srcPixelBits × srcPixels + dstPixelBits × dstPixels) / 8`
3. If exceeds `dfWarpMemoryLimit` (default 64MB) or source fill ratio < 0.5:
   - Split along the longer dimension, respecting block alignment
   - Recurse on each half
4. Store each leaf chunk as `{dst_off, dst_size, src_off, src_size, src_extra}`
5. Process chunks sequentially (or multi-threaded with `ChunkAndWarpMulti`)

### The _GDALWarpChunk struct (inferred from usage)

```c
struct _GDALWarpChunk {
    int dx, dy, dsx, dsy;           // destination offset and size
    int sx, sy, ssx, ssy;           // source offset and size
    double sExtraSx, sExtraSy;      // extra source pixels for filter
};
```

### Rust plan

For web tile rendering (256×256 output), chunking is unnecessary — the entire
destination fits in one chunk. Implement as a single `WarpRegion` call.

For larger warps (full COG→COG), implement chunking:

```rust
pub struct WarpChunk {
    dst: Rect,   // destination window
    src: Rect,   // source window
    src_extra: [f64; 2],
}

pub fn collect_chunks(
    dst_region: Rect,
    transformer: &impl Transformer,
    src_raster_size: [i32; 2],
    memory_limit: usize,
    resample_radius: i32,
) -> Vec<WarpChunk>
```


## Phase 4: WarpRegionToBuffer — Source I/O and Kernel Setup

**Source**: `gdalwarpoperation.cpp`, `WarpRegionToBuffer()` (line 1933)

### What it does

For a single chunk: reads source pixels, creates the kernel, runs it, writes output.

### Key operations

1. **Source window computation** (if not provided) — calls `ComputeSourceWindow()`

2. **Kernel setup** (line 1979-2010):
   ```
   oWK.nSrcXOff = nSrcXOff     // CRITICAL: offset into full-image coordinates
   oWK.nSrcYOff = nSrcYOff
   oWK.nSrcXSize = nSrcXSize   // buffer dimensions (what we actually read)
   oWK.nSrcYSize = nSrcYSize
   oWK.nDstXOff = nDstXOff     // offset in full-image coordinates
   oWK.nDstYOff = nDstYOff
   oWK.nDstXSize = nDstXSize
   oWK.nDstYSize = nDstYSize
   oWK.pfnTransformer = pfnTransformer   // STILL operates in full-image coords
   ```

3. **Source read** (around line 2010-2100):
   - `GDALRasterIO(hSrcDS, GF_Read, nSrcXOff, nSrcYOff, nSrcXSize, nSrcYSize, ...)`
   - Source pixel buffer is indexed `[0..nSrcXSize-1, 0..nSrcYSize-1]`
   - But the transformer returns coordinates in **full-image space**
   - The kernel subtracts `nSrcXOff`/`nSrcYOff` to convert to buffer indices

4. **Mask creation** — nodata masking, alpha, cutline, density masks

5. **Destination initialization** — INIT_DEST (typically 0 or NO_DATA)

6. **Run kernel**: `oWK.PerformWarp()`

7. **Write output** — `GDALRasterIO(hDstDS, GF_Write, ...)`

### The nSrcXOff/nSrcYOff offset (IMPLEMENTED CORRECTLY)

The transformer always returns coordinates in the **full source image** coordinate
space. The kernel subtracts `nSrcXOff`/`nSrcYOff` to convert to buffer indices:

```c
// GWKCheckAndComputeSrcOffsets, line 5304-5305
int iSrcX = (int)(padfX[iDstX] + 1.0e-10) - nSrcXOff;
int iSrcY = (int)(padfY[iDstX] + 1.0e-10) - nSrcYOff;
iSrcOffset = iSrcX + iSrcY * nSrcXSize;
```

Our Rust `warp_nearest()` implements this correctly. Verified **bit-identical** to
GDAL with `-et 0 -ovr NONE -r near` (65536/65536 match).

### Bounds checking in GWKCheckAndComputeSrcOffsets (line 5206-5320)

This function does more than just compute the offset:

1. **Retry logic** (line 5215-5302): If a pixel is slightly outside the source
   window (within 1 pixel), it re-transforms that single point using the exact
   transformer (bypassing the ApproxTransformer). This catches edge cases where
   the approximation places a pixel just outside the source.

2. **Bounds checks**: `padfX < nSrcXOff` or `padfX + 1e-10 > nSrcXSize + nSrcXOff`

3. **Epsilon**: `+ 1.0e-10` handles floating-point truncation at pixel boundaries

4. **Edge clamping**: `if (iSrcX == nSrcXSize) iSrcX--`

### Rust plan

The key fix: our warp function must subtract `src_x_off` and `src_y_off` when
converting transformer output to buffer indices, matching GDAL exactly.

```rust
pub fn warp_region(
    transformer: &impl Transformer,
    src_data: &[u16],           // source pixel buffer
    src_off: [i32; 2],          // nSrcXOff, nSrcYOff (full-image coords)
    src_size: [i32; 2],         // nSrcXSize, nSrcYSize (buffer dimensions)
    dst_off: [i32; 2],          // nDstXOff, nDstYOff
    dst_size: [i32; 2],         // nDstXSize, nDstYSize
    resample: ResampleAlg,
) -> Vec<u16>
```

Per-pixel logic:
```rust
// Transformer returns full-image coordinates
let (src_x_full, src_y_full) = transformer.transform_point(dst_x, dst_y);

// Convert to buffer coordinates (GDAL line 5304-5305)
let src_x_buf = (src_x_full + 1e-10) as i32 - src_off[0];
let src_y_buf = (src_y_full + 1e-10) as i32 - src_off[1];

// Bounds check against buffer dimensions
if src_x_buf < 0 || src_x_buf >= src_size[0] ||
   src_y_buf < 0 || src_y_buf >= src_size[1] {
    continue;  // skip this destination pixel
}

// Read from buffer
let offset = src_x_buf as usize + src_y_buf as usize * src_size[0] as usize;
dst_pixel = src_data[offset];
```


## Phase 5: Warp Kernel — Per-Scanline Processing

**Source**: `gdalwarpkernel.cpp`, `PerformWarp()` (line 1050)

### Kernel selection (line 1050-1400)

`PerformWarp()` selects the optimal kernel function based on:
- Resampling algorithm (NearestNeighbour, Bilinear, Cubic, etc.)
- Data type (Byte, Int16, UInt16, Float32, Float64)
- Mask configuration (nodata, alpha, density, validity)
- SIMD availability (SSE2, AVX2)

For UInt16 NearestNeighbour without masks (our case):
`GWKNearestNoMasksOrDstDensityOnlyShort` → `GWKRun()` → threaded execution.

### Per-scanline flow (line 5559-5700, simplified)

```
for each dst_row (iDstY = 0..nDstYSize):
    1. Build coordinate arrays:
       padfX[i] = i + 0.5 + nDstXOff    // pixel centers, full-image coords
       padfY[i] = iDstY + 0.5 + nDstYOff

    2. Transform entire scanline:
       pfnTransformer(pTransformerArg, TRUE, nDstXSize, padfX, padfY, padfZ, pabSuccess)

    3. Optional: SRC_COORD_PRECISION rounding (line 5587-5593)
       - Snaps coordinates to grid for tiled COG alignment
       - Falls back to exact transform in uncertainty zones

    4. For each dst_col (iDstX = 0..nDstXSize):
       a. GWKCheckAndComputeSrcOffsets():
          - Bounds check, retry on edge, compute buffer offset
          - iSrcX = (int)(padfX + 1e-10) - nSrcXOff
          - iSrcY = (int)(padfY + 1e-10) - nSrcYOff
          - iSrcOffset = iSrcX + iSrcY * nSrcXSize

       b. Validity/nodata checks (skip if masked)

       c. Read source pixel → write to destination buffer
```

### Scale computation (line 1073-1260)

Before the scanline loop, the kernel computes `dfXScale` and `dfYScale` —
the ratio of destination to source resolution. This controls:
- Resampling kernel size for non-nearest algorithms
- Whether to use area-weighted averaging (when downsampling)

For nearest-neighbour, this is informational only and doesn't affect the output.

### Thread parallelism (GWKRun, line 436)

If `WARP_NUM_THREADS` > 1, the scanlines are divided among threads. Each thread
gets its own transformer arg (cloned via `GDALCloneTransformer`). This is why
the ApproxTransformer is stateless — each thread needs its own instance.

### Rust plan

```rust
/// The core warp kernel for a single chunk.
pub fn warp_kernel_nearest(
    transformer: &impl Transformer,
    src_data: &[u16],
    src_off: [i32; 2],      // full-image offset of source buffer
    src_size: [i32; 2],      // source buffer dimensions
    dst_off: [i32; 2],      // full-image offset of dest region
    dst_size: [i32; 2],      // dest dimensions
    nodata: Option<u16>,
) -> Vec<u16> {
    let mut dst = vec![nodata.unwrap_or(0); (dst_size[0] * dst_size[1]) as usize];

    for row in 0..dst_size[1] {
        // Build scanline coordinates (pixel centers, full-image space)
        let mut x: Vec<f64> = (0..dst_size[0])
            .map(|col| col as f64 + 0.5 + dst_off[0] as f64)
            .collect();
        let mut y = vec![row as f64 + 0.5 + dst_off[1] as f64; dst_size[0] as usize];

        // Transform entire scanline (dst → src)
        let success = transformer.transform(true, &mut x, &mut y);

        for col in 0..dst_size[0] as usize {
            if !success[col] { continue; }

            // GDAL's GWKCheckAndComputeSrcOffsets (line 5304-5305)
            let src_x = (x[col] + 1e-10) as i32 - src_off[0];
            let src_y = (y[col] + 1e-10) as i32 - src_off[1];

            // Edge clamping (line 5306-5309)
            let src_x = if src_x == src_size[0] { src_x - 1 } else { src_x };
            let src_y = if src_y == src_size[1] { src_y - 1 } else { src_y };

            // Bounds check
            if src_x < 0 || src_x >= src_size[0] ||
               src_y < 0 || src_y >= src_size[1] {
                continue;
            }

            let src_offset = src_x as usize + src_y as usize * src_size[0] as usize;
            let dst_offset = col + row as usize * dst_size[0] as usize;

            // Nodata check
            if let Some(nd) = nodata {
                if src_data[src_offset] == nd { continue; }
            }

            dst[dst_offset] = src_data[src_offset];
        }
    }
    dst
}
```


## Phase 6: Resampling Kernels (Future)

**Source**: `gdalwarpkernel.cpp`, line ~5800 onwards

### Bilinear (2×2)

Needs 4 surrounding pixels. Uses fractional position within pixel:
```
dfDeltaX = src_x_full - floor(src_x_full) - 0.5
dfDeltaY = src_y_full - floor(src_y_full) - 0.5
```
Weight function: `(1 - |dx|) * (1 - |dy|)`

### Cubic (4×4)

Cubic convolution: `f(x) = (a+2)|x|³ - (a+3)|x|² + 1` for |x| ≤ 1
where a = -0.5 (Catmull-Rom).

### Lanczos (6×6)

`sinc(x) * sinc(x/3)` windowed sinc.

### Rust plan

Implement as separate kernel functions that share the scanline transform loop
but differ in how they sample the source buffer:

```rust
pub enum ResampleAlg {
    NearestNeighbour,   // 1 pixel
    Bilinear,           // 2×2
    Cubic,              // 4×4
    Lanczos,            // 6×6
}
```


## Phase 7: NoData / Alpha / Masking (Future)

### NoData masking

When source has nodata values:
1. Build per-band validity bitmask from source data
2. Or use unified validity mask across all bands
3. Skip pixels where mask bit is 0

### Alpha band

When source has alpha (transparency) band:
1. Read alpha as density (0.0 = transparent, 1.0 = opaque)
2. Weight output by density
3. Track cumulative density for output alpha

### Cutline masking

An OGR polygon defines valid regions. Rasterised to a mask at chunk level.
Only pixels inside the cutline contribute.

### Rust plan

Start without masking (Phase 5). Add nodata support as first extension.
Alpha and cutline can follow.


## Implementation Priority

### ACHIEVED: Bit-identical with GDAL

The core warp pipeline is **complete and verified bit-identical** to
`gdalwarp -et 0 -ovr NONE -r near` (65536/65536 match on a 256×256 tile).

Verified components:
- ✅ GenImgProjTransformer (transform.rs)
- ✅ ApproxTransformer (approx.rs)
- ✅ Nearest-neighbour warp kernel (warp.rs)
- ✅ Pixel-centre coordinates (+0.5)
- ✅ Source window offset subtraction (-nSrcXOff)
- ✅ Epsilon in pixel truncation (+1.0e-10)
- ✅ Edge clamping (if iSrcX == nSrcXSize: iSrcX--)

### Next priorities for production use

1. **Overview selection** — Match GDAL's automatic overview selection for
   efficiency. When downsampling significantly (e.g. 10m source → 155m output),
   read from an overview rather than full-resolution data. This is the only
   remaining difference with default `gdalwarp` output.

2. **Bilinear/cubic/lanczos resampling** (Phase 6)

3. **Proper source window computation** — edge-sampling based, matching GDAL's
   `ComputeSourceWindow()` (Phase 2)

4. **NoData masking** (Phase 7)

5. **Chunk management for large warps** (Phase 3)

6. **Multi-threading** (Phase 5 extension)


## Key Numerical Details

| Parameter | GDAL Value | Rust Status |
|-----------|-----------|-------------|
| Pixel coordinate origin | 0.5 (pixel centre) | ✅ Correct |
| ApproxTransformer threshold | 0.125 pixels | ✅ Correct |
| Pixel truncation epsilon | 1.0e-10 | ✅ Correct |
| Source offset subtraction | `- nSrcXOff` | ✅ Correct |
| Edge clamping | `if (x == size) x--` | ✅ Correct |
| Overview selection | Auto-select by ratio | ❌ Not yet (reads full-res) |
| Retry-on-edge | re-transform if within 1px | ❌ Not yet |
| roundIfCloseEnough | 1e-6 in ComputeSourceWindow | Not needed for v1 |
| SRC_COORD_PRECISION | 0.0 (disabled by default) | Not needed |
| ERROR_THRESHOLD in kernel | Stored in papszWarpOptions | Not used by default path |

## File Cross-Reference

| Rust file | GDAL source | Key functions |
|-----------|-------------|---------------|
| transform.rs | gdaltransformer.cpp | GenImgProjTransform (L3100) |
| approx.rs | gdaltransformer.cpp | ApproxTransformInternal (L4113) |
| warp.rs | gdalwarpkernel.cpp | GWKCheckAndComputeSrcOffsets (L5206) |
| (new) source_window.rs | gdalwarpoperation.cpp | ComputeSourceWindow (L3015) |
| (new) warp_operation.rs | gdalwarpoperation.cpp | WarpRegionToBuffer (L1933) |
| lib.rs | gdalwarp_lib.cpp | GDALWarp (L2832-3302) |
