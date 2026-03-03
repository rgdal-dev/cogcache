## Bilinear kernel scaling for downsampled warps

**Root cause found: 2026-03-03**

### The problem

When warping GEBCO bathymetry (86400×43200) to a 256×256 tile at 2km
resolution with bilinear resampling, our output differed from GDAL's
by up to 1669 depth-metres at steep gradients. Correlation was 0.9996
instead of 1.0.

### What we ruled out

Through systematic elimination:

- **Kernel formula**: manual bilinear in R exactly reproduces our Rust output ✅
- **Coordinates**: Rust PROJ = sf PROJ = reproj, zero difference at 100 sampled pixels ✅
- **Geotransform**: our tile_gt = GDAL output geotransform, zero difference ✅
- **NoData path**: GDAL with/without nodata produces identical output ✅
- **ApproxTransformer**: et=0 (exact) gives same result as et=0.125 ✅
- **Source data**: GDAL from synthetic small-tif = GDAL from COG, zero difference ✅
- **Overview selection**: -ovr NONE, DISABLE_READDIR_ON_OPEN, all same result ✅
- **Chunking/memory**: -wm 2000 (huge memory), SOURCE_EXTRA=100, all same result ✅
- **Nearest-neighbour**: 65536/65536 match between us and GDAL ✅

### The root cause

The source/destination ratio was 4.6×4.5 — each destination pixel covers
~4.5 source pixels in each direction. At this ratio, **GDAL scales the
bilinear kernel** to cover the appropriate source area (~9×9 instead of
2×2). This is the antialiasing behaviour in `gdalwarpkernel.cpp`:

```
dfXScale = (double)poWK->nDstXSize / (dfSrcXMax - dfSrcXMin);
dfYScale = (double)poWK->nDstYSize / (dfSrcYMax - dfSrcYMin);
if (dfXScale < 1.0) nXRadius = (int)ceil(nRadius / dfXScale);
if (dfYScale < 1.0) nYRadius = (int)ceil(nRadius / dfYScale);
```

When dfXScale < 1 (downsampling), the kernel radius is expanded by 1/scale.
For bilinear (radius=1) at 4.5:1, this means radius ≈ 4.5 → a ~9×9 kernel.
The weights within this expanded kernel are computed using the standard
bilinear triangle function, but stretched to cover the larger area.

### Verification

At ~1:1 source/destination ratio (500m output resolution vs GEBCO's ~430m):
**65536/65536 bit-identical match** between our Rust bilinear and GDAL's.

This confirms:
- Our bilinear kernel formula is correct
- Our coordinate transform pipeline is correct
- The only difference is kernel scaling for downsampled warps

### Impact

Our bilinear/cubic/lanczos kernels are correct for:
- **1:1 ratio** (matching source resolution): bit-identical to GDAL ✅
- **Upsampling** (finer output than source): bit-identical to GDAL ✅
- **Downsampling** (coarser output than source): produces aliased output ✗

For the vwarp use case (warp-to-tiles for web maps or analysis), the target
resolution is typically chosen to match or exceed the source resolution.
Downsampled warps are less common — you'd normally read from an overview
level instead.

### Fix options

1. **Document the limitation**: our bilinear is point-sampled, not
   area-integrated. For downsampled warps, use overviews or pre-filter.
   This is a valid choice for a lightweight warp library.

2. **Implement kernel scaling** (match GDAL): scale the kernel radius by
   1/ratio when downsampling. The `source_window.rs` already computes
   the scale factor for padding; the kernel just needs to use it.
   
   In `warp.rs`, `warp_interpolated` would need to receive `dx_scale`
   and `dy_scale` and pass them to the sample functions. The bilinear
   kernel would switch from 2×2 to (2/scale)×(2/scale) with stretched
   weights.

3. **Pre-downsample the source**: read the source at the overview level
   closest to the destination resolution, then warp at ~1:1. This is
   what GDAL's overview selection (`-ovr AUTO`) does. It's simpler than
   kernel scaling and produces good results.

Option 3 is the most practical for the vwarp architecture — the plan
step already knows the source/destination ratio and can select the
appropriate overview level. Kernel scaling (option 2) is worth adding
later for completeness.

### GDAL source references

- `gdalwarpkernel.cpp` ~L3800: `dfXScale`/`dfYScale` computation
- `gdalwarpkernel.cpp` ~L2460: `GWKResampleNoMasksOrDstDensityOnlyThread`
  general case with scaled kernels
- `gdalwarpkernel.cpp` ~L5510: `GWKNearestThread` — no scaling (always 1 pixel)
- `gdalwarp_lib.cpp` ~L1455: `-ovr` handling, overview selection logic
