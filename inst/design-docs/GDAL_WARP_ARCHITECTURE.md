# GDAL Warp Pipeline — Rust Implementation Status

Based on GDAL source commit `a3b7b01d3e` (3.13.0dev, 2026-02-23).

## The call chain

```
gdalwarp
  └─ GDALWarp()                              gdalwarp_lib.cpp
      ├─ GDALCreateGenImgProjTransformer2()  → transform.rs     ✅ DONE
      ├─ GDALCreateApproxTransformer()       → approx.rs        ✅ DONE
      └─ GDALWarpOperation
          ├─ Initialize()                    (pole detection: TODO)
          ├─ CollectChunkListInternal()
          │   └─ ComputeSourceWindow()       → source_window.rs ✅ DONE
          └─ WarpRegionToBuffer()
              ├─ GDALRasterIO (read)         (caller's job in our arch)
              └─ GWKNearestThread()          → warp.rs          ✅ DONE
```

## Completed components

### Transformer (transform.rs)

Three-step composition: dst pixel → dst geo → src geo → src pixel.
`GenImgProjTransformer` wraps a PROJ transform between two geotransforms.
Implements the `Transformer` trait (batch `transform(inverse, &mut x, &mut y)`).

**Validated**: bit-identical to GDAL's `GDALGenImgProjTransform`.

### ApproxTransformer (approx.rs)

Recursive bisection with linear interpolation, matching GDAL 3.13.0dev's
`GDALApproxTransformInternal`. Checks error at midpoint and quarter-points.
Default threshold: 0.125 source pixels.

**Validated**: bit-identical to GDAL across Lambert Conformal Conic.

### ComputeSourceWindow (source_window.rs)

Pure function on `&impl Transformer`. Algorithm:

1. Corner check: if any dest corner fails → ALL mode (every edge pixel)
2. Edge sampling: `DEFAULT_STEP_COUNT (21) × 4` perimeter points
3. Grid fallback: if any edge point fails → `(step_count+2)²` interior grid
4. Snap-if-close (1e-6), resampling padding, scale adjustment
5. Antimeridian heuristic: >90% source width → full width
6. Fill ratio for chunk-splitting decisions

Returns `SourceWindow { xoff, yoff, xsize, ysize, fill_ratio, n_failed, n_samples }`.

**Validated**: correct on all 48 Fiji LCC tiles including antimeridian cases.

### Warp kernel (warp.rs)

Per-scanline nearest-neighbour: build coords at pixel centres (+0.5),
batch transform, truncate with epsilon (+1e-10), subtract source offset,
bounds check, sample.

**Validated**: 65536/65536 match vs `gdalwarp -et 0 -ovr NONE -r near`.

## Key numerical details

| Parameter | Value | Status |
|-----------|-------|--------|
| Pixel coordinate origin | +0.5 (centre) | ✅ |
| ApproxTransformer threshold | 0.125 px | ✅ |
| Pixel truncation epsilon | +1e-10 | ✅ |
| Source offset subtraction | `- nSrcXOff` | ✅ |
| Edge clamping | `if (x == size) x--` | ✅ |
| roundIfCloseEnough | 1e-6 in source_window | ✅ |
| Overview selection | auto by ratio | ❌ TODO |
| Retry-on-edge | re-transform if within 1px | ❌ TODO |
| Pole detection | aDstXYSpecialPoints | ❌ TODO |
| Bilinear/cubic/lanczos | filter kernel | ❌ TODO |

## File cross-reference

| Rust file | GDAL source | Key functions |
|-----------|-------------|---------------|
| transform.rs | gdaltransformer.cpp ~L2800 | GDALGenImgProjTransform |
| approx.rs | gdaltransformer.cpp ~L4113 | GDALApproxTransformInternal |
| source_window.rs | gdalwarpoperation.cpp ~L2754,3015 | ComputeSourceWindow(TransformPoints) |
| warp.rs | gdalwarpkernel.cpp ~L5206,5510 | GWKCheckAndComputeSrcOffsets, GWKNearestThread |
| lib.rs | gdalwarp_lib.cpp ~L2832-3302 | GDALWarp orchestration |
