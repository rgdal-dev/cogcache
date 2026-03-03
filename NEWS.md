# cogcache 0.0.1.9010

## Overview selection

New R function `select_overview()` picks the coarsest overview level whose
resolution is still finer than the destination, keeping the warp kernel
operating at ~1:1 source/destination ratio. At this ratio, bilinear and cubic
are bit-identical to GDAL (65536/65536 match on GEBCO test tile). Lanczos shows
minor differences (max_diff=80, corr=0.999991) due to residual kernel width
scaling at ratios slightly above 1.0.

The overview selection logic matches GDAL's `-ovr AUTO` behaviour: walk the
overview pyramid, select the level where the overview-relative source window is
at least as large as the destination in both dimensions. Uses `OVERVIEW_LEVEL`
open option in gdalraster for clean access to overview bands.

## Antimeridian source-read characterisation

Diagnostic analysis of antimeridian-crossing tiles (Fiji LCC, GEBCO) shows the
source X distribution is cleanly bimodal with a 98% gap. The existing
`compute_source_window` antimeridian heuristic correctly flags these tiles but
reads the entire latitude band (97.9M pixels). A split read of the two actual
source strips requires 1.3M pixels — a 70× reduction. Detection and split-read
strategy documented in `inst/design-docs/warp-plan-efficiency-notes.md`.

## Source window padding

Added +1 extra pixel padding in `compute_source_window` for interpolation
kernels, ensuring bilinear/cubic/lanczos neighbourhoods are fully available at
buffer edges. This matches GDAL's `nExtraSrcPixels` logic.




# cogcache 0.0.1.9009

Cogcache is a stepping stone toward a decomposed, inspectable reimplementation
of GDAL's warp pipeline in Rust, callable from R via extendr. This version
establishes the core warp architecture with validated resampling kernels.

## Warp pipeline

Four Rust modules implement GDAL's warp architecture as independent, composable
components:

* `GenImgProjTransformer`: the three-step pixel→geo→CRS→pixel coordinate transform, 
 matching GDAL's `GDALGenImgProjTransform` (cross-referenced to GDAL commit a3b7b01d3e).
* `ApproxTransformer`: adaptive coordinate interpolation wrapping any transformer, matching 
 GDAL's `GDALApproxTransform` with configurable max_error (default 0.125 source pixels).
* `ComputeSourceWindow`: determines which source pixels are needed for a given destination 
 window, including edge sampling with grid fallback, antimeridian heuristic, and 
 resampling kernel padding.
* Resampling kernels: nearest-neighbour, bilinear, cubic, and Lanczos, implemented as 
 separable filters with scanline-oriented processing.

## Resampling validation

Nearest-neighbour is bit-identical to GDAL across all test cases (Sentinel-2,
GEBCO, IBCSO).

Bilinear, cubic, and Lanczos are bit-identical to GDAL at source/destination
ratios near 1:1. At higher downsampling ratios (e.g. 4.5:1), results diverge
because GDAL scales the interpolation kernel width by the local downsampling
factor as an antialiasing measure. Our kernels use fixed filter widths. This is
a known, understood difference — not a bug. Overview selection to keep the
kernel operating at ~1:1 is the next development target and belongs in cogcache
as part of the resampling story.

## Grid arithmetic

Grid operations (inverse geotransform, extent/dimension conversions) use the
vaster-rs crate from hypertidy/vaster-rs, replacing an earlier inline `grid.rs`
module. This is a light dependency that separates logic that is calculated mostly
inline in the GDAL implementation. 

## COG tile access

Direct HTTP fetch and DEFLATE decoding of COG tiles via `rust_fetch_decode_tile`
and `rust_decode_tile`, bypassing GDAL's I/O layer. Currently assumes UInt16
data type. This is a prototype convenience for development — production I/O
strategy (direct-decode vs GDAL-via-vsicurl) is a future design decision.

## Error handling

No panics cross the Rust→R boundary. All PROJ and transform failures return
`Result` types that map to R errors via extendr. The ApproxTransformer
propagates inner transformer failures without masking them.

## R interface

Exported functions: `rust_warp_resample` (main entry point with selectable
resampling algorithm and ApproxTransformer), `rust_compute_source_window`,
`rust_gen_img_proj_transform`, `rust_warp_scanline`, `rust_warp_approx`
(nearest-neighbour with approx), `rust_decode_tile`, `rust_fetch_decode_tile`.
Legacy functions `rust_warp_map` and `rust_apply_warp` remain from the
pre-scanline architecture and will be removed at the vwarp rename.

## Documentation

Design rationale and architecture decisions are in
`inst/design-docs/VWARP_VISION.md` and `src/rust/src/DESIGN_NOTES.md`. Test and
diagnostic and demonstration scripts are in `inst/[test|demo|diag]/-scripts/`.
