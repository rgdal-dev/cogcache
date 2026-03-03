# cogcache 0.0.1.9011

## Chunk list subdivision (antimeridian efficiency)

New Rust function `collect_chunk_list` implements GDAL's
`CollectChunkListInternal` recursive destination subdivision, with a key
improvement: **discontinuity-aware splitting**.

### The problem

GDAL's `ComputeSourceWindow` antimeridian heuristic correctly snaps to full
source width when the backward probe spans >90% of the source raster. This
prevents missed pixels but reads the entire latitude band (97.9M pixels for a
256x256 GEBCO tile over Fiji). GDAL mitigates this through its memory budget
(`dfWarpMemoryLimit`) which forces subdivision of large source reads. Our Rust
implementation had no memory budget path, so the fill-ratio was the only
subdivision trigger — and it was masked by the snap (full-width read looks like
100% fill).

### The fix (Rust)

Three changes to `source_window.rs`:

1. **Forced low fill ratio on antimeridian snap.** When the >90% width heuristic
   fires, force `fill_ratio = 0.01` instead of computing from post-snap
   dimensions. This correctly signals that the read is wasteful.

2. **Discontinuity range detection** (`find_discontinuity_range`). Before blind
   binary subdivision, scan multiple rows through the destination chunk,
   transform each to source coordinates, and find the largest source-X gap. If
   it exceeds 50% of source width, return the column range where the
   discontinuity occurs across all sampled rows. This is a per-pixel scanline
   probe — cheap (256 transforms for a 256-wide tile) and precise.

3. **Three-way split.** When a discontinuity range is found, split the
   destination into: left chunk (entirely east side, compact read), narrow middle
   strip (straddles the discontinuity, accepts full-width read), right chunk
   (entirely west side, compact read). The middle strip is emitted without
   further recursion.

New exported R function: `rust_collect_chunk_list`.

### The fix (R)

New R function `plan_warp_reads` (in `R/plan.R`) wraps `rust_collect_chunk_list`
and post-processes the middle strip:

- `detect_split_strips` transforms the strip's destination pixels to source
  coordinates, finds the bimodal gap, and returns two compact source read
  windows (one per side of the antimeridian).
- The warp pipeline executes split-read chunks by warping each strip
  independently and merging non-nodata values.

### Results (Fiji LCC tile [5,3], GEBCO)

```
Chunk 1: dst=(0,0,80,256)   src=(86021,24871,379,1134)  lon=[178.4,180.0]
Chunk 2: dst=(80,0,5,256)   split=TRUE
  read[1]: src=(0,24877,19,1120)      lon=[-180.0,-179.9]
  read[2]: src=(86384,24877,16,1120)  lon=[179.9,180.0]
Chunk 3: dst=(85,0,171,256) src=(0,24856,808,1147)      lon=[-180.0,-176.6]

Total: 1,395,762 pixels (1.4% of full-width 99,273,600)
```

70x reduction in source I/O. Adjacent non-straddling tiles are unaffected
(single compact chunk, fill_ratio ~1.0).

### GDAL source cross-references

| Our code | GDAL equivalent | gdalwarpoperation.cpp |
|---|---|---|
| `collect_chunk_list` | `CollectChunkListInternal` | lines 1456-1624 |
| `find_discontinuity_range` | (no equivalent — GDAL uses memory budget) | — |
| `build_source_grid` | `ComputeSourceWindowStartingFromSource` setup | lines 2656-2703 |
| `refine_from_source` | `ComputeSourceWindowStartingFromSource` per-call | lines 2720-2747 |
| forced fill_ratio | `dfWarpMemoryLimit` check | lines 1528-1531 |

### New files

- `src/rust/src/source_window.rs`: `collect_chunk_list`, `find_discontinuity_range`,
  `build_source_grid`, `refine_from_source`, `ChunkPlan`, `SourceGrid`
- `R/plan.R`: `plan_warp_reads`, `detect_split_strips`
- `inst/test-scripts/test_chunk_list.R`, `inst/test-scripts/test_plan.R`

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
