# Warp plan efficiency notes

## The problem

Coordinate transforms between arbitrary CRS pairs can produce source windows much larger than what the destination actually needs. The antimeridian is the most common singularity for geographic sources, but the same class of problem arises at the poles (polar stereographic → any equatorial CRS), at projection zone boundaries (UTM zone edges), and with any transform where the Jacobian is singular or the domain wraps.

## The antimeridian case (characterised)

For Fiji LCC (lon_0=178) warping from GEBCO (global, -180 to +180), tile [5,3] straddles the antimeridian. Source X coordinates wrap from ~86029 (178.5°E) to ~800 (-176.7°E = 183.3°E). The bounding box spans 86400 columns — the entire latitude band.

`compute_source_window` catches this with the >90% width heuristic and snaps to full width. This is correct (no pixels are missed) but reads 97.9M pixels for a 256×256 tile.

The source X distribution is cleanly bimodal with a 98% gap (columns 864 to 85536 are empty). A split read of the two actual strips requires 1.3M pixels — 1.4% of the full-width read, a 70× reduction.

Adjacent tiles [4,3] (173.5°–178.5°E) and [6,3] (-176.8°–-171.8°) are entirely on one side and produce compact source windows (~1200 columns each). Only the straddling tile is affected.

See also: opendatacube/odc-geo#208 and opendatacube/odc-stac#172 for the same problem in the Python ecosystem, where the round-trip UTM → WGS84 → Web Mercator flips coordinates to the wrong side of the antimeridian, producing pixel indices offset by one full world width.

## What GDAL does

GDAL addresses singularities through heuristics accumulated over 20+ years:

* **Antimeridian snap** (`ComputeSourceWindow` ~L3325): if the source window covers >90% of source width, snap to full width. Simple, reliable, wasteful.
* **Recursive destination subdivision** (`CollectChunkListInternal`): when the source window exceeds the warp memory budget, subdivide the destination recursively. Each sub-window gets its own `ComputeSourceWindow`. Eventually sub-windows are small enough to not straddle the singularity. This implicitly splits the antimeridian case without special-casing it.
* **Forward probing** (`ComputeSourceWindowStartingFromSource`): transform the source extent forward to destination coords to check if the source actually contributes. Catches cases where the backward transform wraps but the forward transform shows the source is on one side only.
* **Pole detection** (`aDstXYSpecialPoints`): special points for polar singularities.

A full catalogue of GDAL's singularity-handling strategies and which CRS pairs trigger them would be a valuable analysis — there are likely more heuristics embedded in specific transform paths.

## Position for cogcache/vwarp

**The kernel stays simple.** It warps a source buffer to a destination buffer. It does not know about CRS singularities. This is by design — the kernel's contract is: given correct source data and correct coordinates, produce correct output.

**The planning layer handles singularities.** Detection of inefficient source reads, split-read strategies, overview selection, and transform caching all live in R (or Python). This code is inspectable, testable in isolation, and doesn't need to be fast — it runs once per tile, not once per pixel.

**Detection should be empirical, not hardcoded.** Rather than maintaining a list of "problem CRS pairs," detect singularities by examining the actual transformed coordinates: bimodal source distributions, source windows much larger than destination extent × local scale, large gaps in the source coordinate field. These signals are CRS-agnostic.

**Common transforms deserve cached plans.** For production pipelines (Sentinel-2 UTM → EPSG:3857, Landsat → WGS84, GEBCO → polar stereo), the transform plan — overview level, source strips, split strategy — is the same for every tile in the same grid position. MGRS tiles share transform geometry across all 60 UTM zones (modulo a longitude offset), reducing the unique plan space to ~20 zone templates × a handful of target CRSes. Precomputing and caching these plans is a natural optimisation.

**Adaptive techniques can refine detection.** The boundary-sampling approach (as in `compute_source_window`) works for most cases but can miss interior singularities in highly nonlinear transforms. Progressively more robust techniques:

1. **Gap detection** (implemented): sort source X, find largest gap, split if >50% of source width. Handles the antimeridian. Cheap.
2. **Adaptive boundary densification** (cf. hypertidy/bigcurve, D3/topojson adaptive resampling): bifurcate boundary segments until the transformed curve is well-represented. Handles nonlinear distortion near singularities without requiring a full grid. The bifurcation criterion can use the segment length in source coords vs the chord, similar to how D3 densifies projected linestrings.
3. **Local distortion analysis** (cf. tissot, proj_factors): use the Jacobian of the transform to identify where distortion is extreme and where adaptive densification is needed. High scale factors or angular distortion flag regions where the boundary sampling is unreliable. This connects directly to the tissot package's distortion indicators.
4. **Field analysis**: treat the full grid of transformed coordinates as a 2D field. Discontinuities appear as large gradients relative to local scale. More expensive (requires transforming a grid) but catches interior singularities that boundary methods miss. Relevant for polar projections where the boundary is well-behaved but interior points sweep through all longitudes.

These techniques layer naturally — start with gap detection, escalate to adaptive densification for flagged tiles, use field analysis for pathological projections. The planning layer's job is to try the cheapest method first and escalate only when needed.

## Connection to the wider ecosystem

This is the "grey area" between planning and execution. The planning code is quick to write in R or Python (especially with structured tools like `compute_source_window` and `rust_gen_img_proj_transform` available), but the decisions it makes determine the efficiency of the entire pipeline. Keeping it in the scripting layer means it's easy to experiment with different strategies, profile them against real workloads, and harden the ones that matter.

The cached-plan approach also connects to the broader question of interoperability with tiling systems (MGRS, Mercator tile pyramids, polar grids). Each tiling system defines a finite set of destination geometries, and each source CRS defines the transform characteristics. The cross product is large but highly redundant.
