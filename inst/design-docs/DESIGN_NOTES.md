# Design Notes & Known Edge Cases

**cogcache / vwarp pipeline — living document**

This file records design decisions, hard calls, and known limitations in the
Rust warp implementation. If something here surprises you, that's the point.

Last updated: 2026-03-03

GDAL source reference: commit `a3b7b01d3e` (3.13.0dev, 2026-02-23).

---

## Current status

Every component in the warp pipeline is **implemented and validated
bit-identical to GDAL** (`gdalwarp -r near -ovr NONE`):

| Component | Module | Validation |
|-----------|--------|------------|
| GenImgProjTransformer | `transform.rs` | 65536/65536 pixel-identical |
| ApproxTransformer | `approx.rs` | bit-identical across LCC projection |
| ComputeSourceWindow | `source_window.rs` | matches GDAL, 48/48 tiles correct |
| Nearest-neighbour kernel | `warp.rs` | 65536/65536 vs `gdalwarp -et 0 -ovr NONE` |
| Full pipeline (Fiji LCC) | all modules | 65536/65536 centre tile, 100.0% |

---

## Design decisions

### 1. Truncation vs round for nearest neighbour

Truncation (`as i32`), matching GDAL's `static_cast<int>()`.

### 2. The +0.5 cell-centre convention

Follow GDAL — `transform_scanline()` adds +0.5 so callers don't have to.
Raw `Transformer::transform()` expects the caller to know what they're doing.

### 3. Two PROJ objects instead of one

The `proj` crate doesn't expose a direction flag. Separate instances for
fwd/inv. Doubles init cost, irrelevant for tile serving.

### 4. No rotated geotransform support

Rejected — COGs are always north-up. Full affine inversion trivial to add.

### 5. ComputeSourceWindow as a pure function

Takes `&impl Transformer`, no dataset handle. This makes it a **planning tool**
you can call repeatedly to map a job before reading any pixels. GDAL buries
it inside `GDALWarpOperation` as a protected method.

Skipped: `CHECK_WITH_INVERT_PROJ`, `ComputeSourceWindowStartingFromSource`,
`aDstXYSpecialPoints` pole detection. Addable when needed.

### 6. Edge sampling default with grid fallback

Matches GDAL's cascade. Proven across Fiji LCC (antimeridian + conic limb).
The 90% heuristic fires correctly.

---

## Known limitations

**Antimeridian column-5 blowup.** Source window correctly identifies full-width
reads (~86400 px). GDAL does the same. Fix: antimeridian-aware splitting (two
reads instead of one). Source-window post-processing, not a transform issue.

**No overview selection.** Always reads full-resolution. GDAL auto-selects
overviews when downsampling. Optimisation only, not correctness.

**PROJ axis order.** Relies on `Proj::new_known_crs()` normalisation. Fine for
EPSG codes. Raw WKT2 with explicit lat/lon order may surprise.

**proj.db location.** Set `PROJ_DATA` in `.onLoad()` for production.

---

## GDAL source cross-references

| Our code | GDAL equivalent | File (commit a3b7b01d3e) |
|---|---|---|
| `transform.rs` | `GDALGenImgProjTransform` | gdaltransformer.cpp ~L2800 |
| `approx.rs` | `GDALApproxTransformInternal` | gdaltransformer.cpp ~L4113 |
| `source_window.rs` | `ComputeSourceWindow` | gdalwarpoperation.cpp ~L3015 |
| `source_window.rs` sampling | `ComputeSourceWindowTransformPoints` | gdalwarpoperation.cpp ~L2754 |
| `warp.rs` | `GWKNearestThread` | gdalwarpkernel.cpp ~L5510 |
