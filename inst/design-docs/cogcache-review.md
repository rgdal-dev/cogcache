# cogcache Code Review — March 2026

Review of the cogcache Rust warp pipeline (v0.0.1.9008), with attention to idiomatic usage, dependency choices, the known bilinear edge bug, and fit within the hypertidy ecosystem.

## Overall assessment

This is genuinely impressive work. The four Rust modules (transform, approx, source_window, warp) map cleanly onto GDAL's warp architecture, the GDAL source cross-referencing is meticulous, and the bit-identical nearest-neighbour results against production `gdalwarp` are a strong foundation. The VWARP_VISION document is one of the clearest articulations I've seen of why "reproject(data, crs) is a mistake" and what a better architecture looks like.

The code is ready for the next phase (crate extraction and vwarp rename). What follows is a catalogue of things to fix, think about, or tighten up before that happens.


## 1. The bilinear "bug" — kernel scaling for downsampled warps

The conversation transcript from Part 1 identified a max_diff of 1672 at pixel [1,157] when comparing our bilinear output against GDAL's on a GEBCO tile. Through systematic elimination of every other hypothesis (kernel formula, coordinate transform, geotransform, nodata handling, overview selection, chunking, source data), the root cause was identified:

**The source/destination ratio was 4.5:1.** At this downsampling ratio, GDAL scales the bilinear kernel from 2×2 to ~9×9, acting as an antialiasing filter. Our kernel always uses the textbook 2×2 bilinear regardless of scale. At a steep bathymetric gradient (a -6306m pit surrounded by -3400m shelf), the 2×2 kernel samples the pit bottom while GDAL's wider kernel averages in the surrounding shelf.

**Verification**: at ~1:1 ratio (500m output resolution matching GEBCO's ~430m source), our bilinear produces **65536/65536 bit-identical** results to GDAL.

The kernel formula, coordinate transforms, and all infrastructure are correct. The difference is a missing feature, not a bug. See section 9 of `src/rust/src/DESIGN_NOTES.md` for the full investigation trail, verification table, and fix options.

For the vwarp architecture, the most practical fix is **overview selection in the plan step**: read the source at the overview level closest to the destination resolution, then warp at ~1:1. The plan already knows the source/destination ratio. Kernel scaling (matching GDAL's `dfXScale`/`dfYScale` logic) can be added later for completeness.


## 2. Idiomatic Rust issues

### Type narrowing

The `warp_resample` function takes `src_pixels: &[i32]` for what is actually Int16 (GEBCO) or UInt16 (Sentinel-2) data. This works because the extendr bridge maps R integer to i32, but it means the inner kernel functions do `src[off] as f64` conversions from i32 that could be from i16. For correctness this is fine (i32 → f64 is lossless), but it's worth a comment noting the type widening is intentional and that future multi-type support would use generics here.

### The `SampleFn` type alias

```rust
type SampleFn = fn(&[i32], usize, usize, f64, f64, i32) -> Option<f64>;
```

This is fine for now, but when you extract to a standalone crate, consider making this generic over pixel type. A trait-based approach would let you handle f32/f64/i16/u16 without duplicating kernel code:

```rust
trait Sample: Copy + Into<f64> + PartialEq { ... }
```

### Casting patterns

`(buf_x - 0.5).floor() as i64` — this is correct but Rust's `as` cast on float-to-int is saturating in recent editions (no UB), which is good. Worth a comment noting this matches C's `static_cast<int>(floor(...))` semantics.

`src_col_off as usize` in `warp_resample` — the caller passes these as `usize` already, but `lib.rs` receives them as `i32` from R and converts with `as usize`. If `src_col_off` were ever negative (it shouldn't be), this would silently wrap. Consider adding a bounds check in the extendr wrapper:

```rust
if src_col_off < 0 || src_row_off < 0 {
    return Err(Error::Other("source offsets must be non-negative".to_string()));
}
```

### Error handling in transform.rs

The PROJ `convert()` call on line 143 returns `Err` for failed transforms, but the error is discarded. In batch transforms (many points), individual point failures are expected (pole singularities, out-of-domain), so this is correct. But consider using `match` with a more specific pattern:

```rust
Err(proj::ProjError::Conversion(_)) => { success[i] = false; }
Err(e) => { return Err(...); }  // unexpected errors should propagate
```

This distinguishes "this point can't be transformed" (expected) from "PROJ is misconfigured" (a real error).

### The `Vec<bool>` success flags

Every `transform()` call allocates a new `Vec<bool>` for success flags. For scanline processing this is fine (256-512 bools per scanline), but for the approx transformer's recursion, the entry point allocates one master `Vec<bool>` and the internal recursion writes into slices of it — good. The trait signature forces a new allocation per call though. If performance matters later, consider an alternative trait method that takes `&mut [bool]`.

### Module visibility

`transform.rs` and `approx.rs` have `pub struct` and `pub fn` but the modules themselves are `pub mod` in lib.rs. This is correct for the extendr use case (lib.rs needs access). When extracting to standalone crates, the public API surface should be tightened — `GenImgProjTransformer` fields should not be public (they aren't currently, good), and `ApproxTransformer::approx_entry`/`approx_internal` should remain private.


## 3. Dependency review

### Current Rust dependencies

| Crate | Version | Assessment |
|-------|---------|------------|
| `extendr-api` | 0.7 | Current stable for R bindings. Fine. |
| `flate2` | 1 | Standard DEFLATE. Fine. |
| `ureq` | 2 | Blocking HTTP client. Fine for single-tile fetch. For parallel tile fetching (vwarp vision item 11), consider `reqwest` with async, or `object_store` for S3/GCS. |
| `proj` | 0.28 + bundled_proj | Bundles libproj. The key requirement is that PROJ is full-featured with tight access to proj.db — as long as that holds, bundled is fine. The `.onLoad` fix (now applied) handles proj.db discovery across platforms. |
| `vaster` (git) | hypertidy/vaster-rs | Only used for `inv_geotransform`. Correct — vaster is the foundation. Git dependency is fine for development; publish to crates.io before the vwarp rename for reproducibility. |

### No GDAL Rust dependency

The architecture is clear: GDAL (via gdalraster, vapour, or direct CLI) handles I/O and serves as the validation oracle. The Rust code handles planning, coordinate transformation, and warp kernels. Using the `gdal` Rust crate internally would undermine the independence that makes this project interesting — the whole point is that the warp logic is decomposed and inspectable, not tunnelled through GDAL's monolith. GDAL is used to wire up workflows from R, not as an internal dependency.

Any reference to the GDAL C++ source for algorithm comparison should use recent trunk (currently `a3b7b01d3e`, 3.13.0dev). The code isn't pinned to this version but the cross-references should track the latest.

### R-side dependencies

The DESCRIPTION has no R package dependencies — all the R-side code is extendr-generated wrappers. The test scripts use `gdalraster` and `vaster` (R), which is correct. For CRAN, DESCRIPTION would need `Suggests: gdalraster, vaster` for examples/tests.


## 4. Hypertidy ecosystem fit

### vaster overlap

The Rust code uses `vaster::inv_geotransform` — that's it. The R-side test scripts use `vaster::extent_dim_to_gt()`. This is the right level of coupling. The grid arithmetic lives in vaster, the warp logic lives here.

One thing to check: `vaster-rs` and the R `vaster` package should agree on geotransform conventions. The R vaster uses `extent_dim_to_gt()` to produce GDAL-style [xmin, xres, 0, ymax, 0, -yres] geotransforms. The Rust `inv_geotransform` should invert these correctly. The test suite validates this implicitly (bit-identical results), but an explicit unit test of `inv_geotransform(extent_dim_to_gt(...))` → identity round-trip would be good documentation.

### grout relationship

The VWARP_VISION document mentions grout for "tile index ↔ pixel coords" but the current code doesn't use grout at all — the tile indexing in the demo scripts is done manually:

```r
dst_col_off <- (tx - 1L) * tile_size
tile_xmin <- out_xmin + dst_col_off * pixel_res
```

This is grout's job. When the vwarp R package is structured, the tile iteration should use grout's tile index, not manual arithmetic. This would also handle the common pattern of "give me the extent and offset for tile [tx, ty] in this grid" which appears in every demo script.

### vapour / gdalraster usage

The test scripts use `gdalraster::GDALRaster` for reading source data and `gdalraster::warp()` as the reference oracle. This is correct. The production pipeline should continue using GDAL (via gdalraster or vapour) for I/O and use the Rust code for planning and kernel execution.

However, there's a design tension: `rust_fetch_decode_tile` in lib.rs does its own HTTP range request + DEFLATE decode, bypassing GDAL entirely. This is fine for COGs with known DEFLATE+predictor2 encoding, but it's a parallel I/O path that doesn't benefit from GDAL's block cache or support other compression methods (LZW, JPEG, ZSTD, etc.). For the vwarp rename, consider whether this direct-decode path is a prototype convenience or a long-term architecture choice. If the latter, it needs to handle more codecs. If the former, the production path should read via GDAL (passing byte ranges to /vsicurl/ or /vsimem/).

### ndr relationship

ndr provides xarray-like lazy arrays. The vwarp warp plan ("for target chunk X, read source chunks Y and Z, warp them") maps directly onto ndr's lazy computation model. The plan is the recipe; ndr is the executor. This is complementary but not yet integrated. The key interface point will be: "given a vwarp plan and an ndr lazy array, materialise the warped output."

### tissot overlap

Tissot computes Tissot indicatrices (map projection distortion ellipses). It uses the same grid → geo → CRS pipeline as GenImgProjTransformer, but for a different purpose (distortion analysis vs. pixel reprojection). The `Transformer` trait in transform.rs could be shared with tissot's internals if both are extracted to a common crate. Consider whether `vwarp-transform` should be the shared foundation.


## 5. Code quality observations

### The `.onLoad` PROJ_DATA handling

```r
.onLoad <- function(libname, pkgname) {
  PROJ_DATA <- Sys.getenv("PROJ_DATA")
  proj_data <- "/usr/share/proj"
  Sys.setenv(PROJ_DATA = proj_data)
```

This unconditionally overwrites PROJ_DATA with a hardcoded path. On macOS (Homebrew), PROJ data lives in `/opt/homebrew/share/proj`. On Conda environments, it's somewhere else entirely. Use `sf::sf_proj_search_paths()` or `terra::gdal_proj_search_paths()` if available, otherwise detect the system path. At minimum:

```r
.onLoad <- function(libname, pkgname) {
  if (Sys.getenv("PROJ_DATA") == "") {
    candidates <- c("/usr/share/proj", "/opt/homebrew/share/proj",
                     "/usr/local/share/proj")
    for (p in candidates) {
      if (dir.exists(p)) { Sys.setenv(PROJ_DATA = p); break }
    }
  }
}
```

### The `decode_deflate_pixels` UInt16 handling

```rust
let mut pixels: Vec<i32> = decompressed
    .chunks_exact(2)
    .map(|b| u16::from_le_bytes([b[0], b[1]]) as i32)
    .collect();
```

This reads little-endian UInt16 and widens to i32. For GEBCO (Int16 with nodata=-32768), this produces values 0-65535 for negative depths because it reads as u16, not i16. The predictor-2 undo then produces correct relative differences, but the absolute values are wrong (they're in 0..65535 instead of -32768..32767).

Wait — actually, `u16::from_le_bytes` followed by `as i32` gives 0..65535. For Int16 data like GEBCO, you want `i16::from_le_bytes` followed by `as i32` to get -32768..32767. The predictor-2 modular arithmetic (`% 65536`) happens to produce the correct relative values, but the absolute pixel values returned to R would be unsigned-reinterpreted. If this is used with the warp kernel's nodata check (`val != nodata` where nodata is -32768 passed from R as i32), it would never match because the decoded values are 0..65535.

This suggests `rust_fetch_decode_tile` + `rust_decode_tile` are only used for UInt16 data (Sentinel-2), not Int16 (GEBCO). For the GEBCO demos, `gdalraster::read()` is used instead. Document this limitation or fix the decode to accept a data type parameter.

### Test coverage

The Rust unit tests cover kernel weight functions and flat-field interpolation. Missing:

- Bilinear at known fractional positions (verify exact interpolation value, not just "close to flat")
- Cubic and Lanczos at pixel centres (should exactly reproduce the centre pixel value)
- Edge-case: kernel window partially outside buffer (the bug case)
- ApproxTransformer: verify that `max_error=0` produces identical results to direct GenImgProjTransformer
- source_window: verify padding is sufficient for bilinear/cubic/lanczos

The R test scripts in `inst/test-scripts/` cover these cases empirically against GDAL output, which is great for integration testing but doesn't help when running `cargo test` without GDAL/R.


## 6. Naming and API surface

### The cogcache → vwarp rename

The DESCRIPTION still says "COG Tile Cache with Rust Decode." The code has evolved far beyond that. The rename is overdue. For the crate workspace:

- `vwarp-transform` → Transformer trait + GenImgProjTransformer + ApproxTransformer
- `vwarp-plan` → ComputeSourceWindow + chunk planning
- `vwarp-kernel` → warp kernels (nearest, bilinear, cubic, lanczos)
- `vwarp` → R/Python bindings, orchestration

Keep `cogcache` as a git history artifact. Don't try to preserve the package name.

### Legacy function pruning

`rust_warp_map` and `rust_apply_warp` are the pre-scanline architecture. They allocate full coordinate grids (O(n²) memory) instead of working scanline-by-scanline. They should be deprecated before the rename and removed after.

`rust_warp_scanline` does exact per-pixel PROJ (no approx). It's useful for validation but shouldn't be in the main API. Move to a `_debug` or `_exact` suffix.

`rust_warp_approx` is now a subset of `rust_warp_resample` (it's `rust_warp_resample` with `resample = "near"`). Consider removing it and having `rust_warp_resample` with a default `resample = "near"`.


## 7. Summary of actions

**Done in this review:**

1. ✅ Identified the bilinear discrepancy as kernel scaling for downsampled warps (4.5:1 ratio), not a bug — verified bit-identical at 1:1 ratio
2. ✅ Fixed `.onLoad` PROJ_DATA to detect system path instead of hardcoding `/usr/share/proj`
3. ✅ Added +1 extra padding in source_window.rs for interpolation kernel edge safety (still useful, just not the cause of the big diffs)
4. ✅ Documented the full investigation trail in BILINEAR_SCALING_NOTES.md

**At vwarp rename:**

5. Remove legacy `rust_warp_map`, `rust_apply_warp` from Rust, R wrappers, NAMESPACE, man pages
6. Add input validation for negative offsets in extendr wrappers
7. Publish vaster-rs to crates.io
8. Hard rename cogcache → vwarp
9. Document that `rust_decode_tile` assumes UInt16 (or fix for Int16)

**Short-term (kernel completeness):**

10. Implement overview selection in the plan step — read source at the overview level closest to destination resolution, then warp at ~1:1
11. Consider implementing kernel scaling (GDAL's `dfXScale`/`dfYScale` logic) for downsampled warps

**Medium-term (architecture):**

12. Integrate grout for tile indexing instead of manual arithmetic in demo scripts
13. Define the ndr → vwarp plan interface
14. Decide on direct-decode vs GDAL-via-vsicurl for production I/O
15. Consider shared `Transformer` crate with tissot
