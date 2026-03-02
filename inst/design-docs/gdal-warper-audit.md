# GDAL Warper Audit: Anatomy of gdalwarp for a Rust-Native Reimplementation

**Michael Sumner — March 2026 (draft, not for public circulation)**

---

## Executive Summary

The GDAL warper is ~25,000 lines of C++ across 7 files, plus ~5,000 lines of coordinate
transformer code. It looks intimidating but decomposes cleanly into layers with very
different levels of essential complexity. A Rust-native warper built on async-tiff and
object_store could start with ~2,000 lines of core logic and reach practical utility for
the cloud-native COG tile-serving use case surprisingly quickly.

This document maps out what's **foundational** (must be reimplemented), what's **valuable
but separable** (can be layered on later), and what's **legacy/GDAL-coupling** (should be
rethought or dropped entirely).

---

## File Inventory (GDAL trunk `85a6167`, 3.13.0dev — Feb 2026)

| File | Lines | Role |
|---|---|---|
| `alg/gdalwarpkernel.cpp` | 9,283 | Resampling kernels (the actual pixel math) |
| `apps/gdalwarp_lib.cpp` | 6,630 | CLI/library entry point, option parsing, output creation |
| `alg/gdaltransformer.cpp` | 5,151 | Coordinate transformers (GenImgProj, Approx, GCP, RPC, etc.) |
| `alg/gdalwarpoperation.cpp` | 3,388 | Chunking, source window computation, I/O orchestration |
| `frmts/vrt/vrtwarped.cpp` | 2,589 | VRT integration (on-demand warping as a virtual dataset) |
| `alg/gdalwarper.cpp` | 2,505 | Masking functions (nodata, alpha, cutline, density) |
| `alg/gdalwarper.h` | 644 | Data structures: GDALWarpOptions, GDALWarpKernel, GDALWarpOperation |
| `alg/gdalsimplewarp.cpp` | 444 | Simplified warp (rarely used, mostly legacy) |

---

## Layer 1: The Core Algorithm (FOUNDATIONAL)

The actual warping algorithm lives in `GWKGeneralCaseThread()` (lines 5510–5777 of
gdalwarpkernel.cpp) and is remarkably simple:

```
for each output scanline (iDstY):
    build array of destination pixel coordinates (x + 0.5, y + 0.5)
    call pfnTransformer(TRUE, ...) to map dst→src coordinates
    for each output pixel (iDstX):
        check if transform succeeded
        check validity/density masks at source location
        for each band:
            resample source pixels around the transformed coordinate
            write result to destination buffer
        update destination validity/density masks
```

**This is ~270 lines of C++ including all the mask/density handling.** The core loop itself
is maybe 80 lines. In Rust this would be a straightforward generic function.

### What a Rust version needs

1. **Grid specification** — extent, dimensions, CRS, affine transform. This maps directly
   to the `vaster` crate's grid logic. The GDAL warper represents this implicitly through
   dataset handles; a Rust version should make it an explicit, lightweight struct.

2. **Coordinate transform trait** — `fn transform_points(&self, dst_to_src: bool, x: &mut [f64], y: &mut [f64], z: &mut [f64]) -> Vec<bool>`. Wraps PROJ, proj4rs, or hardcoded fast-paths for common cases (WebMercator↔WGS84).

3. **Resampling kernels** — Pure numerical functions. No I/O, no state, no GDAL coupling.
   See Layer 2 below.

4. **Source pixel accessor** — Given fractional source coordinates, read pixel values.
   This is where the async-tiff / cloud-native I/O plugs in. See Layer 3.

---

## Layer 2: Resampling Kernels (FOUNDATIONAL, very portable)

The kernel functions are pure math. GDAL defines these filter functions:

| Method | Kernel radius | Function |
|---|---|---|
| NearestNeighbour | 0 | floor() |
| Bilinear | 1 | 1 - abs(x) |
| Cubic | 2 | Mitchell-Netravali variant |
| CubicSpline | 2 | B-spline |
| Lanczos | 3 | sinc(x) * sinc(x/3) |

Plus the statistical resamplers (Average, Mode, Min, Max, Med, Q1, Q3, Sum, RMS) which
use a different code path (`GWKAverageOrMode`) that computes the set of contributing
source pixels per destination pixel rather than interpolating at a point.

### The specialization explosion

`PerformWarp()` (line 1050) dispatches to **~30 specialized kernel functions** based on
combinations of:
- Data type: Byte, Int16, UInt16, Float32, Float64
- Resampling method: Nearest, Bilinear, Cubic, CubicSpline
- Mask state: "NoMasksOrDstDensityOnly" vs full mask handling

This produces names like `GWKCubicNoMasksOrDstDensityOnlyUShort`. Each is a separate
function with the type and masking logic baked in at the source level for performance.

**In Rust, this entire specialization tree is replaced by generics + monomorphisation.**
A single generic function `warp_kernel<T: Sample, R: Resample, M: MaskPolicy>()` would
compile to equivalent specialized code. This alone eliminates ~4,000 lines of C++.

### SSE2 optimizations

There are SSE2-optimized paths for bilinear/cubic resampling without masks
(`GWKResampleNoMasks_SSE2_T`, line 4894). In Rust, SIMD can be handled via:
- `std::simd` (nightly) or `packed_simd2`
- Compiler autovectorization (often sufficient for these patterns)
- Architecture-specific intrinsics where needed

### Key subtlety: XScale/YScale computation

Before dispatching, `PerformWarp()` computes the local scale ratio between source and
destination pixels (lines 1073–1277). This is done by transforming a grid of unit squares
from destination to source space and measuring how much they distort. The scale determines
the filter kernel width — crucial for correct downsampling.

GDAL computes this **per chunk** and uses a single scale for the whole chunk. The code
itself notes "per-pixel scale computation would probably be best" (line 1078). A Rust
implementation could do per-scanline or per-tile scale computation more easily.

The polar-discontinuity outlier rejection (lines 1226–1243) is directly relevant to
Antarctic work — it discards scale samples that are <10% of the maximum, which handles
the singularity at the pole.

---

## Layer 3: Source Window Computation (FOUNDATIONAL, tricky)

`ComputeSourceWindow()` (gdalwarpoperation.cpp, line 3015) determines which rectangle of
source pixels is needed for a given output region. This is where decades of edge-case
fixes live.

### The basic algorithm

1. Sample points along the **edges** of the destination region (default: 21 per edge)
2. Transform them back to source pixel coordinates
3. Take the bounding box of the successfully transformed points
4. Add padding for the resampling kernel radius

### The edge cases (where the hard-won knowledge lives)

- **Grid sampling** (`SAMPLE_GRID`): When edge-only sampling misses interior extrema.
  Critical for "inside out" projections like WGS84→Polar Stereographic around the pole.
  The code auto-detects this when corner transforms fail (line 3050–3074).

- **Special points** (`aDstXYSpecialPoints`): Pre-computed points in destination space
  where the transform is singular or extreme (e.g., the pole in a polar projection).
  If any fall within the current chunk, grid sampling is forced (lines 3094–3116).

- **Bogus inverse transform detection** (line 3127): Checks for unreasonably large source
  windows that indicate PROJ returned garbage. Retries with `CHECK_WITH_INVERT_PROJ`.
  This catches old PROJ bugs but remains useful as a safety net.

- **Forward transform fallback** (`ComputeSourceWindowStartingFromSource`, line 2620):
  When the inverse transform is unreliable, also samples from source→destination to find
  which source pixels land in the output region. This is the fix for
  [OSGeo/gdal#862](https://github.com/OSGeo/gdal/issues/862).

- **Orthographic projection handling** (line 3053): Detects when destination corners
  project into space (off-Earth) and switches to `bAll=true` (transform every pixel).
  Referenced: [OSGeo/gdal#9056](https://github.com/OSGeo/gdal/issues/9056).

### For a Rust warper

This logic is essential but can be cleaner. The key insight is that for cloud-native
tile serving, the output region is often a **single web mercator tile** — small enough
that the source window computation is cheap, and the edge cases around massive warps
with heterogeneous distortion mostly don't apply.

For a minimal viable implementation: transform the 4 corners + midpoints of edges,
take the bounding box, add kernel padding. Add grid sampling as a configurable fallback
for non-linear transforms.

---

## Layer 4: Chunking & I/O Orchestration (RETHINK for cloud-native)

`GDALWarpOperation` (gdalwarpoperation.cpp) manages:
- Breaking the output into memory-sized chunks (`CollectChunkList`)
- Loading source data for each chunk
- Setting up validity/density masks
- Coordinating I/O and computation (optionally on separate threads via `ChunkAndWarpMulti`)

### What's legacy here

The chunking exists because GDAL assumes you're warping a **whole raster** that's too
large for memory. `dfWarpMemoryLimit` controls how big each chunk can be. The chunks are
computed by subdividing the output until each chunk's memory footprint fits.

In a tile-server context, the "chunk" is the output tile (typically 256×256 or 512×512
pixels). There's no need for the subdivision machinery.

### What's valuable

The I/O pipeline concept — read source data, warp, write output — maps naturally to an
async pipeline:
```rust
async fn warp_tile(output_tile: TileSpec, source: &CogReader) -> Result<RasterBuffer> {
    let src_window = compute_source_window(&output_tile, &transform)?;
    let src_tiles = source.fetch_tiles(src_window).await?;  // async I/O
    let src_buffer = decode_and_assemble(src_tiles).await?;  // async decode
    let dst_buffer = warp_kernel(&output_tile, &src_buffer, &transform, resampling)?;
    Ok(dst_buffer)
}
```

The async source tile fetching is a **major potential advantage** over GDAL's synchronous
block cache. Multiple source tiles can be fetched concurrently from object storage.

---

## Layer 5: Coordinate Transformers (SEPARATE CONCERN)

`gdaltransformer.cpp` (5,151 lines) implements several transformer types:

| Transformer | Purpose | Complexity |
|---|---|---|
| GenImgProjTransformer | The workhorse: pixel↔georef↔CRS pipeline | Medium |
| ApproxTransformer | Caches/interpolates transform results | Medium |
| GCPTransformer | Ground control point based | High |
| RPCTransformer | Rational polynomial coefficients (satellite) | High |
| TPSTransformer | Thin plate spline | High |
| GeoLocTransformer | Geolocation arrays (swath data) | High |

### For a Rust warper

**GenImgProjTransformer** is the only one needed initially. It composes:
1. Source pixel→source georef (affine transform — trivial)
2. Source CRS→destination CRS (PROJ — use `proj` crate or `proj4rs`)
3. Destination georef→destination pixel (inverse affine — trivial)

**ApproxTransformer** is a performance optimization: instead of calling PROJ for every
pixel, it evaluates PROJ at a subset of points and linearly interpolates between them.
The error threshold controls accuracy. This is worth implementing for large warps but
not needed for the tile-serving case where you're transforming 256–512 points per scanline.

GCP, RPC, TPS, and GeoLoc transformers are specialized — defer to later phases.

---

## Layer 6: Masking (IMPORTANT but overengineered)

`gdalwarper.cpp` (2,505 lines) implements masking functions:

| Masker | Purpose |
|---|---|
| `GDALWarpNoDataMasker` | Marks source pixels matching nodata value |
| `GDALWarpSrcAlphaMasker` | Reads source alpha band as validity |
| `GDALWarpDstAlphaMasker` | Reads destination alpha band |
| `GDALWarpSrcMaskMasker` | Uses .msk sidecar mask files |
| `GDALWarpCutlineMasker` | Applies polygon cutline geometry |

### GDAL's mask model is complex

Three separate mask types propagate through the warper:
- **Per-band validity** (bit masks): per-pixel, per-band valid/invalid
- **Unified validity** (bit mask): per-pixel, all-bands valid/invalid
- **Density** (float): per-pixel blending weight (0.0–1.0)

The density concept enables soft-edge blending at cutlines and partial pixel coverage.

### For a Rust warper

Start with just nodata masking and alpha bands. The density/cutline machinery is
specialized for mosaicking and can come later. Represent validity as a simple bitmask
or `Option<Vec<bool>>`.

---

## Layer 7: The WarpOptions Configuration Surface

GDAL's warper has ~40 string-based options passed via `papszWarpOptions` (a NULL-terminated
string list). The ones that matter for a Rust implementation:

### Essential

| Option | Purpose |
|---|---|
| `NUM_THREADS` | Parallelism |
| `INIT_DEST` | Initialize destination to value or NO_DATA |
| `UNIFIED_SRC_NODATA` | Treat nodata consistently across bands |

### Performance tuning

| Option | Purpose |
|---|---|
| `XSCALE` / `YSCALE` | Override auto-computed resampling scale |
| `SRC_COORD_PRECISION` | Round source coords (reproducibility) |
| `ERROR_THRESHOLD` | Approximate transform tolerance (pixels) |
| `OPTIMIZE_SIZE` | Align chunks to output blocks |
| `SAMPLE_GRID` / `SAMPLE_STEPS` | Source window sampling strategy |

### Specialized / Deferrable

| Option | Purpose |
|---|---|
| `APPLY_VERTICAL_SHIFT` | Vertical datum correction |
| `CUTLINE` / `CUTLINE_BLEND_DIST` | Polygon masking with soft edges |
| `EXCLUDED_VALUES` | Multi-band pixel exclusion (clouds, etc.) |
| `SKIP_NOSOURCE` | Skip output regions with no source coverage |
| `WRITE_FLUSH` | Force disk flush per chunk |
| `STREAMABLE_OUTPUT` | Write output in order |
| `SOURCE_EXTRA` | Extra source pixels beyond computed window |

### The config-in-DuckDB/Parquet idea

These options are currently scattered across:
- Command-line flags in `gdalwarp`
- String list in `GDALWarpOptions::papszWarpOptions`
- XML in VRT files
- Hard-coded defaults

A Parquet/DuckDB configuration store could unify these as a queryable table:

```sql
CREATE TABLE warp_config (
    config_id    TEXT,         -- "sentinel2_to_webmerc_256"
    src_crs      TEXT,         -- "EPSG:32755"
    dst_crs      TEXT,         -- "EPSG:3857"
    resampling   TEXT,         -- "bilinear"
    dst_tile_size INTEGER,     -- 256
    nodata_value DOUBLE,       -- NaN
    num_threads  INTEGER,      -- 4
    error_threshold DOUBLE,    -- 0.125
    extra_options JSON         -- extensible
);
```

This would be consumable by both the Rust warper and GDAL (via a GDAL driver or
configuration reader). It makes warping configurations portable, versionable, and
inspectable — replacing the current pattern of users juggling command-line flags
and scripts.

---

## Layer 8: The DuckDB Cache Architecture

A DuckDB-backed byte-reference store for source data could serve as both:

1. **Tile cache** — persist decoded source tiles keyed by (source_path, tile_x, tile_y, ifd_level)
2. **Warp result cache** — store warped output tiles keyed by (config_id, z, x, y)
3. **Metadata store** — source file inventories, STAC-like catalog entries

```sql
-- Source tile cache
CREATE TABLE tile_cache (
    source_uri  TEXT,
    tile_col    INTEGER,
    tile_row    INTEGER,
    ifd_index   INTEGER,
    fetched_at  TIMESTAMP,
    data        BLOB,          -- compressed tile bytes
    decoded     BLOB           -- or decoded pixel data
);

-- Warped output cache
CREATE TABLE warp_cache (
    config_id   TEXT,
    z           INTEGER,
    x           INTEGER,
    y           INTEGER,
    created_at  TIMESTAMP,
    data        BLOB
);
```

Advantages over GDAL's block cache:
- **Persistent** — survives process restart, shareable between processes
- **Queryable** — "which source tiles have I already fetched?", "how stale is my cache?"
- **Flat** — no complex LRU eviction logic, just `DELETE WHERE fetched_at < ?`
- **Interoperable** — DuckDB reads/writes Parquet natively, so the cache is also a dataset
- **Reusable** — the same cache serves both the Rust warper and GDAL (if GDAL is
  configured to use the same source paths)

---

## Proposed Minimal Viable Rust Warper

### Phase 1: Tile-server MVP (~2,000 lines Rust)

Target: replace the rasterio/GDAL dependency in a titiler-like tile server for
COGs in WebMercator.

- `GridSpec` struct (from vaster concepts)
- `warp_tile()` function: single output tile, single source COG
- Resampling: nearest, bilinear, cubic (generic over data type)
- Transform: affine + PROJ (via `proj` crate)
- Source I/O: async-tiff tiles
- Nodata handling: simple value-based masking
- No chunking needed (output is one tile)

### Phase 2: General warping (~5,000 lines)

- ApproxTransformer equivalent
- Full source window computation with grid sampling
- Statistical resamplers (average, mode, etc.)
- Alpha band handling
- Multi-source mosaicking
- DuckDB cache integration
- Parquet config store

### Phase 3: Advanced features

- GCP/RPC/TPS transformers
- Cutline masking with blend distance
- Vertical datum corrections
- Sum-preserving resampling
- VRT-equivalent lazy composition

---

## Key Design Decisions for Discussion

1. **Trait vs concrete for source I/O?** — `trait SourceReader { async fn read_tiles(...) }`
   allows plugging in async-tiff, Zarr, Parquet byte-refs, or even GDAL itself as a
   fallback.

2. **Per-chunk vs per-pixel scale?** — GDAL's per-chunk scale causes visible
   discontinuities (documented in gdalwarp FAQ Q3). Per-scanline scale is cheap and
   better. Per-pixel is ideal but may not be worth the cost.

3. **How tightly to couple with PROJ?** — `proj4rs` exists as a pure-Rust alternative
   but doesn't cover all coordinate operations. For WASM deployment, pure-Rust is
   essential. For server-side, the `proj` crate (FFI to PROJ C) is fine.

4. **Where does this live?** — New crate in the georust org? Standalone? Needs to be
   somewhere the async-tiff ecosystem and the broader GeoRust community
   can depend on it.

5. **Test oracle** — GDAL's warper output can serve as the ground truth for regression
   testing. Generate reference tiles with gdalwarp, compare with Rust output. The
   `SRC_COORD_PRECISION` option helps make GDAL output reproducible.

---

## Addendum: proj_factors() and the Jacobian

GDAL's `PerformWarp()` computes the local source/destination scale ratio (the Jacobian
of the coordinate transform) by brute force: transforming unit squares and measuring
distortion. But **PROJ itself already exposes this analytically via `proj_factors()`**,
which returns the full set of Tissot indicatrix parameters including meridian/parallel
scale factors, angular distortion, and areal scale.

This means a Rust warper can get exact, analytical per-pixel scale factors directly from
PROJ rather than the sampled/averaged approximation GDAL uses. This would:
- Eliminate the per-chunk scale discontinuities (GDAL FAQ Q3)
- Be cheaper than transforming grids of unit squares
- Be more correct, especially near singularities (poles, projection boundaries)

The `proj` crate would need to expose `proj_factors()` (currently it may not), but this
is a straightforward FFI addition. For a pure-Rust path, the scale factors can be
computed from finite differences of the transform, which is what GDAL does anyway.

This directly connects to the `tissot` package work — the same Tissot parameters that
describe cartographic distortion are exactly the resampling scale factors the warper needs.

---

## Addendum: PROJ is now WASMable (Feb 2026)

Javier Jimenez Shaw announced on 2026-02-23 that PROJ (current master) compiles to
WebAssembly via Emscripten and runs in the browser:

- Transform: <https://jjimenezshaw.github.io/wasm-proj/transform.html>
- Projinfo: <https://jjimenezshaw.github.io/wasm-proj/projinfo.html>

Key details:
- Uses the WASM compilation already in PROJ CI
- Grid files download on demand at runtime (using PROJ's existing network capability)
- The proj.db can be embedded in the binary (since PROJ 9.6.0)
- Runs entirely client-side, no server needed

**This changes the WASM story significantly for a Rust warper.** Previously, the argument
for pure-Rust coordinate transforms (proj4rs) in WASM was that PROJ's C++ couldn't easily
target WASM. Now you could have:

1. **Server-side**: Rust warper + PROJ via FFI (`proj` crate) — full accuracy
2. **WASM/browser**: Rust warper + PROJ via Emscripten WASM — same accuracy, runs in browser
3. **Lightweight WASM**: Rust warper + proj4rs (pure Rust) — smaller bundle, less coverage

Option 2 is particularly interesting because it means **the exact same coordinate
transform pipeline** works on server and in browser. Existing client-side tools like
deck.gl-raster currently use proj4js for browser reprojection; replacing it with actual
PROJ would eliminate the accuracy differences between server and client rendering.

The GeoRust community should know about this — it removes what was
previously considered a major blocker for PROJ in browser-based geospatial tools.

---

## Addendum: DuckDB Cache with Byte-Referenced Blocks

The DuckDB cache idea extends beyond simple key-value tile storage. DuckDB can hold
**byte-referenced compressed blocks** — rows that contain both the raw compressed tile
data and metadata about where it came from:

```sql
CREATE TABLE block_store (
    source_uri   TEXT,          -- "s3://bucket/scene.tif"
    ifd_index    SMALLINT,      -- overview level
    tile_col     INTEGER,
    tile_row     INTEGER,
    byte_offset  BIGINT,        -- offset in source file
    byte_length  INTEGER,       -- compressed size
    compression  TEXT,          -- "deflate", "jpeg", "webp"
    raw_bytes    BLOB,          -- the actual compressed bytes
    fetched_at   TIMESTAMP,
    etag         TEXT           -- for cache invalidation
);
```

The key insight: this is a **flat, queryable, portable representation** of what GDAL's
block cache holds in an opaque in-memory hash table. You can:

- **Export it as Parquet** for archival or sharing (`COPY block_store TO 'cache.parquet'`)
- **Query what you have** before fetching (`SELECT ... WHERE source_uri = ? AND tile_col BETWEEN ? AND ?`)
- **Share between processes** (DuckDB supports concurrent readers)
- **Warm the cache** from a separate process while the tile server is running
- **Inspect cache hit rates** with SQL aggregations
- **Prune intelligently** (`DELETE WHERE fetched_at < NOW() - INTERVAL '7 days'`)

The compressed blocks can be decoded lazily — only when the warper actually needs the
pixels. This means the cache stores data in its most compact form and decompression
happens at the last possible moment.

For the warp config store, the same DuckDB database can hold both the cache and the
configurations, making a single `.duckdb` file a complete, self-contained, portable
representation of a warping pipeline's state.

---

## References

- GDAL source: <https://github.com/OSGeo/gdal> (trunk as of 2026-02-28)
- async-tiff: <https://github.com/developmentseed/async-tiff>
- warp-resample-profiling: <https://developmentseed.org/warp-resample-profiling/>
- GeoRust: <https://georust.org/>
- proj4rs: <https://github.com/3liz/proj4rs>
- vaster: hypertidy ecosystem (grid logic for raster operations)
- PROJ in WASM: <https://lists.osgeo.org/pipermail/proj/2026-February/011970.html>
- PROJ WASM issue: <https://github.com/OSGeo/PROJ/issues/4327>
- proj_factors: <https://proj.org/en/stable/development/reference/functions.html#c.proj_factors>
