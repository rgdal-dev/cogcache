# vwarp — Virtual Warp Planning and Execution

**Michael Sumner — March 2026**

> "reproject(data, crs) is a mistake"

## The core insight

Every reprojection function in R — `terra::project()`, `sf::st_transform()`,
`stars::st_warp()`, `raster::projectRaster()` — takes **data** and a **target CRS**,
then guesses the rest: the extent, the resolution, the dimensions. Those heuristics
are doing `ComputeSourceWindow` backwards: inferring the target grid from the source,
when the user always really knows (or should know) what output grid they want.

The defaults are convenient. They're also the source of almost every "why does my
reprojected raster have weird dimensions / wrong resolution / unexpected extent"
question on Stack Overflow.

**vwarp's philosophy: you declare the target grid explicitly.**

A target grid is four things: `bbox(extent)`, `dim`, `crs`, `chunking`. That's the
minimum, even when heuristics fill in the blanks. Once you have a target grid, everything
flows from it:

```
target_chunk_spec → {source_chunk_refs} → read → warp → write
```

The *mapping* from target chunks to source chunks is where all the insight lives.
The actual warping is mechanical.

## What we've proven

The cogcache prototype demonstrates the full pipeline end-to-end, with
**bit-identical results to GDAL** (`gdalwarp -et 0 -ovr NONE -r near`):

| Component | Rust module | GDAL equivalent | Status |
|-----------|------------|-----------------|--------|
| Coordinate transform | `transform.rs` | `GDALGenImgProjTransform` | ✅ bit-identical |
| Approximate transform | `approx.rs` | `GDALApproxTransform` | ✅ bit-identical |
| Source window computation | `source_window.rs` | `ComputeSourceWindow` | ✅ matches GDAL |
| Nearest-neighbour kernel | `warp.rs` | `GWKNearestThread` | ✅ bit-identical |
| Tile decode (DEFLATE+pred2) | `lib.rs` | GDAL block cache | ✅ byte-identical |

**Validated on real data:**

- Sentinel-2 B04 over Macquarie Island (256×256, EPSG:3857 → EPSG:3031): 65536/65536 match
- IBCSO v2 bathymetry, Mawson Station (EPSG:9354 → LCC): 65536/65536 match
- GEBCO 2024 global (86400×43200, WGS84 → Fiji LCC): 48/48 tiles, 65536/65536 centre tile

The Fiji demo exercises the hardest case: a conic projection centred on the
antimeridian, warping from a global source. The source window computation handles
the 90% heuristic correctly and the per-tile architecture resolves each tile's
source footprint independently.

## Architecture

### The monolith we're decomposing

GDAL's warp pipeline (`gdalwarp` → `GDALWarpOperation`) bundles:

1. **Planning** — what source data is needed? (`ComputeSourceWindow`)
2. **Chunking** — how to subdivide into memory-fitting pieces? (`CollectChunkList`)
3. **I/O** — read source, write dest (`GDALRasterIO`)
4. **Transform** — coordinate mapping (`GenImgProjTransformer`)
5. **Approximation** — interpolate coordinates (`ApproxTransformer`)
6. **Kernel** — resample pixels (`GWKNearestThread` etc.)

These are tangled inside `GDALWarpOperation`. You can't ask "what would you read?"
without actually reading. You can't compare two target grid specs without running two
full warps. The planning step is invisible.

### vwarp separates them

```
┌─────────────────────────────────────────────────┐
│  vwarp                                          │
│                                                 │
│  ┌──────────────────┐   ┌────────────────────┐  │
│  │  Target grid     │   │  Source index       │  │
│  │  spec            │   │  (byte-ref table)   │  │
│  │  (extent,dim,    │   │                     │  │
│  │   crs,chunking)  │   │  COG tile offsets   │  │
│  └────────┬─────────┘   │  Zarr chunk refs    │  │
│           │             │  NetCDF byte ranges  │  │
│           ▼             │  VirtualiZarr store  │  │
│  ┌──────────────────┐   │  Parquet catalogue   │  │
│  │  Plan            │   └─────────┬──────────┘  │
│  │                  │             │              │
│  │  For each target │◄────────────┘              │
│  │  chunk, compute  │                            │
│  │  source window   │   Uses: Transformer        │
│  │  → {src refs}    │         ComputeSourceWindow │
│  └────────┬─────────┘                            │
│           │                                      │
│           ▼         (review / adjust / iterate)  │
│  ┌──────────────────┐                            │
│  │  Execute         │                            │
│  │                  │   Uses: ApproxTransformer   │
│  │  read bytes →    │         warp kernel         │
│  │  decode →        │                            │
│  │  warp →          │                            │
│  │  write           │                            │
│  └──────────────────┘                            │
└─────────────────────────────────────────────────┘
```

**The plan is the product.** Before touching any pixels, you can:

- See how many source blocks each target tile needs
- Identify duplication (the same source block read by multiple tiles)
- Spot antimeridian or pole problems (source window = full raster width)
- Compare fill ratios to find tiles worth splitting
- Estimate total I/O cost
- Try a different target grid and compare instantly

### The source index

Source data is just byte references. A source chunk is:

```
(byte_offset, byte_length, tile_x, tile_y, band, time, spatial_footprint)
```

These come from anywhere:

- **COG**: IFD scanning gives tile offsets/lengths + geotransform gives footprints
- **Zarr**: chunk manifest (key → byte range)
- **NetCDF**: variable/chunk layout from metadata
- **VirtualiZarr**: the `.json` or `.parq` store
- **Parquet catalogue**: our enriched byte-ref table (rustycogs / wholecat)

The format doesn't matter. What matters is that every chunk is a lazy rectangle
with a spatial footprint and a pointer to encoded bytes. That's a vector layer.

### The meta-vector-driver insight

The mapping `target_chunk → {source_chunks}` is a spatial join between two sets
of rectangles across CRS boundaries. Target chunks have footprints in the target
CRS. Source chunks have footprints in the source CRS. The join is mediated by the
coordinate transformer.

This is why a byte-ref parquet with a geometry column isn't just a storage format — 
it's a spatial index of array data. You can query it spatially, join it against
target specs, filter by time/band, and the result is a work plan.

### N-dimensional generality

Warping is inherently 2D — it operates on the spatial axes. But chunks are n-dimensional.
A NetCDF chunk indexed by `(time, depth, lat, lon)` has a 2D spatial footprint and
n-2 slice coordinates. The warp plan maps the spatial part; the slice part passes
through. "x,y" and "lon,lat" are just conventions about which axis pair carries
the spatial role. Nothing prevents the spatial axes being axes 3 and 4 of a 5D array — 
it's just less common.

## Rust crate structure

The cogcache prototype has everything in one R package. The Rust modules are designed
to be extractable as standalone crates:

| Module | Responsibility | Dependencies |
|--------|---------------|-------------|
| `transform` | `Transformer` trait, `GenImgProjTransformer` | `proj` crate |
| `approx` | `ApproxTransformer` (adaptive interpolation) | `transform` |
| `source_window` | `ComputeSourceWindow` (planning) | `transform` |
| `warp` | Nearest-neighbour kernel | `transform` |
| `vaster` | Grid arithmetic (already a separate crate) | none |

The `vaster` Rust crate is published. The others live in cogcache for now but have
no cross-dependencies except through the `Transformer` trait.

**Future crate organisation:**

```
vwarp/                          # workspace
├── vaster/                     # grid arithmetic (published)
├── vwarp-transform/            # Transformer trait + GenImgProjTransformer
├── vwarp-approx/               # ApproxTransformer
├── vwarp-plan/                 # ComputeSourceWindow, chunk planning
├── vwarp-kernel/               # warp kernels (nearest, bilinear, cubic, lanczos)
└── vwarp/                      # top-level: orchestration, R bindings
```

## Relationship to hypertidy packages

| Package | Role | Overlap with vwarp |
|---------|------|--------------------|
| `vaster` | Grid logic (R + Rust) | Foundation — extent/dim/gt math |
| `grout` | Tile index ↔ pixel coords | Source index construction |
| `vapour` | GDAL feature/raster access | Data reader (consumed by vwarp) |
| `gdalraster` | GDAL raster API | Data reader, validation oracle |
| `ndr` | xarray-like lazy arrays | Consumer of vwarp outputs |
| `tissot` | Projection distortion | Complementary — same grid math |
| `geographiclib` | Geodesic calculations | Independent |
| `somap2` | Southern Ocean mapping | Consumer |

**vwarp does not replace vapour or gdalraster.** It replaces the invisible planning
step inside `gdalwarp` with an explicit, inspectable, composable alternative. You
still use GDAL to read pixels. You use vwarp to decide *which* pixels to read.

## What's next

### Immediate (proven, needs polish)

1. **Clean up cogcache → vwarp rename**: the R package, Rust workspace, documentation
2. **Archive organic development docs** to `inst/doc-archive/`
3. **Write the vignette**: "reproject(data, crs) Is a Mistake"
4. **CRAN-ready vaster** (R crate): already close, needs final checks

### Short-term (architecture is clear, implementation straightforward)

5. **Bilinear/cubic/lanczos kernels**: the scanline loop is identical, just different
   source pixel sampling. The `source_window` already takes `resample_padding`.
6. **Overview selection**: match GDAL's automatic zoom-level selection when downsampling
7. **Antimeridian-aware source splitting**: detect when source window wraps and
   issue two reads instead of reading the full raster width
8. **Pole detection**: the `aDstXYSpecialPoints` mechanism from GDAL — transform
   pole coordinates to dest pixel space and check if they fall in the dest window

### Medium-term (design validated, needs prototyping)

9. **Byte-ref parquet as source index**: enriched with geometry column, queryable
   with DuckDB spatial extension. Proof of concept exists in `wholecat/`.
10. **Multi-band support**: the warp kernel operates on one band at a time, but
    the planning step is band-independent. Read all bands for a source window,
    warp each.
11. **Parallel tile execution**: the per-tile architecture is embarrassingly parallel.
    Each tile's plan is independent. `rayon` in Rust or `future` in R.
12. **Cloud-native I/O**: HTTP range requests for COG tiles, direct Zarr chunk
    reads. The byte-ref table tells you exactly which ranges to fetch.

### Long-term (vision, needs design work)

13. **VirtualiZarr interop**: consume `.json`/`.parq` virtual stores as source indices
14. **The meta-vector-driver**: expose chunk footprints as an OGR-compatible vector
    layer. Query: "give me all GHRSST chunks that intersect this Antarctic polygon
    for January 2024" → work plan → warp.
15. **Streaming/incremental warps**: for the estinel pipeline — new source data arrives,
    the plan identifies which target tiles are affected, only those get reprocessed.
16. **GPU kernels**: the warp kernel is pure arithmetic on arrays. The Rust crate
    structure makes it possible to swap in a GPU implementation without touching
    the planning layer.

## GDAL source lineage

All Rust code in this project is written from scratch with reference to GDAL's
algorithms and data structures. The GDAL source commit used for analysis is
`a3b7b01d3e` (GDAL 3.13.0dev, released 2026-02-23). The specific files referenced:

| GDAL source file | Lines | Our implementation |
|-----------------|-------|-------------------|
| `gdaltransformer.cpp` | ~L2800, ~L3500, ~L4113 | `transform.rs`, `approx.rs` |
| `gdalwarpoperation.cpp` | ~L2754, ~L3015 | `source_window.rs` |
| `gdalwarpkernel.cpp` | ~L5206, ~L5510 | `warp.rs` |
| `gdalwarp_lib.cpp` | ~L2832, ~L3260 | `lib.rs` (orchestration) |

GDAL is MIT-licensed. Our code follows the same algorithms but is an independent
Rust implementation, not a translation or port of the C++ source.
