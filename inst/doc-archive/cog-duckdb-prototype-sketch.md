# COG → DuckDB Cache Prototype Sketch

**Goal**: Prove the concept that a COG's internal tile structure can be faithfully catalogued in DuckDB, compressed tile bytes fetched and stored as blobs, and the results verified against what GDAL reads internally. This is the evidence layer for the Rust warper architecture.

## What we're building

A self-contained R script (or Quarto doc) that:

1. Reads a COG's IFD/tile layout into a DuckDB table (the "catalog")
2. Fetches raw compressed tile bytes via HTTP range requests
3. Stores them as BLOBs in DuckDB (the "cache")
4. Reads the same pixels via GDAL (vapour / terra / gdalwarp)
5. Decodes our cached bytes and compares — proving byte-level correspondence

The deliverable is a `.duckdb` file you can hand someone and say "this *is* that COG's tile structure, and here are the actual bytes."

## Prerequisites

```r
# Packages we need
library(duckdb)       # cache store
library(vapour)       # GDAL access (info, raster reads, VSICURL)
library(httr2)        # HTTP range requests (or curl)
library(jsonlite)     # parsing GDAL info output

# Maybe:
library(terra)        # for comparison reads and visual checks
library(wk)           # if we want to express tile footprints as geometries
```

### Why these packages

- **vapour** over terra for the core work: we want raw access to GDAL's view of the file, not a high-level raster object. `vapour_raster_info()` gives us the structural metadata. `vapour_read_raster_raw()` gives us pixel values from specific windows for comparison.
- **duckdb** over RSQLite: native Parquet export, BLOB handling, better analytical queries for inspecting the cache later.
- **httr2** over base `download.file`: we need HTTP Range headers for byte-range fetches. Could also use `curl` package directly.
- **terra** only for visual sanity checks — plotting tiles to confirm they look right.

## Step 0: Pick a COG

Use something you already work with in estinel. A Sentinel-2 L2A COG from Element 84's open bucket is ideal:

```r
# Example — a Sentinel-2 B04 (red) band COG on AWS
cog_url <- "/vsicurl/https://sentinel-cogs.s3.us-west-2.amazonaws.com/sentinel-s2-l2a-cogs/55/H/FA/2024/1/S2B_55HFA_20240115_0_L2A/B04.tif"

# Or a local COG if you have one — avoids network latency during dev
# cog_url <- "path/to/local.tif"
```

Pick a tile that covers somewhere familiar (Hobart-adjacent Sentinel-2 granule?) so you can visually sanity-check.

## Step 1: Build tile catalog with grout + gdalraster

**grout** gives us the tile grid geometry from dimension + block size — this is the
catalog of *which* tiles exist and their pixel/geo extents. **gdalraster** gives us
the byte offsets via the GTiff driver's `BLOCK_OFFSET` metadata items.

```r
library(grout)
library(gdalraster)
library(duckdb)

# Open via gdalraster for structural metadata + byte offsets
ds <- new(GDALRaster, cog_url)

dim_xy <- c(ds$getRasterXSize(), ds$getRasterYSize())
block_size <- ds$getBlockSize(band = 1)  # e.g. c(512, 512)
gt <- ds$getGeoTransform()
crs_wkt <- ds$getProjectionRef()

# grout builds the tile grid — handles edge tiles, geo extents, everything
tile_spec <- grout(dim_xy, block_size)
# tile_spec already knows tile_col, tile_row, pixel offsets, geo extents
```

## Step 2: Get byte offsets from GDAL's TIFF metadata domain

The GTiff driver exposes per-block byte offsets as band-level metadata in the
`"TIFF"` domain. Each tile's offset is at `BLOCK_OFFSET_{xblock}_{yblock}` and
byte count at `BLOCK_SIZE_{xblock}_{yblock}`. gdalraster's `$getMetadataItem()`
gives us direct access.

```r
# Get byte offset and size for a specific block
get_block_byte_ref <- function(ds, band, tile_col, tile_row) {
  offset_key <- sprintf("BLOCK_OFFSET_%d_%d", tile_col, tile_row)
  size_key <- sprintf("BLOCK_SIZE_%d_%d", tile_col, tile_row)
  
  offset <- ds$getMetadataItem(band = band, mdi_name = offset_key, domain = "TIFF")
  bcount <- ds$getMetadataItem(band = band, mdi_name = size_key, domain = "TIFF")
  
  # These return character strings, or "" if not available
  list(
    byte_offset = if (nzchar(offset)) as.numeric(offset) else NA_real_,
    byte_length = if (nzchar(bcount)) as.integer(bcount) else NA_integer_
  )
}

# Test on first tile
get_block_byte_ref(ds, band = 1, tile_col = 0, tile_row = 0)
# Should return something like list(byte_offset = 32768, byte_length = 148576)
```

Now build the full catalog with byte refs:

```r
# Build catalog from grout's tile grid + gdalraster byte offsets
n_tiles_x <- ceiling(dim_xy[1] / block_size[1])
n_tiles_y <- ceiling(dim_xy[2] / block_size[2])

tile_catalog <- expand.grid(
  tile_col = seq_len(n_tiles_x) - 1L,
  tile_row = seq_len(n_tiles_y) - 1L
)

# Pixel extents (grout does this — use grout's output if it gives a data.frame)
tile_catalog$pixel_x_off <- tile_catalog$tile_col * block_size[1]
tile_catalog$pixel_y_off <- tile_catalog$tile_row * block_size[2]
tile_catalog$pixel_x_size <- pmin(block_size[1], dim_xy[1] - tile_catalog$pixel_x_off)
tile_catalog$pixel_y_size <- pmin(block_size[2], dim_xy[2] - tile_catalog$pixel_y_off)

# Geo extents from geotransform
tile_catalog$geo_xmin <- gt[1] + tile_catalog$pixel_x_off * gt[2]
tile_catalog$geo_xmax <- gt[1] + (tile_catalog$pixel_x_off + tile_catalog$pixel_x_size) * gt[2]
tile_catalog$geo_ymax <- gt[4] + tile_catalog$pixel_y_off * gt[6]
tile_catalog$geo_ymin <- gt[4] + (tile_catalog$pixel_y_off + tile_catalog$pixel_y_size) * gt[6]

# Byte offsets from GDAL's TIFF metadata domain
byte_refs <- mapply(
  function(tc, tr) get_block_byte_ref(ds, band = 1, tc, tr),
  tile_catalog$tile_col, tile_catalog$tile_row,
  SIMPLIFY = FALSE
)
tile_catalog$byte_offset <- vapply(byte_refs, `[[`, numeric(1), "byte_offset")
tile_catalog$byte_length <- vapply(byte_refs, `[[`, integer(1), "byte_length")

tile_catalog$source_uri <- cog_url
tile_catalog$ifd_index <- 0L

# Into DuckDB
con <- dbConnect(duckdb())
dbWriteTable(con, "tile_catalog", tile_catalog)

# Quick sanity check
dbGetQuery(con, "
  SELECT count(*) as n_tiles,
         count(byte_offset) as n_with_offsets,
         sum(byte_length) as total_compressed_bytes,
         min(byte_offset) as first_offset,
         max(byte_offset + byte_length) as last_byte
  FROM tile_catalog
")
```

**Note**: The `BLOCK_OFFSET` / `BLOCK_SIZE` metadata items are documented in the
GDAL GTiff driver docs. If gdalraster's `$getMetadataItem()` with domain `"TIFF"`
doesn't work, try domain `""` — or check if vapour exposes this. The key names
are definitely `BLOCK_OFFSET_{x}_{y}` and `BLOCK_SIZE_{x}_{y}`.

## Step 3: Fetch and store compressed tile bytes

```r
library(blob)

# Fetch a specific tile's raw compressed bytes via HTTP range request
fetch_tile_bytes <- function(url, byte_offset, byte_length) {
  # Strip /vsicurl/ prefix for HTTP access
  http_url <- sub("^/vsicurl/", "", url)
  
  req <- httr2::request(http_url) |>
    httr2::req_headers(
      Range = sprintf("bytes=%d-%d", byte_offset, byte_offset + byte_length - 1)
    )
  
  resp <- httr2::req_perform(req)
  httr2::resp_body_raw(resp)
}

# Add blob column to the catalog
dbExecute(con, "ALTER TABLE tile_catalog ADD COLUMN raw_bytes BLOB")

# Fetch a subset of tiles (not all — just enough to prove the concept)
sample_tiles <- dbGetQuery(con, "
  SELECT tile_col, tile_row, byte_offset, byte_length
  FROM tile_catalog
  WHERE byte_offset IS NOT NULL
  ORDER BY random()
  LIMIT 20
")

# Fetch all tiles, collect into a list
fetched <- lapply(seq_len(nrow(sample_tiles)), function(i) {
  raw <- fetch_tile_bytes(
    cog_url,
    sample_tiles$byte_offset[i],
    sample_tiles$byte_length[i]
  )
  list(tile_col = sample_tiles$tile_col[i],
       tile_row = sample_tiles$tile_row[i],
       raw_bytes = raw)
})

# Build a data.frame with as_blob() — note: blob() won't work, must use as_blob()
raw_list <- lapply(fetched, function(x) x$raw_bytes)

blob_df <- data.frame(
  tile_col = vapply(fetched, `[[`, integer(1), "tile_col"),
  tile_row = vapply(fetched, `[[`, integer(1), "tile_row")
)
blob_df$raw_bytes <- as_blob(raw_list)

# Bulk write to DuckDB and update the catalog
dbWriteTable(con, "tile_blobs", blob_df, overwrite = TRUE)

dbExecute(con, "
  UPDATE tile_catalog SET raw_bytes = tile_blobs.raw_bytes
  FROM tile_blobs
  WHERE tile_catalog.tile_col = tile_blobs.tile_col
    AND tile_catalog.tile_row = tile_blobs.tile_row
    AND tile_catalog.ifd_index = 0
")

# Verify: blob sizes must match byte_length exactly
dbGetQuery(con, "
  SELECT tile_col, tile_row, octet_length(raw_bytes) as blob_size, byte_length
  FROM tile_catalog
  WHERE raw_bytes IS NOT NULL
  LIMIT 5
")
# blob_size == byte_length for every row confirms the bytes are intact
```

**Key gotchas discovered**:
- Use `as_blob()` not `blob()` — `blob()` rejects a list of raw vectors, `as_blob()` accepts it
- DuckDB uses `octet_length()` for BLOB size, not `length()` (which only works on VARCHAR/BIT/arrays)
- The raw vectors from `httr2::resp_body_raw()` are plain `"raw"` class, which `as_blob()` handles fine
```

## Step 4: Decode and compare with GDAL

This is the payoff — proving that our cached bytes decode to the same pixels GDAL reads.

**`vrt://?block=i,j`** (GDAL 3.13) is the clean way to read a single block. No manual
window arithmetic, no edge-tile size calculations — GDAL handles it all via
`GDALRasterBand::ReadBlock()`.

```r
# Read a specific block via vrt://?block= (GDAL 3.13)
read_block_via_gdal <- function(cog_url, tile_col, tile_row, block_size) {
  block_url <- sprintf("vrt://%s?block=%d,%d", cog_url, tile_col, tile_row)
  
  # vapour reads the block as a subwindow
  vapour_read_raster_raw(block_url, 
    window = c(0, 0, block_size[1], block_size[2]),
    dimension = block_size)
  
  # Or with gdalraster:
  # ds_block <- new(GDALRaster, block_url)
  # ds_block$read(band = 1, xoff = 0, yoff = 0,
  #               xsize = ds_block$getRasterXSize(),
  #               ysize = ds_block$getRasterYSize(),
  #               out_xsize = ds_block$getRasterXSize(),
  #               out_ysize = ds_block$getRasterYSize())
}

# Decode our cached bytes
# For DEFLATE-compressed tiles:
decode_tile <- function(raw_bytes, compression, pixel_count, nbytes_per_pixel = 2L) {
  if (compression == "deflate") {
    memDecompress(raw_bytes, type = "gzip")
  } else if (compression == "none") {
    raw_bytes
  } else {
    stop("Unsupported compression: ", compression)
  }
}

# Compare!
compare_tile <- function(con, tile_col, tile_row, cog_url, block_size) {
  # Our cached bytes
  cached <- dbGetQuery(con, sprintf(
    "SELECT raw_bytes, compression FROM tile_catalog
     WHERE tile_col = %d AND tile_row = %d AND raw_bytes IS NOT NULL",
    tile_col, tile_row
  ))
  
  decoded_ours <- decode_tile(cached$raw_bytes[[1]], cached$compression)
  
  # GDAL's read
  gdal_pixels <- read_tile_via_gdal(cog_url, tile_col, tile_row, block_size)
  
  # Compare raw bytes
  identical(decoded_ours, gdal_pixels[[1]])
}
```

If `identical()` returns TRUE, we've proven the correspondence. The bytes in DuckDB *are* the pixels GDAL reads.

## Step 5: The DuckDB cache as a portable artifact

```r
# Export to Parquet — the cache is now a shareable dataset
dbExecute(con, "
  COPY (
    SELECT source_uri, ifd_index, tile_col, tile_row,
           pixel_x_off, pixel_y_off, pixel_x_size, pixel_y_size,
           geo_xmin, geo_ymax, geo_xmax, geo_ymin,
           byte_offset, byte_length, compression,
           length(raw_bytes) as cached_bytes_size
    FROM tile_catalog
  ) TO 'tile_catalog.parquet' (FORMAT PARQUET)
")

# The .duckdb file itself is the complete artifact
dbDisconnect(con, shutdown = TRUE)

# This file can be:
# - Opened by anyone with DuckDB (R, Python, CLI, Rust)
# - Queried: "which tiles cover this geographic extent?"
# - Used to warm a Rust warper's cache without touching S3
# - Exported as Parquet for archival
# - Extended with warp_config table (Phase 2)
```

## Step 6: Visual comparison (the demo slide)

```r
# Side by side: GDAL read vs our cache decode, for a single tile
par(mfrow = c(1, 2))

# GDAL
gdal_pix <- read_tile_via_gdal(cog_url, 5, 5, block_size)
image(matrix(gdal_pix[[1]], block_size[1], block_size[2]),
      main = "GDAL read", col = grey.colors(256))

# Our cache
cached <- dbGetQuery(con, "SELECT raw_bytes, compression FROM tile_catalog
                           WHERE tile_col = 5 AND tile_row = 5")
our_pix <- decode_tile(cached$raw_bytes[[1]], cached$compression)
image(matrix(our_pix, block_size[1], block_size[2]),
      main = "DuckDB cache", col = grey.colors(256))

# They should be visually identical — and byte-identical
```

## What this proves

When this works, we've demonstrated:

1. **COG structure is transparent** — the tile grid is a queryable catalog, not opaque internal state
2. **Bytes are portable** — compressed tile data moves from S3 → DuckDB → decoded pixels with no GDAL in the middle (for the decode step)
3. **Cache is inspectable** — SQL queries answer "what's cached?", "how much?", "which tiles cover this bbox?"
4. **Cache is shareable** — `.duckdb` or `.parquet` file replaces "you need access to this S3 bucket"
5. **GDAL agreement** — our independent byte-level access produces identical pixels to GDAL's internal pipeline

## Known unknowns / things to figure out

- **Compression codecs**: Sentinel-2 COGs use DEFLATE. JPEG-compressed COGs (common for RGB) need libjpeg decode. WEBP needs libwebp. For the prototype, start with DEFLATE-only.
- **Predictor**: TIFF DEFLATE often uses horizontal differencing predictor (predictor=2). Need to undo this after decompression. Check the TIFF tag. gdalraster should expose this via metadata too.
- **Byte order**: TIFF stores data in the byte order indicated by the header. Sentinel-2 is typically little-endian UInt16.
- **Overview IFDs**: COGs have multiple IFDs (full res, 2x, 4x, ...). The tile catalog should cover all levels. For overviews, open via `vrt://{url}?overview_level=N` and repeat the `BLOCK_OFFSET` queries.
- **`BLOCK_OFFSET` metadata domain**: The GDAL docs say these are in the `"TIFF"` metadata domain, but test with `""` too. The key names are definitely `BLOCK_OFFSET_{x}_{y}` and `BLOCK_SIZE_{x}_{y}`.
- **`vrt://?block=` with vapour**: Need to confirm vapour handles the `vrt://` connection string. If not, gdalraster definitely does. Could also use `gdal_translate` on the command line.
- **Nodata handling**: Sentinel-2 L2A uses 0 as nodata. Check that our decode handles this the same way GDAL does (it should — nodata is a metadata concept, not a byte-level transform).
- **grout output format**: Check what grout returns — if it's a data.frame or tbl with tile extents already computed, use that directly instead of the manual `expand.grid` + geotransform arithmetic.

## Stretch goals for the prototype

- [ ] Parse multiple IFD levels (overviews) into the catalog
- [ ] Compute which output web tiles (z/x/y) need which source tiles
- [ ] Store a warp_config row: source CRS, target EPSG:3857, resampling=bilinear
- [ ] Use PROJ (via sf::st_transform or direct proj4 call) to transform tile corner coordinates
- [ ] Compare our tile→webmerc mapping with what GDAL computes

## File layout

```
cog-cache-prototype/
├── 01-catalog.R          # Build tile catalog from COG
├── 02-fetch-bytes.R      # HTTP range fetch, store in DuckDB
├── 03-compare-gdal.R     # Decode and compare with GDAL reads
├── 04-visualize.R        # Side-by-side plots
├── cache.duckdb          # The output artifact
├── tile_catalog.parquet  # Exported catalog (no blobs)
└── README.md
```

Or as a single Quarto doc: `cog-duckdb-cache.qmd` with all steps in sequence.
