## Full pipeline: Rust warp map + Rust tile decode = warped output, no GDAL in pixel path
##
## This connects:
##   1. rust_warp_map()          — which source pixels do we need?
##   2. tile catalog (DuckDB)    — which tiles contain those pixels?
##   3. rust_fetch_decode_tile() — fetch & decode those tiles
##   4. rust_apply_warp()        — sample decoded pixels into dest grid
##
## Run after: rextendr::document() and devtools::load_all()

library(gdalraster)
library(reproj)
library(vaster)
library(slippymath)
library(duckdb)

## ============================================================================
## Step 1: Grid specs (same as before)
## ============================================================================

cog_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"

ds <- new(GDALRaster, cog_url)

src_ext <- ds$bbox()[c(1, 3, 2, 4)]
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
block_size <- ds$getBlockSize(band = 1)

tile_extent_3857 <- function(z, x, y) {
  n <- 2^z
  lon <- c(x, x + 1) / n * 360 - 180
  lat <- atan(sinh(pi * (1 - 2 * c(y, y + 1) / n))) * 180 / pi
  reproj_extent(c(lon, rev(lat)), "EPSG:3857", source = "EPSG:4326")
}

llex <- reproj_extent(src_ext, "EPSG:4326", source = src_crs)
z <- 10L
tile <- lonlat_to_tilenum(mean(llex[1:2]), mean(llex[3:4]), z)

dest_ext <- tile_extent_3857(z, tile$x, tile$y)
dest_dim <- c(256L, 256L)
dest_crs <- srs_to_wkt("EPSG:3857")
dest_gt  <- extent_dim_to_gt(dest_ext, dest_dim)

## ============================================================================
## Step 2: Compute warp map (Rust)
## ============================================================================

cat("=== Computing warp map ===\n")
wmap <- rust_warp_map(
  src_crs = src_crs,
  src_gt  = src_gt,
  src_dim = as.integer(src_dim),
  dst_crs = dest_crs,
  dst_gt  = dest_gt,
  dst_dim = as.integer(dest_dim)
)



valid <- !is.na(wmap$src_cols)
cat("Valid pixels:", sum(valid), "/", length(valid), "\n")

## ============================================================================
## Step 3: Figure out which source tiles we need
## ============================================================================

cat("\n=== Tile lookup ===\n")

## Source pixel ranges (1-based from warp map)
src_col_range <- range(wmap$src_cols[valid])
src_row_range <- range(wmap$src_rows[valid])

## Which tile indices do these pixels fall in? (0-based tile indices)
tile_col_range <- c(
  floor((src_col_range[1] - 1) / block_size[1]),
  floor((src_col_range[2] - 1) / block_size[1])
)
tile_row_range <- c(
  floor((src_row_range[1] - 1) / block_size[2]),
  floor((src_row_range[2] - 1) / block_size[2])
)

needed_tiles <- expand.grid(
  tile_col = seq(tile_col_range[1], tile_col_range[2]),
  tile_row = seq(tile_row_range[1], tile_row_range[2])
)

cat("Source pixel range: cols", src_col_range, "rows", src_row_range, "\n")
cat("Block size:", block_size, "\n")
cat("Tiles needed:", nrow(needed_tiles), "\n")
cat("Tile col range:", tile_col_range, "\n")
cat("Tile row range:", tile_row_range, "\n")

## ============================================================================
## Step 4: Get byte offsets from GDAL metadata + fetch/decode tiles
## ============================================================================

cat("\n=== Fetching and decoding tiles ===\n")

## Get predictor
predictor_str <- ds$getMetadataItem(band = 0, mdi_name = "PREDICTOR",
                                     domain = "IMAGE_STRUCTURE")
predictor <- if (nzchar(predictor_str)) as.integer(predictor_str) else 1L
cat("Predictor:", predictor, "\n")

## Strip /vsicurl/ for HTTP access
http_url <- sub("^/vsicurl/", "", cog_url)

## Build a source pixel buffer covering all needed tiles
## This is the "assembled tiles" approach: stitch tiles into a
## contiguous buffer, then apply_warp indexes into it

## Pixel extent of all needed tiles (0-based)
buf_col_start <- tile_col_range[1] * block_size[1]  # 0-based pixel
buf_row_start <- tile_row_range[1] * block_size[2]
buf_ncol <- (tile_col_range[2] - tile_col_range[1] + 1) * block_size[1]
buf_nrow <- (tile_row_range[2] - tile_row_range[1] + 1) * block_size[2]

## Clamp to actual image dimensions
buf_ncol <- min(buf_ncol, src_dim[1] - buf_col_start)
buf_nrow <- min(buf_nrow, src_dim[2] - buf_row_start)

cat("Tile buffer:", buf_ncol, "x", buf_nrow, "pixels\n")

## Allocate the buffer (fill with 0 = nodata)
nodata <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- 0L
tile_buf <- rep(as.integer(nodata), buf_ncol * buf_nrow)

t0 <- proc.time()

for (i in seq_len(nrow(needed_tiles))) {
  tc <- needed_tiles$tile_col[i]
  tr <- needed_tiles$tile_row[i]

  ## Get byte offset from GDAL metadata
  offset_key <- sprintf("BLOCK_OFFSET_%d_%d", tc, tr)
  size_key   <- sprintf("BLOCK_SIZE_%d_%d", tc, tr)

  offset_str <- ds$getMetadataItem(band = 1, mdi_name = offset_key, domain = "TIFF")
  size_str   <- ds$getMetadataItem(band = 1, mdi_name = size_key, domain = "TIFF")

  if (!nzchar(offset_str) || !nzchar(size_str)) {
    cat("  Tile", tc, tr, "- no data (empty tile)\n")
    next
  }

  byte_offset <- as.numeric(offset_str)
  byte_length <- as.integer(size_str)

  ## Actual tile dimensions (handle edge tiles)
  tw <- min(block_size[1], src_dim[1] - tc * block_size[1])
  th <- min(block_size[2], src_dim[2] - tr * block_size[2])

  ## Fetch and decode via Rust
  pixels <- rust_fetch_decode_tile(
    url = http_url,
    byte_offset = byte_offset,
    byte_length = byte_length,
    tile_width = as.integer(tw),
    tile_height = as.integer(th),
    predictor = predictor
  )

  ## Place into the buffer
  ## Tile pixel origin relative to buffer (0-based)
  tile_px_col <- (tc - tile_col_range[1]) * block_size[1]
  tile_px_row <- (tr - tile_row_range[1]) * block_size[2]

  ## Copy row by row into the flat buffer
  for (row in seq_len(th)) {
    src_start <- (row - 1) * tw + 1
    src_end   <- src_start + tw - 1

    dst_row <- tile_px_row + row - 1
    dst_start <- dst_row * buf_ncol + tile_px_col + 1
    dst_end   <- dst_start + tw - 1

    if (dst_end <= length(tile_buf)) {
      tile_buf[dst_start:dst_end] <- pixels[src_start:src_end]
    }
  }
}

t1 <- proc.time()
cat("Fetch + decode:", (t1 - t0)[3], "seconds\n")

## ============================================================================
## Step 5: Apply warp map to the tile buffer
## ============================================================================

cat("\n=== Applying warp map ===\n")

## Adjust warp map indices: convert from source-image 1-based to buffer 1-based
wmap_cols_buf <- wmap$src_cols - buf_col_start  # already 1-based from warp map,
wmap_rows_buf <- wmap$src_rows - buf_row_start  # buf starts at buf_col_start (0-based)

t2 <- proc.time()
warped <- rust_apply_warp(
  src_pixels = tile_buf,
  src_ncol = as.integer(buf_ncol),
  src_nrow = as.integer(buf_nrow),
  src_cols = as.integer(wmap_cols_buf),
  src_rows = as.integer(wmap_rows_buf),
  nodata = as.integer(nodata)
)
t3 <- proc.time()
cat("Apply warp:", (t3 - t2)[3], "seconds\n")

## ============================================================================
## Step 6: Compare against GDAL pixel_extract (our proven benchmark)
## ============================================================================

cat("\n=== Comparison: full Rust pipeline vs R pixel_extract ===\n")

dest_xy <- xy_from_cell(dest_dim, dest_ext, seq_len(prod(dest_dim)))
r_pixels <- pixel_extract(ds, bands = 1L, dest_xy, xy_srs = dest_crs, interp = "near")

warped[warped == nodata] <- NA
r_pixels[r_pixels == nodata] <- NA

both_valid <- !is.na(warped) & !is.na(r_pixels)
diff <- abs(warped - r_pixels)

cat("Pixels compared:", sum(both_valid), "\n")
cat("Identical:      ", sum(diff[both_valid] == 0), "\n")
cat("Off by ≤1:      ", sum(diff[both_valid] <= 1), "\n")
cat("Max difference: ", max(diff[both_valid]), "\n")

## ============================================================================
## Step 7: Visual
## ============================================================================

cat("\n=== Visual comparison ===\n")
par(mfrow = c(1, 2))
zlim <- quantile(c(warped, r_pixels), c(0.01, 0.99), na.rm = TRUE)

ximage::ximage(
  matrix(warped, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "Full Rust pipeline\n(warp map + fetch/decode + sample)",
  col = hcl.colors(256), zlim = zlim)
ximage::ximage(
  matrix(r_pixels, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "R pixel_extract\n(benchmark)",
  col = hcl.colors(256), zlim = zlim)

ds$close()

cat("\n=== Done ===\n")
cat("The full pipeline uses NO GDAL for pixel reading.\n")
cat("GDAL is only used for: CRS strings, geotransform, block metadata.\n")
cat("Pixels flow: HTTP range request → Rust DEFLATE decode → Rust warp sample\n")
