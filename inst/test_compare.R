## Test: compare Rust decode with R decode
##
## Run this after building the cogcache package:
##   rextendr::document()  # or devtools::load_all()
##
## Uses the same COG and tile data from the qmd prototype.

library(cogcache)
library(gdalraster)
library(httr2)

cog_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"

ds <- new(GDALRaster, cog_url)
block_size <- ds$getBlockSize(band = 1)
dim_xy <- c(ds$getRasterXSize(), ds$getRasterYSize())

## --- R decode function (from the qmd) ----------------------------------------
decode_deflate_tile_r <- function(raw_bytes, tile_width, tile_height,
                                   nbytes = 2L, predictor = 2L) {
  decompressed <- memDecompress(raw_bytes, type = "gzip")
  pix <- readBin(decompressed, what = "integer", n = tile_width * tile_height,
                 size = nbytes, signed = FALSE, endian = "little")
  if (predictor == 1L) return(pix)
  modulus <- 2L^(nbytes * 8L)
  mat <- matrix(pix, nrow = tile_height, ncol = tile_width, byrow = TRUE)
  for (i in seq_len(nrow(mat))) {
    mat[i, ] <- cumsum(mat[i, ]) %% modulus
  }
  as.integer(t(mat))
}

## --- Fetch a tile with real data ---------------------------------------------
tc <- 5L; tr <- 5L
actual_w <- min(block_size[1], dim_xy[1] - tc * block_size[1])
actual_h <- min(block_size[2], dim_xy[2] - tr * block_size[2])

## Get byte offset
offset_str <- ds$getMetadataItem(band = 1,
  mdi_name = sprintf("BLOCK_OFFSET_%d_%d", tc, tr), domain = "TIFF")
size_str <- ds$getMetadataItem(band = 1,
  mdi_name = sprintf("BLOCK_SIZE_%d_%d", tc, tr), domain = "TIFF")
byte_offset <- as.numeric(offset_str)
byte_length <- as.integer(size_str)

## Fetch raw bytes via R
http_url <- sub("^/vsicurl/", "", cog_url)
req <- request(http_url) |>
  req_headers(Range = sprintf("bytes=%d-%d", byte_offset, byte_offset + byte_length - 1))
resp <- req_perform(req)
raw_bytes <- resp_body_raw(resp)

cat("Tile:", tc, tr, "\n")
cat("Dimensions:", actual_w, "x", actual_h, "\n")
cat("Compressed bytes:", length(raw_bytes), "\n\n")

## --- Test 1: Rust decode from raw bytes (same bytes, different decoder) ------
cat("Test 1: rust_decode_tile() vs R decode ... ")
r_pix <- decode_deflate_tile_r(raw_bytes, actual_w, actual_h, predictor = 2L)
rust_pix <- rust_decode_tile(raw_bytes, actual_w, actual_h, predictor = 2L)

if (identical(r_pix, rust_pix)) {
  cat("PASS — identical\n")
} else {
  cat("FAIL\n")
  cat("  R first 10:   ", head(r_pix, 10), "\n")
  cat("  Rust first 10:", head(rust_pix, 10), "\n")
  cat("  Mismatches:   ", sum(r_pix != rust_pix), "of", length(r_pix), "\n")
}

## --- Test 2: Rust fetch+decode (full pipeline in Rust) -----------------------
cat("\nTest 2: rust_fetch_decode_tile() vs R decode ... ")
rust_pix2 <- rust_fetch_decode_tile(
  http_url, byte_offset, byte_length,
  actual_w, actual_h, predictor = 2L
)

if (identical(r_pix, rust_pix2)) {
  cat("PASS — identical\n")
} else {
  cat("FAIL\n")
  cat("  R first 10:   ", head(r_pix, 10), "\n")
  cat("  Rust first 10:", head(rust_pix2, 10), "\n")
  cat("  Mismatches:   ", sum(r_pix != rust_pix2), "of", length(r_pix), "\n")
}

## --- Test 3: Compare with GDAL -----------------------------------------------
cat("\nTest 3: Rust decode vs GDAL vrt://?block= ... ")
block_url <- sprintf("vrt://%s?block=%d,%d", cog_url, tc, tr)
ds_block <- new(GDALRaster, block_url)
gdal_pix <- ds_block$read(band = 1, xoff = 0, yoff = 0,
                           xsize = ds_block$getRasterXSize(),
                           ysize = ds_block$getRasterYSize(),
                           out_xsize = ds_block$getRasterXSize(),
                           out_ysize = ds_block$getRasterYSize())
ds_block$close()

## Apply nodata
nodata <- ds$getNoDataValue(band = 1)
rust_pix_na <- rust_pix
if (!is.na(nodata)) rust_pix_na[rust_pix_na == nodata] <- NA_integer_

if (identical(as.numeric(rust_pix_na), as.numeric(gdal_pix))) {
  cat("PASS — identical to GDAL\n")
} else {
  cat("FAIL\n")
}

## --- Benchmark ---------------------------------------------------------------
cat("\nBenchmark: R decode vs Rust decode (1000 iterations on same tile)\n")
cat("R:    ")
print(system.time(for (i in 1:1000) decode_deflate_tile_r(raw_bytes, actual_w, actual_h)))
cat("Rust: ")
print(system.time(for (i in 1:1000) rust_decode_tile(raw_bytes, actual_w, actual_h)))

ds$close()
