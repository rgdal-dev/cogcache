## Test suite for Rust warp pipeline
##
## Baseline: gdalwarp -et 0 -ovr NONE -r near
## (exact transform, full-resolution source, nearest neighbour)
##
## RESULT: 65536 / 65536 bit-identical match achieved.

library(gdalraster)
library(reproj)
library(vaster)
library(slippymath)

## --- Setup ---

cog_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"

ds <- new(GDALRaster, cog_url)

src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())

tile_extent_3857 <- function(z, x, y) {
  n <- 2^z
  lon <- c(x, x + 1) / n * 360 - 180
  lat <- atan(sinh(pi * (1 - 2 * c(y, y + 1) / n))) * 180 / pi
  reproj_extent(c(lon, rev(lat)), "EPSG:3857", source = "EPSG:4326")
}

llex <- reproj_extent(ds$bbox()[c(1, 3, 2, 4)], "EPSG:4326", source = src_crs)
z <- 10L
tile <- lonlat_to_tilenum(mean(llex[1:2]), mean(llex[3:4]), z)

dest_ext <- tile_extent_3857(z, tile$x, tile$y)
dest_dim <- c(256L, 256L)
dest_crs <- srs_to_wkt("EPSG:3857")
dest_gt  <- extent_dim_to_gt(dest_ext, dest_dim)

## --- Compute source window ---
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = dest_crs, dst_gt = dest_gt,
  dst_dim = as.integer(dest_dim)
)

valid <- !is.na(coords$src_x)
src_col_range <- range(as.integer(coords$src_x[valid]))
src_row_range <- range(as.integer(coords$src_y[valid]))

src_xoff  <- src_col_range[1]
src_yoff  <- src_row_range[1]
src_xsize <- src_col_range[2] - src_col_range[1] + 1L
src_ysize <- src_row_range[2] - src_row_range[1] + 1L

src_window <- ds$read(band = 1,
  xoff = src_xoff, yoff = src_yoff,
  xsize = src_xsize, ysize = src_ysize,
  out_xsize = src_xsize, out_ysize = src_ysize)

nodata <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- 0L

## --- GDAL baseline: exact, full-res, nearest ---
tf <- tempfile(tmpdir = "/vsimem", fileext = ".vrt")
warp(ds, tf, dest_crs,
     cl_arg = c("-ts", dest_dim, "-te", dest_ext[c(1, 3, 2, 4)],
                "-r", "near", "-et", "0", "-ovr", "NONE"))
wds <- new(GDALRaster, tf)
gdal_pixels <- read_ds(wds)
wds$close()

## --- Test 1: Exact warp ---
cat("=== Test 1: Exact warp vs GDAL (-et 0 -ovr NONE) ===\n")
t0 <- proc.time()
exact_pixels <- rust_warp_scanline(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = dest_crs, dst_gt = dest_gt,
  dst_dim = as.integer(dest_dim),
  src_pixels = src_window,
  src_ncol = as.integer(src_xsize), src_nrow = as.integer(src_ysize),
  src_col_off = as.integer(src_xoff), src_row_off = as.integer(src_yoff),
  nodata = as.integer(nodata)
)
t_exact <- (proc.time() - t0)[3]
cat(sprintf("  Match: %d / %d  (time: %.3fs)\n\n",
            sum(exact_pixels == gdal_pixels), prod(dest_dim), t_exact))

## --- Test 2: ApproxTransformer at different thresholds ---
cat("=== Test 2: ApproxTransformer thresholds ===\n\n")
cat(sprintf("max_error | time(s)  | vs exact: identical | vs GDAL: identical\n"))
cat(sprintf("----------|----------|--------------------|-----------------\n"))

for (me in c(0.0, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0)) {
  t0 <- proc.time()
  approx_px <- rust_warp_approx(
    src_crs = src_crs, src_gt = src_gt,
    dst_crs = dest_crs, dst_gt = dest_gt,
    dst_dim = as.integer(dest_dim),
    src_pixels = src_window,
    src_ncol = as.integer(src_xsize), src_nrow = as.integer(src_ysize),
    src_col_off = as.integer(src_xoff), src_row_off = as.integer(src_yoff),
    nodata = as.integer(nodata),
    max_error = me
  )
  t_approx <- (proc.time() - t0)[3]

  vs_exact <- sum(approx_px == exact_pixels)
  vs_gdal  <- sum(approx_px == gdal_pixels)

  cat(sprintf("  %6.4f  |  %6.4f  | %5d / %d     | %5d / %d\n",
              me, t_approx, vs_exact, prod(dest_dim), vs_gdal, prod(dest_dim)))
}

cat(sprintf("\nExact warp time: %.4fs\n", t_exact))

ds$close()
