## Test GenImgProjTransformer implementation
##
## Three comparisons:
##   1. rust_gen_img_proj_transform vs R pixel coords (are the fractional
##      source pixel coordinates correct?)
##   2. rust_warp_scanline vs pixel_extract (does per-scanline warp match?)
##   3. rust_warp_scanline vs rust_warp_map + rust_apply_warp (consistency)

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

## ============================================================================
## Test 1: GenImgProjTransformer - check fractional source pixel coords
## ============================================================================

cat("=== Test 1: GenImgProjTransformer pixel coords ===\n")

t0 <- proc.time()
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs,
  src_gt  = src_gt,
  dst_crs = dest_crs,
  dst_gt  = dest_gt,
  dst_dim = as.integer(dest_dim)
)
t1 <- proc.time()
cat("Transform time:", (t1 - t0)[3], "seconds\n")

## Compare against R: compute dest pixel centres → transform → src pixel coords
dest_xy <- xy_from_cell(dest_dim, dest_ext, seq_len(prod(dest_dim)))
src_xy  <- reproj_xy(dest_xy, target = src_crs, source = dest_crs)

## Convert src geo coords to fractional pixel coords using the inverse gt
## (same as what GenImgProjTransformer does internally)
inv_gt <- c(-src_gt[1], 1, 0, -src_gt[4], 0, 1) / c(src_gt[2], src_gt[2], 1, src_gt[6], 1, src_gt[6])
## Actually let's just use col_from_x / row_from_y from vaster thinking:
r_src_col <- (src_xy[, 1] - src_gt[1]) / src_gt[2] - 0.5
r_src_row <- (src_xy[, 2] - src_gt[4]) / src_gt[6] - 0.5

## Compare
col_diff <- abs(coords$src_x - r_src_col)
row_diff <- abs(coords$src_y - r_src_row)

cat("Src col: max diff =", max(col_diff, na.rm = TRUE),
    "mean diff =", mean(col_diff, na.rm = TRUE), "\n")
cat("Src row: max diff =", max(row_diff, na.rm = TRUE),
    "mean diff =", mean(row_diff, na.rm = TRUE), "\n")

## ============================================================================
## Test 2: rust_warp_scanline vs pixel_extract
## ============================================================================

cat("\n=== Test 2: rust_warp_scanline vs pixel_extract ===\n")

## Get the source window
valid <- !is.na(coords$src_x)
src_col_range <- range(as.integer(coords$src_x[valid]))
src_row_range <- range(as.integer(coords$src_y[valid]))

## Read source window via GDAL
src_xoff  <- src_col_range[1]
src_yoff  <- src_row_range[1]
src_xsize <- src_col_range[2] - src_col_range[1] + 1L
src_ysize <- src_row_range[2] - src_row_range[1] + 1L

cat("Source window:", src_xoff, src_yoff, src_xsize, "x", src_ysize, "\n")

src_window <- ds$read(band = 1,
  xoff = src_xoff, yoff = src_yoff,
  xsize = src_xsize, ysize = src_ysize,
  out_xsize = src_xsize, out_ysize = src_ysize)

nodata <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- 0L

## Warp using the per-scanline kernel
t0 <- proc.time()
warped_scanline <- rust_warp_scanline(
  src_crs   = src_crs,
  src_gt    = src_gt,
  dst_crs   = dest_crs,
  dst_gt    = dest_gt,
  dst_dim   = as.integer(dest_dim),
  src_pixels = src_window,
  src_ncol  = as.integer(src_xsize),
  src_nrow  = as.integer(src_ysize),
  src_col_off = as.integer(src_xoff),
  src_row_off = as.integer(src_yoff),
  nodata    = as.integer(nodata)
)
t1 <- proc.time()
cat("Scanline warp time:", (t1 - t0)[3], "seconds\n")

## Benchmark: pixel_extract
t0 <- proc.time()
r_pixels <- pixel_extract(ds, bands = 1L, dest_xy, xy_srs = dest_crs, interp = "near")
t1 <- proc.time()
cat("pixel_extract time:", (t1 - t0)[3], "seconds\n")

## Compare
warped_scanline[warped_scanline == nodata] <- NA
r_pixels[r_pixels == nodata] <- NA
both_valid <- !is.na(warped_scanline) & !is.na(r_pixels)
diff <- abs(warped_scanline - r_pixels)

cat("Pixels compared:", sum(both_valid), "\n")
cat("Identical:      ", sum(diff[both_valid] == 0), "\n")
cat("Off by ≤1:      ", sum(diff[both_valid] <= 1), "\n")
cat("Max difference: ", max(diff[both_valid]), "\n")

## ============================================================================
## Test 3: Consistency with legacy rust_warp_map interface
## ============================================================================

cat("\n=== Test 3: rust_warp_scanline vs legacy rust_warp_map ===\n")

wmap <- rust_warp_map(
  src_crs = src_crs,
  src_gt  = src_gt,
  src_dim = as.integer(src_dim),
  dst_crs = dest_crs,
  dst_gt  = dest_gt,
  dst_dim = as.integer(dest_dim)
)

wmap_cols_buf <- wmap$src_cols - src_xoff
wmap_rows_buf <- wmap$src_rows - src_yoff

legacy_warped <- rust_apply_warp(
  src_pixels = src_window,
  src_ncol = as.integer(src_xsize),
  src_nrow = as.integer(src_ysize),
  src_cols = as.integer(wmap_cols_buf),
  src_rows = as.integer(wmap_rows_buf),
  nodata = as.integer(nodata)
)

legacy_warped[legacy_warped == nodata] <- NA
both_valid2 <- !is.na(warped_scanline) & !is.na(legacy_warped)
diff2 <- abs(warped_scanline - legacy_warped)

cat("Pixels compared:", sum(both_valid2), "\n")
cat("Identical:      ", sum(diff2[both_valid2] == 0), "\n")
cat("Max difference: ", max(diff2[both_valid2]), "\n")

## ============================================================================
## Test 4: Compare against GDAL warp (for reference)
## ============================================================================

cat("\n=== Test 4: rust_warp_scanline vs GDAL warp ===\n")

tf <- tempfile(tmpdir = "/vsimem", fileext = ".vrt")
warp(ds, tf, dest_crs,
     cl_arg = c("-ts", dest_dim,
                "-te", dest_ext[c(1, 3, 2, 4)],
                "-r", "near"))
warp_ds <- new(GDALRaster, tf)
gdal_pixels <- read_ds(warp_ds)
warp_ds$close()

gdal_pixels[gdal_pixels == nodata] <- NA
both_valid3 <- !is.na(warped_scanline) & !is.na(gdal_pixels)
diff3 <- abs(warped_scanline - gdal_pixels)

cat("Pixels compared:", sum(both_valid3), "\n")
cat("Identical:      ", sum(diff3[both_valid3] == 0), "\n")
cat("Off by ≤1:      ", sum(diff3[both_valid3] <= 1), "\n")
cat("Max difference: ", max(diff3[both_valid3]), "\n")

## ============================================================================
## Visual
## ============================================================================

cat("\n=== Visual ===\n")
par(mfrow = c(1, 3))
zlim <- quantile(c(warped_scanline, r_pixels, gdal_pixels),
                 c(0.01, 0.99), na.rm = TRUE)

ximage::ximage(
  matrix(warped_scanline, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "Rust GenImgProj\n(per-scanline)",
  col = hcl.colors(256), zlim = zlim)
ximage::ximage(
  matrix(r_pixels, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "R pixel_extract\n(benchmark)",
  col = hcl.colors(256), zlim = zlim)
ximage::ximage(
  matrix(gdal_pixels, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "GDAL warp",
  col = hcl.colors(256), zlim = zlim)

ds$close()
