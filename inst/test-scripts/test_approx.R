## Test ApproxTransformer
##
## Compare:
##   1. rust_warp_approx vs rust_warp_scanline (exact) — pixel differences
##   2. rust_warp_approx vs GDAL warp — do we converge?
##   3. Timing at different max_error values
##   4. Visual comparison

library(gdalraster)
library(reproj)
library(vaster)
library(slippymath)

## --- Setup (same as before) ---

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

## Read source window (same as before)
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

## --- Exact baseline ---

cat("=== Exact warp (baseline) ===\n")
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
cat("Exact time:", t_exact, "seconds\n\n")

## --- GDAL warp baseline ---

tf <- tempfile(tmpdir = "/vsimem", fileext = ".vrt")
warp(ds, tf, dest_crs,
     cl_arg = c("-ts", dest_dim, "-te", dest_ext[c(1, 3, 2, 4)], "-r", "near"))
warp_ds <- new(GDALRaster, tf)
gdal_pixels <- read_ds(warp_ds)
warp_ds$close()

## --- Test at different error thresholds ---

cat("=== ApproxTransformer at different error thresholds ===\n\n")

thresholds <- c(0.0, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0)

results <- data.frame(
  max_error = numeric(),
  time_sec = numeric(),
  vs_exact_identical = integer(),
  vs_exact_max_diff = integer(),
  vs_gdal_identical = integer(),
  vs_gdal_max_diff = integer(),
  stringsAsFactors = FALSE
)

approx_pixels_list <- list()

for (me in thresholds) {
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

  approx_pixels_list[[as.character(me)]] <- approx_px

  ## Compare vs exact
  a <- approx_px; a[a == nodata] <- NA
  e <- exact_pixels; e[e == nodata] <- NA
  g <- gdal_pixels; g[g == nodata] <- NA

  bv_exact <- !is.na(a) & !is.na(e)
  bv_gdal  <- !is.na(a) & !is.na(g)

  diff_exact <- abs(a - e)
  diff_gdal  <- abs(a - g)

  results <- rbind(results, data.frame(
    max_error = me,
    time_sec = round(t_approx, 4),
    vs_exact_identical = sum(diff_exact[bv_exact] == 0),
    vs_exact_max_diff = max(diff_exact[bv_exact]),
    vs_gdal_identical = sum(diff_gdal[bv_gdal] == 0),
    vs_gdal_max_diff = max(diff_gdal[bv_gdal])
  ))
}

## Print results table
cat("max_error | time(s)  | vs exact: identical / max_diff | vs GDAL: identical / max_diff\n")
cat("----------|----------|-------------------------------|-------------------------------\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %6.4f  | %7.4f  | %5d / 65536  max_diff=%2d    | %5d / 65536  max_diff=%2d\n",
    r$max_error, r$time_sec,
    r$vs_exact_identical, r$vs_exact_max_diff,
    r$vs_gdal_identical, r$vs_gdal_max_diff))
}

cat(sprintf("\nExact warp time: %.4f seconds\n", t_exact))
cat(sprintf("Speedup at max_error=0.125: %.1fx\n",
  t_exact / results$time_sec[results$max_error == 0.125]))

## --- Visual comparison ---

cat("\n=== Visual ===\n")
par(mfrow = c(2, 2))
zlim <- quantile(c(exact_pixels, gdal_pixels), c(0.01, 0.99), na.rm = TRUE)

to_mat <- function(v) matrix(v, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE)

ximage::ximage(to_mat(approx_pixels_list[["0.125"]]), dest_ext,
  main = "Approx (max_error=0.125)", col = hcl.colors(256), zlim = zlim)
ximage::ximage(to_mat(gdal_pixels), dest_ext,
  main = "GDAL warp", col = hcl.colors(256), zlim = zlim)

## Difference maps against GDAL
diff_approx_gdal <- approx_pixels_list[["0.125"]] - gdal_pixels
diff_approx_gdal[abs(diff_approx_gdal) > 1000] <- NA
diff_exact_gdal <- exact_pixels - gdal_pixels
diff_exact_gdal[abs(diff_exact_gdal) > 1000] <- NA

dlim <- range(c(diff_approx_gdal, diff_exact_gdal), na.rm = TRUE)

ximage::ximage(to_mat(diff_approx_gdal), dest_ext,
  main = "Approx - GDAL", col = hcl.colors(256), zlim = dlim)
ximage::ximage(to_mat(diff_exact_gdal), dest_ext,
  main = "Exact - GDAL", col = hcl.colors(256), zlim = dlim)

ds$close()
