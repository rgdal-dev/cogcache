## Diagnostic: what does the ApproxTransformer interpolation look like
## for a single scanline?

library(gdalraster)
library(reproj)
library(vaster)
library(slippymath)

cog_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"
ds <- new(GDALRaster, cog_url)

src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()

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

## Get exact source coords for all pixels
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = dest_crs, dst_gt = dest_gt,
  dst_dim = as.integer(dest_dim)
)

## Look at one scanline (row 128, the middle)
row <- 128
idx <- (row - 1) * 256 + 1:256

exact_src_x <- coords$src_x[idx]
exact_src_y <- coords$src_y[idx]

## What would linear interpolation give? (first to last)
interp_src_x <- seq(exact_src_x[1], exact_src_x[256], length.out = 256)
interp_src_y <- seq(exact_src_y[1], exact_src_y[256], length.out = 256)

## Error
err_x <- exact_src_x - interp_src_x
err_y <- exact_src_y - interp_src_y

cat("Scanline", row, "interpolation error (pixels):\n")
cat("  X: range", range(err_x), "max abs", max(abs(err_x)), "\n")
cat("  Y: range", range(err_y), "max abs", max(abs(err_y)), "\n")
cat("  L1 at midpoint:", abs(err_x[128]) + abs(err_y[128]), "\n")
cat("  Max L1 anywhere:", max(abs(err_x) + abs(err_y)), "\n")

par(mfrow = c(2, 2))
plot(exact_src_x, type = "l", main = "Exact src_x across scanline",
     ylab = "source pixel col", xlab = "dest pixel")
lines(interp_src_x, col = "red", lty = 2)
legend("topleft", c("exact", "interpolated"), col = c("black", "red"), lty = c(1, 2))

plot(err_x, type = "l", main = "X interpolation error (pixels)",
     ylab = "error (src pixels)", xlab = "dest pixel")
abline(h = 0, col = "grey")

plot(exact_src_y, type = "l", main = "Exact src_y across scanline",
     ylab = "source pixel row", xlab = "dest pixel")
lines(interp_src_y, col = "red", lty = 2)

plot(err_y, type = "l", main = "Y interpolation error (pixels)",
     ylab = "error (src pixels)", xlab = "dest pixel")
abline(h = 0, col = "grey")

ds$close()
