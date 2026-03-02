## Test rust_warp_map against R prototype
##
## Run after: rextendr::document() and devtools::load_all()

library(gdalraster)
library(reproj)
library(vaster)
library(slippymath)

## --- Setup: same grid specs as the qmd ---

cog_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"

ds <- new(GDALRaster, cog_url)

src_ext <- ds$bbox()[c(1, 3, 2, 4)]
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())

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

## --- R warp (our benchmark) ---

cat("=== R warp via pixel_extract ===\n")
dest_xy <- xy_from_cell(dest_dim, dest_ext, seq_len(prod(dest_dim)))
r_pixels <- pixel_extract(ds, bands = 1L, dest_xy, xy_srs = dest_crs, interp = "near")

## --- Rust warp map ---

cat("\n=== Rust warp map ===\n")
wmap <- rust_warp_map(
  src_crs = src_crs,
  src_gt  = src_gt,
  src_dim = as.integer(src_dim),
  dst_crs = dest_crs,
  dst_gt  = dest_gt,
  dst_dim = as.integer(dest_dim)
)

cat("Mapped pixels:", sum(!is.na(wmap$src_cols)), "/", length(wmap$src_cols), "\n")

## --- Read source window and apply warp ---

## Find the bounding box of source pixels we need
valid <- !is.na(wmap$src_cols)
src_xoff  <- min(wmap$src_cols[valid]) - 1L  # 0-based for gdalraster
src_yoff  <- min(wmap$src_rows[valid]) - 1L
src_xsize <- max(wmap$src_cols[valid]) - src_xoff
src_ysize <- max(wmap$src_rows[valid]) - src_yoff

cat("Source window:", src_xoff, src_yoff, src_xsize, "x", src_ysize, "\n")

src_window <- ds$read(band = 1,
  xoff = src_xoff, yoff = src_yoff,
  xsize = src_xsize, ysize = src_ysize,
  out_xsize = src_xsize, out_ysize = src_ysize)

## Apply warp map (adjust indices to window-relative)
nodata <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- 0L

## Adjust 1-based src_cols/src_rows to be relative to the window
wmap_cols_adj <- wmap$src_cols - src_xoff
wmap_rows_adj <- wmap$src_rows - src_yoff

rust_pixels <- rust_apply_warp(
  src_pixels = src_window,
  src_ncol = src_xsize,
  src_nrow = src_ysize,
  src_cols = wmap_cols_adj,
  src_rows = wmap_rows_adj,
  nodata = as.integer(nodata)
)

## --- Compare Rust warp vs R pixel_extract ---

cat("\n=== Comparison: Rust warp vs R pixel_extract ===\n")

## Mask nodata
rust_pixels[rust_pixels == nodata] <- NA
r_pixels[r_pixels == nodata] <- NA

both_valid <- !is.na(rust_pixels) & !is.na(r_pixels)
diff <- abs(rust_pixels - r_pixels)

cat("Pixels compared:", sum(both_valid), "\n")
cat("Identical:      ", sum(diff[both_valid] == 0), "\n")
cat("Off by ≤1:      ", sum(diff[both_valid] <= 1), "\n")
cat("Max difference: ", max(diff[both_valid]), "\n")

## --- Also compare against GDAL warp ---

cat("\n=== Comparison: Rust warp vs GDAL warp ===\n")
tf <- tempfile(tmpdir = "/vsimem", fileext = ".vrt")
warp(ds, tf, dest_crs,
     cl_arg = c("-ts", dest_dim,
                "-te", dest_ext[c(1, 3, 2, 4)],
                "-r", "near"))
warp_ds <- new(GDALRaster, tf)
gdal_pixels <- read_ds(warp_ds)
warp_ds$close()

gdal_pixels[gdal_pixels == nodata] <- NA
both_valid2 <- !is.na(rust_pixels) & !is.na(gdal_pixels)
diff2 <- abs(rust_pixels - gdal_pixels)

cat("Pixels compared:", sum(both_valid2), "\n")
cat("Identical:      ", sum(diff2[both_valid2] == 0), "\n")
cat("Off by ≤1:      ", sum(diff2[both_valid2] <= 1), "\n")
cat("Max difference: ", max(diff2[both_valid2]), "\n")

## --- Visual ---

cat("\n=== Visual comparison ===\n")
par(mfrow = c(1, 3))
zlim <- quantile(c(rust_pixels, r_pixels, gdal_pixels), c(0.01, 0.99), na.rm = TRUE)

ximage::ximage(matrix(rust_pixels, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "Rust warp", col = hcl.colors(256), zlim = zlim)
ximage::ximage(matrix(r_pixels, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "R pixel_extract", col = hcl.colors(256), zlim = zlim)
ximage::ximage(matrix(gdal_pixels, nrow = dest_dim[2], ncol = dest_dim[1], byrow = TRUE),
  dest_ext, main = "GDAL warp", col = hcl.colors(256), zlim = zlim)

ds$close()
