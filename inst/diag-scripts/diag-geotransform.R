## Diagnostic: does GDAL's output geotransform match ours exactly?
##
## If GDAL rounds or adjusts the geotransform from -te/-ts, then
## "pixel [1,157]" means a different location in GDAL's output vs ours.

library(gdalraster)
library(vaster)

cog_url <- "/vsicurl/https://data.source.coop/alexgleith/gebco-2024/GEBCO_2024.tif"

fiji_lcc <- paste(
  "+proj=lcc +lat_0=-18 +lon_0=178",
  "+lat_1=-10 +lat_2=-25",
  "+datum=WGS84 +units=m +no_defs"
)

ds <- new(GDALRaster, cog_url)
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
out_crs <- srs_to_wkt(fiji_lcc)
ds$close()

tile_size <- 256L
pixel_res <- 2000
out_xmin <- -2000000; out_ymin <- -1500000
tx <- 4; ty <- 3
dst_col_off <- (tx - 1L) * tile_size
dst_row_off <- (ty - 1L) * tile_size
tile_xmin <- out_xmin + dst_col_off * pixel_res
tile_xmax <- tile_xmin + tile_size * pixel_res
tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) - dst_row_off * pixel_res
tile_ymin <- tile_ymax - tile_size * pixel_res
tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
tile_gt   <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))

cat("Our tile_gt:\n")
cat(sprintf("  [%.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
            tile_gt[1], tile_gt[2], tile_gt[3], tile_gt[4], tile_gt[5], tile_gt[6]))
cat(sprintf("  extent: xmin=%.4f xmax=%.4f ymin=%.4f ymax=%.4f\n",
            tile_xmin, tile_xmax, tile_ymin, tile_ymax))

## Create GDAL warp output and read its geotransform
gdal_file <- tempfile(fileext = ".tif")
warp(cog_url, gdal_file, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds <- new(GDALRaster, gdal_file)
gdal_gt <- gds$getGeoTransform()
gdal_xs <- gds$getRasterXSize()
gdal_ys <- gds$getRasterYSize()
gds$close()

cat("\nGDAL output geotransform:\n")
cat(sprintf("  [%.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
            gdal_gt[1], gdal_gt[2], gdal_gt[3], gdal_gt[4], gdal_gt[5], gdal_gt[6]))
cat(sprintf("  dims: %d x %d\n", gdal_xs, gdal_ys))

cat("\nDifference (ours - gdal):\n")
dgt <- tile_gt - gdal_gt
cat(sprintf("  [%e, %e, %e, %e, %e, %e]\n",
            dgt[1], dgt[2], dgt[3], dgt[4], dgt[5], dgt[6]))

## What does this mean for pixel [1,157]?
our_geo_x <- tile_gt[1] + 1.5 * tile_gt[2] + 157.5 * tile_gt[3]
our_geo_y <- tile_gt[4] + 1.5 * tile_gt[5] + 157.5 * tile_gt[6]
gdal_geo_x <- gdal_gt[1] + 1.5 * gdal_gt[2] + 157.5 * gdal_gt[3]
gdal_geo_y <- gdal_gt[4] + 1.5 * gdal_gt[5] + 157.5 * gdal_gt[6]

cat(sprintf("\nPixel [1,157] centre in projected coords:\n"))
cat(sprintf("  Ours: (%.6f, %.6f)\n", our_geo_x, our_geo_y))
cat(sprintf("  GDAL: (%.6f, %.6f)\n", gdal_geo_x, gdal_geo_y))
cat(sprintf("  Diff: (%.6f, %.6f) metres\n", our_geo_x - gdal_geo_x, our_geo_y - gdal_geo_y))

## Also check: what -te does GDAL actually honour?
## The extent GDAL uses = gt[1], gt[1]+xs*gt[2], gt[4]+ys*gt[6], gt[4]
gdal_xmin <- gdal_gt[1]
gdal_xmax <- gdal_gt[1] + gdal_xs * gdal_gt[2]
gdal_ymax <- gdal_gt[4]
gdal_ymin <- gdal_gt[4] + gdal_ys * gdal_gt[6]
cat(sprintf("\nGDAL effective extent: xmin=%.4f xmax=%.4f ymin=%.4f ymax=%.4f\n",
            gdal_xmin, gdal_xmax, gdal_ymin, gdal_ymax))

cat("\nDone.\n")
