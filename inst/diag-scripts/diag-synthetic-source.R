## Diagnostic: eliminate ALL source data ambiguity
##
## Create a tiny GeoTIFF containing just the source neighbourhood,
## warp it with GDAL bilinear, and compare to our Rust bilinear.
## This rules out overviews, block boundaries, IO buffering, etc.
##
## If GDAL and Rust agree on this synthetic source, the bug is in
## how the real COG data reaches the warp kernel in GDAL.

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

## Read source window (same as test_resample-diag.R)
sw <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
src_block <- ds$getBlockSize(band = 1)
src_tw <- src_block[1]; src_th <- src_block[2]
read_xoff  <- (sw$xoff %/% src_tw) * src_tw
read_yoff  <- (sw$yoff %/% src_th) * src_th
read_xend  <- min(ceiling((sw$xoff + sw$xsize) / src_tw) * src_tw, src_dim[1])
read_yend  <- min(ceiling((sw$yoff + sw$ysize) / src_th) * src_th, src_dim[2])
read_xsize <- as.integer(read_xend - read_xoff)
read_ysize <- as.integer(read_yend - read_yoff)

src_pixels <- ds$read(band = 1, xoff = read_xoff, yoff = read_yoff,
                      xsize = read_xsize, ysize = read_ysize,
                      out_xsize = read_xsize, out_ysize = read_ysize)
ds$close()

## Create a small GeoTIFF with exactly these source pixels
## geotransform: shift the original by read_xoff, read_yoff pixels
small_gt <- src_gt
small_gt[1] <- src_gt[1] + read_xoff * src_gt[2]  # xmin shifts
small_gt[4] <- src_gt[4] + read_yoff * src_gt[6]  # ymax shifts

small_tif <- tempfile(fileext = ".tif")
create(dst_filename = small_tif,
        format = "GTiff",
       xsize = read_xsize, ysize = read_ysize,
       nbands = 1L, dataType = "Int16",
       options = c("COMPRESS=NONE"), return_obj = FALSE)

sds <- new(GDALRaster, small_tif, read_only = FALSE)
sds$setGeoTransform(small_gt)
sds$setProjection(src_crs)
sds$write(band = 1, xoff = 0, yoff = 0,
          xsize = read_xsize, ysize = read_ysize,
          data = src_pixels)
sds$close()

## Verify the small file
cat(sprintf("Created source: %d x %d\n", read_xsize, read_ysize))
cat(sprintf("  gt: [%.4f, %.10f, 0, %.4f, 0, %.10f]\n",
            small_gt[1], small_gt[2], small_gt[4], small_gt[6]))

## Now warp the small file with GDAL bilinear
gdal_from_small <- tempfile(fileext = ".tif")
warp(small_tif, gdal_from_small, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds_small <- new(GDALRaster, gdal_from_small)
gdal_small_px <- gds_small$read(band = 1, xoff = 0, yoff = 0,
                                xsize = tile_size, ysize = tile_size,
                                out_xsize = tile_size, out_ysize = tile_size)
gds_small$close()

## And from original COG
gdal_from_cog <- tempfile(fileext = ".tif")
warp(cog_url, gdal_from_cog, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds_cog <- new(GDALRaster, gdal_from_cog)
gdal_cog_px <- gds_cog$read(band = 1, xoff = 0, yoff = 0,
                            xsize = tile_size, ysize = tile_size,
                            out_xsize = tile_size, out_ysize = tile_size)
gds_cog$close()

## Rust bilinear
nodata <- -32768L
our <- rust_warp_resample(
  src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
  src_ncol = read_xsize, src_nrow = read_ysize,
  src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
  nodata = nodata, max_error = 0.0, resample = "bilinear"
)

idx <- 157 * tile_size + 1 + 1

cat("\n=== Rust vs GDAL-from-small-tif ===\n")
both1 <- our != nodata & gdal_small_px != nodata
d1 <- our[both1] - gdal_small_px[both1]
cat(sprintf("  match %d/%d  max_diff %d\n", sum(our == gdal_small_px), tile_size^2, max(abs(d1))))
cat(sprintf("  [1,157]: ours=%d  gdal_small=%d\n", our[idx], gdal_small_px[idx]))

cat("\n=== Rust vs GDAL-from-COG ===\n")
both2 <- our != nodata & gdal_cog_px != nodata
d2 <- our[both2] - gdal_cog_px[both2]
cat(sprintf("  match %d/%d  max_diff %d\n", sum(our == gdal_cog_px), tile_size^2, max(abs(d2))))
cat(sprintf("  [1,157]: ours=%d  gdal_cog=%d\n", our[idx], gdal_cog_px[idx]))

cat("\n=== GDAL-from-small vs GDAL-from-COG ===\n")
both3 <- gdal_small_px != nodata & gdal_cog_px != nodata
d3 <- gdal_small_px[both3] - gdal_cog_px[both3]
cat(sprintf("  match %d/%d  max_diff %d\n",
            sum(gdal_small_px == gdal_cog_px), tile_size^2, max(abs(d3))))
cat(sprintf("  [1,157]: small=%d  cog=%d\n", gdal_small_px[idx], gdal_cog_px[idx]))

cat("\nDone.\n")
