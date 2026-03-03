## Diagnostic: does removing nodata from GEBCO change GDAL's bilinear output?
##
## Hypothesis: GEBCO has nodata=-32768 metadata. This makes GDAL use the
## masked/general bilinear path instead of GWKBilinearResampleNoMasks4SampleT.
## Our Rust code uses the NoMasks formula. If removing nodata from GEBCO
## makes GDAL produce values matching our Rust output, the fix is to match
## GDAL's masked bilinear path.

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

## Source window + read
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

nodata <- -32768L

## Rust bilinear (our code, et=0 for exact comparison)
our <- rust_warp_resample(
  src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
  src_ncol = read_xsize, src_nrow = read_ysize,
  src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
  nodata = nodata, max_error = 0.0, resample = "bilinear"
)

## GDAL bilinear with nodata (default GEBCO)
gdal_file1 <- tempfile(fileext = ".tif")
warp(cog_url, gdal_file1, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds1 <- new(GDALRaster, gdal_file1)
gdal_with_nd <- gds1$read(band = 1, xoff = 0, yoff = 0,
                          xsize = tile_size, ysize = tile_size,
                          out_xsize = tile_size, out_ysize = tile_size)
gds1$close()

## GDAL bilinear WITHOUT nodata — create a VRT that strips nodata
vrt_file <- tempfile(fileext = ".vrt")
translate(cog_url, vrt_file, cl_arg = c("-of", "VRT", "-a_nodata", "none"))
gdal_file2 <- tempfile(fileext = ".tif")
warp(vrt_file, gdal_file2, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds2 <- new(GDALRaster, gdal_file2)
gdal_no_nd <- gds2$read(band = 1, xoff = 0, yoff = 0,
                        xsize = tile_size, ysize = tile_size,
                        out_xsize = tile_size, out_ysize = tile_size)
gds2$close()

## Compare
cat("=== Rust vs GDAL-with-nodata (original) ===\n")
both1 <- our != nodata & gdal_with_nd != nodata
d1 <- our[both1] - gdal_with_nd[both1]
cat(sprintf("  match %d/%d  corr %.6f  mean %.2f  max_diff %d\n",
            sum(our == gdal_with_nd), tile_size^2,
            cor(our[both1], gdal_with_nd[both1]),
            mean(d1), max(abs(d1))))

cat("\n=== Rust vs GDAL-no-nodata (NoMasks path) ===\n")
both2 <- our != nodata & gdal_no_nd != nodata
d2 <- our[both2] - gdal_no_nd[both2]
cat(sprintf("  match %d/%d  corr %.6f  mean %.2f  max_diff %d\n",
            sum(our == gdal_no_nd), tile_size^2,
            cor(our[both2], gdal_no_nd[both2]),
            mean(d2), max(abs(d2))))

cat("\n=== GDAL-with-nodata vs GDAL-no-nodata ===\n")
both3 <- gdal_with_nd != nodata & gdal_no_nd != nodata
d3 <- gdal_with_nd[both3] - gdal_no_nd[both3]
cat(sprintf("  match %d/%d  corr %.6f  mean %.2f  max_diff %d\n",
            sum(gdal_with_nd == gdal_no_nd), tile_size^2,
            cor(gdal_with_nd[both3], gdal_no_nd[both3]),
            mean(d3), max(abs(d3))))

## Spot check [1,157]
idx <- 157 * tile_size + 1 + 1
cat(sprintf("\nPixel [1,157]: ours=%d  gdal_nd=%d  gdal_no_nd=%d\n",
            our[idx], gdal_with_nd[idx], gdal_no_nd[idx]))

cat("\nDone.\n")
