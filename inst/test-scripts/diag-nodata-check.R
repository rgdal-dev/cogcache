## Diagnostic: check for -32767 (actual nodata) in source data
## and test if passing correct nodata to our kernel changes results

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

## Check actual nodata
nd_actual <- ds$getNoDataValue(1)
cat(sprintf("Actual nodata from GDAL: %s\n", nd_actual))

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

## Check for nodata values in source
cat(sprintf("\nSource buffer: %d pixels\n", length(src_pixels)))
cat(sprintf("  min=%d  max=%d\n", min(src_pixels), max(src_pixels)))
cat(sprintf("  count of -32768: %d\n", sum(src_pixels == -32768L)))
cat(sprintf("  count of -32767: %d\n", sum(src_pixels == -32767L)))

## Now test: does the nodata value we pass to our kernel matter?
## Try with -32768 (what test scripts use) vs -32767 (actual) vs NA-equivalent
for (nd in c(-32768L, -32767L, -99999L)) {
  our <- rust_warp_resample(
    src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
    dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
    src_ncol = read_xsize, src_nrow = read_ysize,
    src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
    nodata = nd, max_error = 0.0, resample = "bilinear"
  )
  idx <- 157 * tile_size + 1 + 1
  n_nodata <- sum(our == nd)
  cat(sprintf("  nodata=%d: [1,157]=%d  n_nodata_out=%d\n", nd, our[idx], n_nodata))
}

## KEY TEST: what if there's a signed/unsigned issue in the source read?
## gdalraster reads Int16 as R integer. But our Rust code receives i32.
## Check: are there any values that would differ between Int16 and UInt16
## interpretation?
cat(sprintf("\n  Any values > 32767: %d\n", sum(src_pixels > 32767L)))
cat(sprintf("  Any values < -32768: %d\n", sum(src_pixels < -32768L)))

## Actually, let's check what R thinks the type is
cat(sprintf("  R type: %s\n", typeof(src_pixels)))
cat(sprintf("  Range: [%d, %d]\n", min(src_pixels), max(src_pixels)))

cat("\nDone.\n")
