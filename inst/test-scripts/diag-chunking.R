## Diagnostic: control GDAL's internal chunking and check source window
##
## The -4030 contour is ~2 pixels from our landing point.
## GDAL and our transforms are identical (proven).
## Could GDAL's internal chunking cause source downsampling?
##
## Test: use -wm (warp memory) to force GDAL to process in one chunk

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

idx <- 157 * tile_size + 1 + 1
nodata <- -32768L

## Test 1: default gdalwarp
f1 <- tempfile(fileext = ".tif")
warp(cog_url, f1, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
g1 <- new(GDALRaster, f1)
px1 <- g1$read(band = 1, xoff = 0, yoff = 0,
               xsize = tile_size, ysize = tile_size,
               out_xsize = tile_size, out_ysize = tile_size)
g1$close()

## Test 2: massive warp memory to force single chunk
f2 <- tempfile(fileext = ".tif")
warp(cog_url, f2, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0",
                "-wm", "2000"))
g2 <- new(GDALRaster, f2)
px2 <- g2$read(band = 1, xoff = 0, yoff = 0,
               xsize = tile_size, ysize = tile_size,
               out_xsize = tile_size, out_ysize = tile_size)
g2$close()

## Test 3: use the small synthetic tif (no COG complications)
## but also with huge warp memory
tile_gt <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))
sw <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
ds <- new(GDALRaster, cog_url)
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

## Rust bilinear
our <- rust_warp_resample(
  src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
  src_ncol = read_xsize, src_nrow = read_ysize,
  src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
  nodata = nodata, max_error = 0.0, resample = "bilinear"
)

cat(sprintf("At [1,157]:\n"))
cat(sprintf("  Rust:                   %d\n", our[idx]))
cat(sprintf("  GDAL default:           %d\n", px1[idx]))
cat(sprintf("  GDAL wm=2000:           %d\n", px2[idx]))

cat(sprintf("\nGDAL default vs wm=2000: match=%d/%d  max_diff=%d\n",
            sum(px1 == px2), tile_size^2, max(abs(px1[px1 != nodata & px2 != nodata] - px2[px1 != nodata & px2 != nodata]))))

## Test 4: CRITICAL - what if GDAL is using the source OVERVIEW for
## the source read even with -ovr NONE?
## -ovr NONE means "don't auto-select overview for output resolution"
## but does it affect the SOURCE read?
## Force with explicit overview level 0 (base resolution)
f4 <- tempfile(fileext = ".tif")
## Try --config GDAL_DISABLE_READDIR_ON_OPEN EMPTY to avoid any overview detection
set_config_option("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY")
warp(cog_url, f4, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0",
                "-wm", "2000"))
set_config_option("GDAL_DISABLE_READDIR_ON_OPEN", "")
g4 <- new(GDALRaster, f4)
px4 <- g4$read(band = 1, xoff = 0, yoff = 0,
               xsize = tile_size, ysize = tile_size,
               out_xsize = tile_size, out_ysize = tile_size)
g4$close()
cat(sprintf("  GDAL DISABLE_READDIR:   %d\n", px4[idx]))

## Test 5: what about -wo SOURCE_EXTRA? This adds extra source pixels
f5 <- tempfile(fileext = ".tif")
warp(cog_url, f5, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0",
                "-wo", "SOURCE_EXTRA=100"))
g5 <- new(GDALRaster, f5)
px5 <- g5$read(band = 1, xoff = 0, yoff = 0,
               xsize = tile_size, ysize = tile_size,
               out_xsize = tile_size, out_ysize = tile_size)
g5$close()
cat(sprintf("  GDAL SOURCE_EXTRA=100:  %d\n", px5[idx]))

## Test 6: what if GDAL's internal source window leads to downsampled read?
## Check what GDAL thinks the source window is by using -wo INIT_DEST=NO_DATA
## and checking if output dimensions suggest chunking
cat(sprintf("\nSource window comparison:\n"))
cat(sprintf("  Our sw: xoff=%d yoff=%d xsize=%d ysize=%d\n",
            sw$xoff, sw$yoff, sw$xsize, sw$ysize))
cat(sprintf("  Our read (tile-aligned): %d x %d at (%d, %d)\n",
            read_xsize, read_ysize, read_xoff, read_yoff))
cat(sprintf("  Ratio: src/dst = %.1f x %.1f (>>1 means downsampling)\n",
            sw$xsize / tile_size, sw$ysize / tile_size))

cat("\nDone.\n")
