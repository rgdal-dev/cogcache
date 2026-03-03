## Verification: does our bilinear match GDAL when src/dst ratio is ~1:1?
##
## The 4.5:1 downsampling causes GDAL to scale the bilinear kernel.
## At 1:1 (or upsampling), GDAL should use the standard 2x2 bilinear
## and our results should match closely.

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

## GEBCO pixel size is ~1/240 degree ≈ 463m at equator, ~430m at -18°
## Use 500m output resolution → ~1:1 ratio
pixel_res_1to1 <- 500  # metres

tile_size <- 256L
out_xmin <- -2000000; out_ymin <- -1500000
tx <- 4; ty <- 3
dst_col_off <- (tx - 1L) * tile_size
dst_row_off <- (ty - 1L) * tile_size
tile_xmin <- out_xmin + dst_col_off * pixel_res_1to1
tile_xmax <- tile_xmin + tile_size * pixel_res_1to1
tile_ymax <- (out_ymin + ceiling(750000 / pixel_res_1to1) * pixel_res_1to1) - dst_row_off * pixel_res_1to1
tile_ymin <- tile_ymax - tile_size * pixel_res_1to1
tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
tile_gt   <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))

## Source window
sw <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
cat(sprintf("Source window: %d x %d for %d x %d destination\n",
            sw$xsize, sw$ysize, tile_size, tile_size))
cat(sprintf("Ratio: %.2f x %.2f\n", sw$xsize / tile_size, sw$ysize / tile_size))

## Read source
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

## Rust bilinear
our <- rust_warp_resample(
  src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
  src_ncol = read_xsize, src_nrow = read_ysize,
  src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
  nodata = nodata, max_error = 0.0, resample = "bilinear"
)

## GDAL bilinear
gdal_file <- tempfile(fileext = ".tif")
warp(cog_url, gdal_file, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds <- new(GDALRaster, gdal_file)
gdal_px <- gds$read(band = 1, xoff = 0, yoff = 0,
                    xsize = tile_size, ysize = tile_size,
                    out_xsize = tile_size, out_ysize = tile_size)
gds$close()

both <- our != nodata & gdal_px != nodata
diffs <- our[both] - gdal_px[both]
n_match <- sum(our == gdal_px)
cat(sprintf("\n=== 1:1 ratio bilinear comparison ===\n"))
cat(sprintf("  match: %d/%d (%.1f%%)\n", n_match, tile_size^2, 100*n_match/tile_size^2))
cat(sprintf("  corr: %.8f\n", cor(our[both], gdal_px[both])))
cat(sprintf("  mean diff: %.4f\n", mean(diffs)))
cat(sprintf("  max |diff|: %d\n", max(abs(diffs))))
cat(sprintf("  |diff| > 10: %d\n", sum(abs(diffs) > 10)))
cat(sprintf("  |diff| > 1: %d\n", sum(abs(diffs) > 1)))

cat("\nDone.\n")
