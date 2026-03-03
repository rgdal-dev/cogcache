## Diagnostic: what are the actual coordinate values at [1,157]?
## Run after loading cogcache

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
nodata  <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- -32768L
out_crs <- srs_to_wkt(fiji_lcc)

## Tile [4,3]
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

## Source window
sw <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
cat("Source window from compute_source_window:\n")
cat(sprintf("  xoff=%d yoff=%d xsize=%d ysize=%d\n", sw$xoff, sw$yoff, sw$xsize, sw$ysize))

## Tile-aligned read
src_block <- ds$getBlockSize(band = 1)
src_tw <- src_block[1]; src_th <- src_block[2]
read_xoff  <- (sw$xoff %/% src_tw) * src_tw
read_yoff  <- (sw$yoff %/% src_th) * src_th
read_xend  <- min(ceiling((sw$xoff + sw$xsize) / src_tw) * src_tw, src_dim[1])
read_yend  <- min(ceiling((sw$yoff + sw$ysize) / src_th) * src_th, src_dim[2])
read_xsize <- as.integer(read_xend - read_xoff)
read_ysize <- as.integer(read_yend - read_yoff)
cat(sprintf("  read_xoff=%d read_yoff=%d read_xsize=%d read_ysize=%d\n",
            read_xoff, read_yoff, read_xsize, read_ysize))
cat(sprintf("  left padding = sw$xoff - read_xoff = %d\n", sw$xoff - read_xoff))
cat(sprintf("  top padding  = sw$yoff - read_yoff = %d\n", sw$yoff - read_yoff))

## Transform the problem pixel [col=1, row=157] to see where it lands
## Destination pixel centre: (1 + 0.5, 157 + 0.5) in tile coords
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size)
)

## Problem pixel is at linear index (157 * 256 + 1) + 1 (1-based R)
idx <- 157 * tile_size + 1 + 1  # row 157, col 1, 1-based
src_x_px <- coords$src_x[idx]
src_y_px <- coords$src_y[idx]
cat(sprintf("\nProblem pixel [col=1, row=157]:\n"))
cat(sprintf("  Full-image src coords: x=%.4f y=%.4f\n", src_x_px, src_y_px))
cat(sprintf("  Buffer-relative (read_xoff): buf_x=%.4f buf_y=%.4f\n",
            src_x_px - read_xoff, src_y_px - read_yoff))
cat(sprintf("  Buffer-relative (sw$xoff):   buf_x=%.4f buf_y=%.4f\n",
            src_x_px - sw$xoff, src_y_px - sw$yoff))

## Bilinear kernel needs: ix = floor(buf_x - 0.5), must be >= 0
buf_x <- src_x_px - read_xoff
ix <- floor(buf_x - 0.5)
cat(sprintf("  Bilinear ix = floor(%.4f - 0.5) = %d  (needs ix >= 0 and ix+1 < %d)\n",
            buf_x, ix, read_xsize))

## Also check a few neighbours
cat("\nNearby destination pixels [col, row] -> buf_x:\n")
for (col in 0:3) {
  idx_i <- 157 * tile_size + col + 1
  bx <- coords$src_x[idx_i] - read_xoff
  by <- coords$src_y[idx_i] - read_yoff
  cat(sprintf("  [%d, 157] -> buf_x=%.4f  ix=%d\n", col, bx, floor(bx - 0.5)))
}

ds$close()
cat("\nDone.\n")
