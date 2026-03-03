## Test collect_chunk_list: recursive destination subdivision
##
## For tile [5,3] (antimeridian crosser), compute_source_window returns
## fill_ratio ~0.01 and a full-width source window. collect_chunk_list
## should subdivide the destination until each chunk has a compact source.
##
## For tile [4,3] (no crossing), it should return a single chunk.

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

make_tile_ext <- function(tx, ty) {
  dst_col_off <- (tx - 1L) * tile_size
  dst_row_off <- (ty - 1L) * tile_size
  tile_xmin <- out_xmin + dst_col_off * pixel_res
  tile_xmax <- tile_xmin + tile_size * pixel_res
  tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) - dst_row_off * pixel_res
  tile_ymin <- tile_ymax - tile_size * pixel_res
  c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
}

cat("=== Tile [4,3]: normal (no antimeridian crossing) ===\n")
ext4 <- make_tile_ext(4, 3)
gt4 <- extent_dim_to_gt(ext4, c(tile_size, tile_size))

sw4 <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = gt4,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
cat(sprintf("  source window: %d x %d  fill_ratio=%.3f\n",
            sw4$xsize, sw4$ysize, sw4$fill_ratio))

chunks4 <- rust_collect_chunk_list(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = gt4,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L, min_fill_ratio = 0.5, min_dst_size = 8L
)
cat(sprintf("  chunks: %d\n", length(chunks4)))
for (i in seq_along(chunks4)) {
  ch <- chunks4[[i]]
  cat(sprintf("    [%d] dst=(%d,%d,%d,%d) src=(%d,%d,%d,%d) fill=%.3f\n",
              i, ch$dst_xoff, ch$dst_yoff, ch$dst_xsize, ch$dst_ysize,
              ch$src_xoff, ch$src_yoff, ch$src_xsize, ch$src_ysize, ch$fill_ratio))
}

cat("\n=== Tile [5,3]: antimeridian crosser ===\n")
ext5 <- make_tile_ext(5, 3)
gt5 <- extent_dim_to_gt(ext5, c(tile_size, tile_size))

sw5 <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = gt5,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
cat(sprintf("  source window: %d x %d  fill_ratio=%.3f\n",
            sw5$xsize, sw5$ysize, sw5$fill_ratio))

chunks5 <- rust_collect_chunk_list(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = gt5,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L, min_fill_ratio = 0.5, min_dst_size = 8L
)
cat(sprintf("  chunks: %d\n", length(chunks5)))

total_src_pixels <- 0
for (i in seq_along(chunks5)) {
  ch <- chunks5[[i]]
  src_px <- ch$src_xsize * ch$src_ysize
  total_src_pixels <- total_src_pixels + src_px
  src_lon_min <- src_gt[1] + ch$src_xoff * src_gt[2]
  src_lon_max <- src_gt[1] + (ch$src_xoff + ch$src_xsize) * src_gt[2]
  cat(sprintf("    [%d] dst=(%d,%d,%d,%d) src=(%d,%d,%d,%d) fill=%.3f lon=[%.1f,%.1f]\n",
              i, ch$dst_xoff, ch$dst_yoff, ch$dst_xsize, ch$dst_ysize,
              ch$src_xoff, ch$src_yoff, ch$src_xsize, ch$src_ysize,
              ch$fill_ratio, src_lon_min, src_lon_max))
}

full_width_pixels <- sw5$xsize * sw5$ysize
cat(sprintf("\n  Total source pixels: %d (%.1f%% of full-width %d)\n",
            total_src_pixels, 100 * total_src_pixels / full_width_pixels,
            full_width_pixels))

cat("\nDone.\n")
