## Test plan_warp_reads with split-read detection
library(gdalraster)
library(vaster)

cog_url <- "/vsicurl/https://data.source.coop/alexgleith/gebco-2024/GEBCO_2024.tif"

fiji_lcc <- paste(
  "+proj=lcc +lat_0=-18 +lon_0=178",
  "+lat_1=-10 +lat_2=-25",
  "+datum=WGS84 +units=m +no_defs"
)

ds <- new(GDALRaster, cog_url, TRUE)
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
out_crs <- srs_to_wkt(fiji_lcc)
ds$close()

tile_size <- 256L
pixel_res <- 2000
out_xmin <- -2000000; out_ymin <- -1500000

make_tile_gt <- function(tx, ty) {
  dst_col_off <- (tx - 1L) * tile_size
  dst_row_off <- (ty - 1L) * tile_size
  tile_xmin <- out_xmin + dst_col_off * pixel_res
  tile_xmax <- tile_xmin + tile_size * pixel_res
  tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) -
    dst_row_off * pixel_res
  tile_ymin <- tile_ymax - tile_size * pixel_res
  extent_dim_to_gt(c(tile_xmin, tile_xmax, tile_ymin, tile_ymax),
                   c(tile_size, tile_size))
}

cat("=== Tile [5,3]: antimeridian crosser ===\n\n")
gt5 <- make_tile_gt(5, 3)

plan <- plan_warp_reads(
  src_crs = src_crs, src_gt = src_gt, src_dim = src_dim,
  dst_crs = out_crs, dst_gt = gt5,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)

total_src_pixels <- 0
for (i in seq_along(plan)) {
  p <- plan[[i]]
  cat(sprintf("Chunk %d: dst=(%d,%d,%d,%d) split=%s\n",
              i, p$dst_off[1], p$dst_off[2],
              p$dst_size[1], p$dst_size[2], p$is_split))
  for (j in seq_along(p$src_reads)) {
    sr <- p$src_reads[[j]]
    px <- sr$xsize * sr$ysize
    total_src_pixels <- total_src_pixels + px
    lon_min <- src_gt[1] + sr$xoff * src_gt[2]
    lon_max <- src_gt[1] + (sr$xoff + sr$xsize) * src_gt[2]
    cat(sprintf("  read[%d]: src=(%d,%d,%d,%d) = %s px  lon=[%.1f,%.1f]\n",
                j, sr$xoff, sr$yoff, sr$xsize, sr$ysize,
                format(px, big.mark = ","), lon_min, lon_max))
  }
}

full_width <- 86400L * 1149L
cat(sprintf("\nTotal source pixels: %s (%.1f%% of full-width %s)\n",
            format(total_src_pixels, big.mark = ","),
            100 * total_src_pixels / full_width,
            format(full_width, big.mark = ",")))
