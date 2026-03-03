## Test rust_compute_source_window against the Fiji LCC demo
##
## Compares the Rust ComputeSourceWindow against the R source_footprint()
## function and checks that tile 5's source blowup is resolved.

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

out_xmin <- -2000000; out_xmax <- 2000000
out_ymin <- -1500000; out_ymax <- 1500000
pixel_res <- 2000
tile_size <- 256L

out_ncol <- as.integer(ceiling((out_xmax - out_xmin) / pixel_res))
out_nrow <- as.integer(ceiling((out_ymax - out_ymin) / pixel_res))
out_extent <- c(out_xmin, out_xmin + out_ncol * pixel_res,
                out_ymin, out_ymin + out_nrow * pixel_res)
out_gt <- extent_dim_to_gt(out_extent, c(out_ncol, out_nrow))
out_crs <- srs_to_wkt(fiji_lcc)

n_tiles_x <- ceiling(out_ncol / tile_size)
n_tiles_y <- ceiling(out_nrow / tile_size)

cat(sprintf("Output: %d x %d, %d x %d tiles\n",
            out_ncol, out_nrow, n_tiles_x, n_tiles_y))

## --- Test every tile ---

cat("\nTesting rust_compute_source_window on all tiles...\n\n")

src_block <- ds$getBlockSize(band = 1)
src_tw <- src_block[1]; src_th <- src_block[2]

for (ty in seq_len(n_tiles_y)) {
  for (tx in seq_len(n_tiles_x)) {
    dst_col_off <- (tx - 1L) * tile_size
    dst_row_off <- (ty - 1L) * tile_size
    dst_ncol <- min(tile_size, out_ncol - dst_col_off)
    dst_nrow <- min(tile_size, out_nrow - dst_row_off)

    tile_xmin <- out_extent[1] + dst_col_off * pixel_res
    tile_xmax <- tile_xmin + dst_ncol * pixel_res
    tile_ymax <- out_extent[4] - dst_row_off * pixel_res
    tile_ymin <- tile_ymax - dst_nrow * pixel_res
    tile_extent <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
    tile_gt <- extent_dim_to_gt(tile_extent, c(dst_ncol, dst_nrow))

    sw <- rust_compute_source_window(
      src_crs  = src_crs,
      src_gt   = src_gt,
      src_dim  = as.integer(src_dim),
      dst_crs  = out_crs,
      dst_gt   = tile_gt,
      dst_off  = c(0L, 0L),
      dst_size = as.integer(c(dst_ncol, dst_nrow)),
      resample_padding = 0L
    )


    if (is.null(sw)) {
      cat(sprintf("  [%d,%d] NULL (no valid transforms)\n", tx, ty))
      next
    }

    ## How many source COG blocks does this span?
    n_blk_x <- ceiling(sw$xsize / src_tw)
    n_blk_y <- ceiling(sw$ysize / src_th)
    n_blks  <- n_blk_x * n_blk_y
    src_mpx <- sw$xsize * sw$ysize / 1e6

    cat(sprintf("  [%d,%d] src off=(%d,%d) size=(%d,%d) = %.1f Mpx, %d blks, fill=%.2f, failed=%d/%d\n",
                tx, ty, sw$xoff, sw$yoff, sw$xsize, sw$ysize,
                src_mpx, n_blks, sw$fill_ratio, sw$n_failed, sw$n_samples))
  }
}

ds$close()
cat("\nDone.\n")
