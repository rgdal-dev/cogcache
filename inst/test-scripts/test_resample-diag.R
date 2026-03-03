## Diagnostic: separate ApproxTransformer effect from kernel accuracy
##
## For bilinear, compare:
##   1. Rust bilinear with approx (max_error=0.125) vs GDAL bilinear (default)
##   2. Rust bilinear with exact (max_error=0) vs GDAL bilinear -et 0
##
## If the max_diff drops dramatically with exact transforms,
## the issue is ApproxTransformer coordinate differences amplified
## by interpolation — not a kernel bug.

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

## Source window with bilinear padding
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

cat("=== Bilinear: Approx (default) vs Exact ===\n\n")

for (et in c(0.125, 0.0)) {
  label <- if (et > 0) "approx (et=0.125)" else "exact  (et=0)"

  ## Rust
  our <- rust_warp_resample(
    src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
    dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
    src_ncol = read_xsize, src_nrow = read_ysize,
    src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
    nodata = as.integer(nodata), max_error = et, resample = "bilinear"
  )

  ## GDAL
  gdal_file <- tempfile(fileext = ".tif")
  et_args <- if (et > 0) c() else c("-et", "0")
  warp(ds, gdal_file, out_crs,
       cl_arg = c("-ts", tile_size, tile_size,
                  "-te", tile_ext[c(1, 3, 2, 4)],
                  "-r", "bilinear", "-ovr", "NONE", et_args))
  gds <- new(GDALRaster, gdal_file)
  gdal_px <- gds$read(band = 1, xoff = 0, yoff = 0,
                      xsize = tile_size, ysize = tile_size,
                      out_xsize = tile_size, out_ysize = tile_size)
  gds$close()

  both_valid <- our != nodata & gdal_px != nodata
  n_both <- sum(both_valid)
  n_match <- sum(our == gdal_px)
  diffs <- our[both_valid] - gdal_px[both_valid]
  corr <- if (n_both > 10) cor(our[both_valid], gdal_px[both_valid]) else NA

  cat(sprintf("  %s:  match %5d/%d (%.1f%%)  corr %.6f  mean %.2f  max %d  |diffs|>100: %d\n",
              label, n_match, tile_size^2, 100 * n_match / tile_size^2,
              corr, mean(diffs), max(abs(diffs)),
              sum(abs(diffs) > 100)))

  ## Where are the big diffs?
  if (any(abs(diffs) > 500)) {
    big <- which(abs(our - gdal_px) > 500 & both_valid)
    cat(sprintf("    Big diffs (>500): %d pixels. Examples:\n", length(big)))
    for (i in head(big, 5)) {
      row_i <- (i - 1) %/% tile_size
      col_i <- (i - 1) %% tile_size
      cat(sprintf("      [%d,%d] ours=%d gdal=%d diff=%d\n",
                  col_i, row_i, our[i], gdal_px[i], our[i] - gdal_px[i]))
    }
  }
}

ds$close()
cat("\nDone.\n")
