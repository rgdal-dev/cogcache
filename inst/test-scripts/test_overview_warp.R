## Test overview-aware warp pipeline
##
## Demonstrates: select_overview() -> compute_source_window() at overview
## level -> read from overview -> warp at ~1:1 -> compare to GDAL
##
## The key insight: GDAL internally scales the interpolation kernel when
## downsampling. By reading from an appropriate overview first, we keep
## the kernel at ~1:1 where our bilinear is bit-identical to GDAL's.

library(gdalraster)
library(vaster)
source("R/overview.R")

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
nodata  <- -32768L

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

## Step 1: compute source window at full resolution (for overview selection)
sw_full <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
cat(sprintf("Full-res source window: %d x %d (ratio %.1f x %.1f)\n",
            sw_full$xsize, sw_full$ysize,
            sw_full$xsize / tile_size, sw_full$ysize / tile_size))

## Step 2: select overview
ovr <- select_overview(ds, src_gt = src_gt, src_dim = src_dim,
                       sw = sw_full, dst_dim = c(tile_size, tile_size))
cat(sprintf("\nSelected overview level %d of %d\n", ovr$ovr_level, ovr$n_overviews))
cat(sprintf("  Overview dims: %d x %d\n", ovr$ovr_dim[1], ovr$ovr_dim[2]))
cat(sprintf("  Overview pixel size: %.8f x %.8f\n", ovr$ovr_gt[2], abs(ovr$ovr_gt[6])))
cat(sprintf("  Expected ratio: %.2f x %.2f\n", ovr$ratio_x, ovr$ratio_y))

## Step 3: compute source window at overview resolution
sw_ovr <- rust_compute_source_window(
  src_crs = src_crs, src_gt = ovr$ovr_gt, src_dim = ovr$ovr_dim,
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
cat(sprintf("\nOverview source window: %d x %d (ratio %.2f x %.2f)\n",
            sw_ovr$xsize, sw_ovr$ysize,
            sw_ovr$xsize / tile_size, sw_ovr$ysize / tile_size))

## Step 4: read from overview level
## gdalraster read with overview: we need to read from the overview band.
## For a GeoTIFF with internal overviews, we open the overview subdataset
## or use RasterIO with the overview band.
##
## gdalraster::read() reads from the base band and resamples if
## out_xsize != xsize. To read from a specific overview we can:
## (a) Use GDAL_OVERVIEW_LEVEL open option
## (b) Use ds$read() with appropriate xsize/out_xsize to trigger
##     the overview selection
## (c) Open the overview as a separate dataset via /vsicurl/...?ovr=N
##
## Method (a) is cleanest: reopen with OVERVIEW_LEVEL config

if (ovr$ovr_level > 0) {
  ## Open at the selected overview level
  ## GDAL open option OVERVIEW_LEVEL is 0-based (0 = first overview)
  ovr_open_level <- ovr$ovr_level - 1L
  ovr_url <- cog_url
  ds_ovr <- new(GDALRaster, ovr_url,
                TRUE,
                open_options = paste0("OVERVIEW_LEVEL=", ovr_open_level))

  ovr_xsize_actual <- ds_ovr$getRasterXSize()
  ovr_ysize_actual <- ds_ovr$getRasterYSize()
  cat(sprintf("  Actual overview dims: %d x %d\n",
              ovr_xsize_actual, ovr_ysize_actual))

  ## Adjust the overview geotransform to match actual dimensions
  ## (in case our estimate differs slightly from the real overview)
  ovr_gt_actual <- ds_ovr$getGeoTransform()

  ## Recompute source window with actual overview geotransform
  sw_ovr <- rust_compute_source_window(
    src_crs = src_crs, src_gt = ovr_gt_actual,
    src_dim = as.integer(c(ovr_xsize_actual, ovr_ysize_actual)),
    dst_crs = out_crs, dst_gt = tile_gt,
    dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
    resample_padding = 1L
  )
  cat(sprintf("  Recomputed sw: %d x %d (ratio %.2f x %.2f)\n",
              sw_ovr$xsize, sw_ovr$ysize,
              sw_ovr$xsize / tile_size, sw_ovr$ysize / tile_size))

  ## Tile-aligned read from overview
  ovr_block <- ds_ovr$getBlockSize(band = 1)
  ovr_tw <- ovr_block[1]; ovr_th <- ovr_block[2]
  read_xoff  <- (sw_ovr$xoff %/% ovr_tw) * ovr_tw
  read_yoff  <- (sw_ovr$yoff %/% ovr_th) * ovr_th
  read_xend  <- min(ceiling((sw_ovr$xoff + sw_ovr$xsize) / ovr_tw) * ovr_tw,
                    ovr_xsize_actual)
  read_yend  <- min(ceiling((sw_ovr$yoff + sw_ovr$ysize) / ovr_th) * ovr_th,
                    ovr_ysize_actual)
  read_xsize <- as.integer(read_xend - read_xoff)
  read_ysize <- as.integer(read_yend - read_yoff)

  src_pixels <- ds_ovr$read(band = 1, xoff = read_xoff, yoff = read_yoff,
                            xsize = read_xsize, ysize = read_ysize,
                            out_xsize = read_xsize, out_ysize = read_ysize)
  use_gt <- ovr_gt_actual
  ds_ovr$close()
} else {
  ## Full resolution — existing pipeline
  src_block <- ds$getBlockSize(band = 1)
  src_tw <- src_block[1]; src_th <- src_block[2]
  read_xoff  <- (sw_full$xoff %/% src_tw) * src_tw
  read_yoff  <- (sw_full$yoff %/% src_th) * src_th
  read_xend  <- min(ceiling((sw_full$xoff + sw_full$xsize) / src_tw) * src_tw,
                    src_dim[1])
  read_yend  <- min(ceiling((sw_full$yoff + sw_full$ysize) / src_th) * src_th,
                    src_dim[2])
  read_xsize <- as.integer(read_xend - read_xoff)
  read_ysize <- as.integer(read_yend - read_yoff)

  src_pixels <- ds$read(band = 1, xoff = read_xoff, yoff = read_yoff,
                        xsize = read_xsize, ysize = read_ysize,
                        out_xsize = read_xsize, out_ysize = read_ysize)
  use_gt <- src_gt
}

ds$close()

cat(sprintf("\nRead buffer: %d x %d at (%d, %d)\n",
            read_xsize, read_ysize, read_xoff, read_yoff))

## Step 5: warp at ~1:1 ratio
for (resample in c("near", "bilinear", "cubic", "lanczos")) {
  rp <- switch(resample, near = 0L, bilinear = 1L, cubic = 2L, lanczos = 3L)

  our <- rust_warp_resample(
    src_crs = src_crs, src_gt = use_gt, dst_crs = out_crs, dst_gt = tile_gt,
    dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
    src_ncol = read_xsize, src_nrow = read_ysize,
    src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
    nodata = nodata, max_error = 0.125, resample = resample
  )

  ## GDAL reference (with -ovr AUTO to let GDAL pick the same overview)
  gdal_file <- tempfile(fileext = ".tif")
  warp(cog_url, gdal_file, out_crs,
       cl_arg = c("-ts", tile_size, tile_size,
                  "-te", tile_ext[c(1, 3, 2, 4)],
                  "-r", resample, "-ovr", "AUTO"))
  gds <- new(GDALRaster, gdal_file)
  gdal_px <- gds$read(band = 1, xoff = 0, yoff = 0,
                      xsize = tile_size, ysize = tile_size,
                      out_xsize = tile_size, out_ysize = tile_size)
  gds$close()

  both <- our != nodata & gdal_px != nodata
  diffs <- our[both] - gdal_px[both]
  n_match <- sum(our == gdal_px)

  cat(sprintf("  %-10s match %5d/%d (%5.1f%%)  corr %.6f  mean_diff %6.2f  max_diff %d\n",
              resample, n_match, tile_size^2, 100 * n_match / tile_size^2,
              cor(our[both], gdal_px[both]),
              mean(diffs), max(abs(diffs))))
}

cat("\nDone.\n")
