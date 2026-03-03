## Test resampling algorithms against gdalwarp
##
## For each of near, bilinear, cubic, lanczos:
##   1. Warp a GEBCO tile with our Rust kernel
##   2. Warp the same tile with gdalwarp -r <method>
##   3. Compare pixel values
##
## Bilinear/cubic/lanczos won't be bit-identical (floating-point path
## differences, weight normalisation, edge handling). We check:
##   - Correlation > 0.9999
##   - Max absolute difference small relative to data range
##   - No systematic bias (mean diff ≈ 0)

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

## Pick a single 256×256 tile near the centre (away from edges/antimeridian)
tile_size <- 256L
pixel_res <- 2000
out_xmin <- -2000000
out_ymin <- -1500000

## Tile [4,3] — comfortably in the middle
tx <- 4; ty <- 3
dst_col_off <- (tx - 1L) * tile_size
dst_row_off <- (ty - 1L) * tile_size
tile_xmin <- out_xmin + dst_col_off * pixel_res
tile_xmax <- tile_xmin + tile_size * pixel_res
tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) - dst_row_off * pixel_res
tile_ymin <- tile_ymax - tile_size * pixel_res
tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
tile_gt   <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))

## Map of resampling methods and their source window padding
methods <- list(
  near     = list(gdal_flag = "near",     padding = 0L),
  bilinear = list(gdal_flag = "bilinear", padding = 1L),
  cubic    = list(gdal_flag = "cubic",    padding = 2L),
  lanczos  = list(gdal_flag = "lanczos",  padding = 3L)
)

cat("Testing resampling kernels on GEBCO Fiji LCC tile [4,3]\n")
cat(sprintf("Tile extent: [%.0f, %.0f, %.0f, %.0f]\n\n",
            tile_ext[1], tile_ext[2], tile_ext[3], tile_ext[4]))

for (method_name in names(methods)) {
  m <- methods[[method_name]]

  ## --- Compute source window with appropriate padding ---
  sw <- rust_compute_source_window(
    src_crs  = src_crs,
    src_gt   = src_gt,
    src_dim  = as.integer(src_dim),
    dst_crs  = out_crs,
    dst_gt   = tile_gt,
    dst_off  = c(0L, 0L),
    dst_size = c(tile_size, tile_size),
    resample_padding = m$padding
  )

  ## Snap to COG blocks
  src_block <- ds$getBlockSize(band = 1)
  src_tw <- src_block[1]; src_th <- src_block[2]
  read_xoff  <- (sw$xoff %/% src_tw) * src_tw
  read_yoff  <- (sw$yoff %/% src_th) * src_th
  read_xend  <- min(ceiling((sw$xoff + sw$xsize) / src_tw) * src_tw, src_dim[1])
  read_yend  <- min(ceiling((sw$yoff + sw$ysize) / src_th) * src_th, src_dim[2])
  read_xsize <- as.integer(read_xend - read_xoff)
  read_ysize <- as.integer(read_yend - read_yoff)

  ## --- Read source ---
  src_pixels <- ds$read(
    band = 1,
    xoff = read_xoff, yoff = read_yoff,
    xsize = read_xsize, ysize = read_ysize,
    out_xsize = read_xsize, out_ysize = read_ysize
  )

  ## --- Rust warp ---
  our_pixels <- rust_warp_resample(
    src_crs     = src_crs,
    src_gt      = src_gt,
    dst_crs     = out_crs,
    dst_gt      = tile_gt,
    dst_dim     = c(tile_size, tile_size),
    src_pixels  = src_pixels,
    src_ncol    = read_xsize,
    src_nrow    = read_ysize,
    src_col_off = as.integer(read_xoff),
    src_row_off = as.integer(read_yoff),
    nodata      = as.integer(nodata),
    max_error   = 0.125,
    resample    = method_name
  )

  ## --- GDAL warp ---
  gdal_file <- tempfile(fileext = ".tif")
  warp(ds, gdal_file, out_crs,
       cl_arg = c("-ts", tile_size, tile_size,
                  "-te", tile_ext[c(1, 3, 2, 4)],
                  "-r", m$gdal_flag, "-ovr", "NONE"))
  gds <- new(GDALRaster, gdal_file)
  gdal_pixels <- gds$read(band = 1, xoff = 0, yoff = 0,
                           xsize = tile_size, ysize = tile_size,
                           out_xsize = tile_size, out_ysize = tile_size)
  gds$close()

  ## --- Compare ---
  n_total <- tile_size * tile_size
  n_match <- sum(our_pixels == gdal_pixels)

  ## Only compare where both have valid data
  both_valid <- our_pixels != nodata & gdal_pixels != nodata
  n_both     <- sum(both_valid)
  diffs      <- our_pixels[both_valid] - gdal_pixels[both_valid]

  if (n_both > 0) {
    corr <- cor(our_pixels[both_valid], gdal_pixels[both_valid])
    cat(sprintf("  %-8s  exact match: %5d/%d (%.1f%%)  |  valid: %d  corr: %.6f  mean_diff: %.2f  max_diff: %d\n",
                method_name, n_match, n_total, 100 * n_match / n_total,
                n_both, corr, mean(diffs), max(abs(diffs))))
  } else {
    cat(sprintf("  %-8s  no valid overlapping pixels!\n", method_name))
  }
}

ds$close()
cat("\nDone.\n")
