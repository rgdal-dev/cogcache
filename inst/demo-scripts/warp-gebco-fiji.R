## Warp GEBCO 2024 bathymetry to a Fiji-centred Lambert Conformal Conic
##
## Uses rust_compute_source_window for source footprint computation
## and rust_warp_approx for the warp kernel.
##
## Two phases:
##   Phase 1: map each dest tile to its source window (Rust)
##   Phase 2: read source pixels + warp + write output tiles

library(gdalraster)
library(vaster)

## --- Configuration ---

cog_url <- "/vsicurl/https://data.source.coop/alexgleith/gebco-2024/GEBCO_2024.tif"

fiji_lcc <- paste(
  "+proj=lcc",
  "+lat_0=-18",
  "+lon_0=178",
  "+lat_1=-10",
  "+lat_2=-25",
  "+datum=WGS84",
  "+units=m",
  "+no_defs"
)

out_xmin <- -2000000;  out_xmax <- 2000000
out_ymin <- -1500000;  out_ymax <- 1500000
tile_size  <- 256L
pixel_res  <- 2000

## --- Open source dataset ---

ds <- new(GDALRaster, cog_url)
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
nodata  <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- -32768L

src_block <- ds$getBlockSize(band = 1)
src_tw <- src_block[1]
src_th <- src_block[2]

cat(sprintf("Source: %d x %d, block %d x %d\n",
            src_dim[1], src_dim[2], src_tw, src_th))

## --- Compute output grid ---

out_ncol <- as.integer(ceiling((out_xmax - out_xmin) / pixel_res))
out_nrow <- as.integer(ceiling((out_ymax - out_ymin) / pixel_res))
out_extent <- c(out_xmin, out_xmin + out_ncol * pixel_res,
                out_ymin, out_ymin + out_nrow * pixel_res)
out_crs <- srs_to_wkt(fiji_lcc)

n_tiles_x <- ceiling(out_ncol / tile_size)
n_tiles_y <- ceiling(out_nrow / tile_size)
cat(sprintf("Output: %d x %d, %d x %d = %d tiles\n\n",
            out_ncol, out_nrow, n_tiles_x, n_tiles_y,
            n_tiles_x * n_tiles_y))

## --- Phase 1: Map dest tiles to source windows (Rust) ---

cat("Phase 1: computing source windows...\n")

tiles <- list()
total_src_blks <- 0

for (ty in seq_len(n_tiles_y)) {
  for (tx in seq_len(n_tiles_x)) {
    dst_col_off <- (tx - 1L) * tile_size
    dst_row_off <- (ty - 1L) * tile_size
    dst_ncol <- min(tile_size, out_ncol - dst_col_off)
    dst_nrow <- min(tile_size, out_nrow - dst_row_off)

    ## Tile geographic extent + geotransform
    tile_xmin <- out_extent[1] + dst_col_off * pixel_res
    tile_xmax <- tile_xmin + dst_ncol * pixel_res
    tile_ymax <- out_extent[4] - dst_row_off * pixel_res
    tile_ymin <- tile_ymax - dst_nrow * pixel_res
    tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
    tile_gt   <- extent_dim_to_gt(tile_ext, c(dst_ncol, dst_nrow))

    ## Rust: compute source window
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

    if (is.null(sw) || sw$xsize == 0 || sw$ysize == 0) next

    ## Snap read window to COG block boundaries for efficient I/O
    read_xoff  <- (sw$xoff %/% src_tw) * src_tw
    read_yoff  <- (sw$yoff %/% src_th) * src_th
    read_xend  <- min(ceiling((sw$xoff + sw$xsize) / src_tw) * src_tw, src_dim[1])
    read_yend  <- min(ceiling((sw$yoff + sw$ysize) / src_th) * src_th, src_dim[2])
    read_xsize <- as.integer(read_xend - read_xoff)
    read_ysize <- as.integer(read_yend - read_yoff)

    n_blks <- ceiling(read_xsize / src_tw) * ceiling(read_ysize / src_th)
    total_src_blks <- total_src_blks + n_blks

    tiles[[length(tiles) + 1L]] <- list(
      tx = tx, ty = ty,
      dst_ncol = dst_ncol, dst_nrow = dst_nrow,
      tile_ext = tile_ext, tile_gt = tile_gt,
      read_xoff = read_xoff, read_yoff = read_yoff,
      read_xsize = read_xsize, read_ysize = read_ysize,
      n_blks = n_blks, fill_ratio = sw$fill_ratio,
      n_failed = sw$n_failed
    )
  }
}

cat(sprintf("  %d active tiles, %d total source block reads\n\n",
            length(tiles), total_src_blks))

## --- Output directory ---

out_dir <- "tiles/tiles_fiji_lcc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## --- Phase 2: Read + Warp + Write ---

cat("Phase 2: warping tiles...\n")

tile_count <- 0
skip_count <- 0
t_total <- proc.time()
read_total <- 0
warp_total <- 0

for (tm in tiles) {
  ## Read source pixels
  t0 <- proc.time()
  src_pixels <- ds$read(
    band = 1,
    xoff = tm$read_xoff, yoff = tm$read_yoff,
    xsize = tm$read_xsize, ysize = tm$read_ysize,
    out_xsize = tm$read_xsize, out_ysize = tm$read_ysize
  )
  t_read <- (proc.time() - t0)[3]
  read_total <- read_total + t_read

  ## Warp using ApproxTransformer
  t0 <- proc.time()
  tile_pixels <- rust_warp_approx(
    src_crs     = src_crs,
    src_gt      = src_gt,
    dst_crs     = out_crs,
    dst_gt      = tm$tile_gt,
    dst_dim     = as.integer(c(tm$dst_ncol, tm$dst_nrow)),
    src_pixels  = src_pixels,
    src_ncol    = as.integer(tm$read_xsize),
    src_nrow    = as.integer(tm$read_ysize),
    src_col_off = as.integer(tm$read_xoff),
    src_row_off = as.integer(tm$read_yoff),
    nodata      = as.integer(nodata),
    max_error   = 0.125
  )
  ximage(matrix(tile_pixels, as.integer(tm$dst_nrow), byrow = TRUE), tm$tile_ext, add = T)
  t_warp <- (proc.time() - t0)[3]
  warp_total <- warp_total + t_warp

  n_valid <- sum(tile_pixels != nodata)
  if (n_valid == 0) {
    skip_count <- skip_count + 1
    next
  }

  ## Write tile
  tile_file <- file.path(out_dir, sprintf("tile_%03d_%03d.tif", tm$tx, tm$ty))
  drv <- gdalraster::create(
    dst_filename = tile_file,
    format = "GTiff",
    xsize = tm$dst_ncol, ysize = tm$dst_nrow, nbands = 1,
    dataType = "Int16",
    options = c("COMPRESS=DEFLATE", "TILED=YES"),
    return_obj = TRUE
  )
  drv$setProjection(out_crs)
  drv$setGeoTransform(tm$tile_gt)
  drv$setNoDataValue(band = 1, as.double(nodata))
  drv$write(band = 1, xoff = 0, yoff = 0,
            xsize = tm$dst_ncol, ysize = tm$dst_nrow,
            as.integer(tile_pixels))
  drv$close()

  tile_count <- tile_count + 1
  if (tile_count %% 10 == 0 || tile_count <= 5) {
    cat(sprintf("  tile %d [%d,%d] src %.1fMpx (%d blks, %.3fs) warp %.3fs, %d valid\n",
                tile_count, tm$tx, tm$ty,
                tm$read_xsize * tm$read_ysize / 1e6,
                tm$n_blks, t_read, t_warp, n_valid))
  }
}

t_elapsed <- (proc.time() - t_total)[3]
cat(sprintf("\n%d tiles written, %d skipped in %.1fs\n",
            tile_count, skip_count, t_elapsed))
cat(sprintf("  read: %.1fs  warp: %.1fs  overhead: %.1fs\n",
            read_total, warp_total,
            t_elapsed - read_total - warp_total))

## --- Validation ---

cat("\n=== Validation: centre tile vs gdalwarp ===\n")

ctx <- ceiling(n_tiles_x / 2)
cty <- ceiling(n_tiles_y / 2)
centre_file <- file.path(out_dir, sprintf("tile_%03d_%03d.tif", ctx, cty))

if (file.exists(centre_file)) {
  cds <- new(GDALRaster, centre_file)
  cgt <- cds$getGeoTransform()
  cext <- c(cgt[1], cgt[1] + cds$getRasterXSize() * cgt[2],
            cgt[4] + cds$getRasterYSize() * cgt[6], cgt[4])
  cdim <- c(cds$getRasterXSize(), cds$getRasterYSize())
  our <- cds$read(band = 1, xoff = 0, yoff = 0,
                  xsize = cdim[1], ysize = cdim[2],
                  out_xsize = cdim[1], out_ysize = cdim[2])
  cds$close()

  gdal_file <- tempfile(fileext = ".tif")
  warp(ds, gdal_file, out_crs,
       cl_arg = c("-ts", cdim, "-te", cext[c(1, 3, 2, 4)],
                  "-r", "near", "-ovr", "NONE"))
  gds <- new(GDALRaster, gdal_file)
  gdal_px <- gds$read(band = 1, xoff = 0, yoff = 0,
                      xsize = cdim[1], ysize = cdim[2],
                      out_xsize = cdim[1], out_ysize = cdim[2])
  gds$close()

  cat(sprintf("  %d / %d match (%.1f%%)\n",
              sum(our == gdal_px), prod(cdim),
              100 * sum(our == gdal_px) / prod(cdim)))
}

## --- VRT mosaic ---

tile_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(tile_files) > 1) {
  vrt_file <- file.path(out_dir, "fiji_lcc_mosaic.vrt")
  buildVRT(vrt_file, tile_files)
  cat(sprintf("\nVRT mosaic: %s (%d tiles)\n", vrt_file, length(tile_files)))
}

ds$close()
cat("Done.\n")
