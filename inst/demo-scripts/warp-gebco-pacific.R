## Warp GEBCO 2024 bathymetry to a Fiji-centred Lambert Conformal Conic
##
## Architecture: per-dest-tile source tile mapping
##
## For each destination tile:
##   1. Transform dest pixel coords → source pixel coords (via gdalraster)
##   2. Determine which COG internal tiles overlap that source footprint
##   3. Read only those source tiles, assemble into local buffer
##   4. Warp from that local buffer using ApproxTransformer
##
## This is the core indexing step for a caching pipeline:
##   dest_tile → {source_tile_ids}
##
## Handles antimeridian correctly because each dest tile resolves its own
## source footprint independently — no monolithic source window.

library(gdalraster)
library(vaster)

## --- Configuration ---

cog_url <- "/vsicurl/https://data.source.coop/alexgleith/gebco-2024/GEBCO_2024.tif"

## Fiji-centred Lambert Conformal Conic
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

## Output extent in LCC metres — ~4000 x 3000 km region
## Covers Vanuatu, New Caledonia, Tonga, Samoa, Cook Islands
out_xmin <- -2000000
out_xmax <-  2000000
out_ymin <- -1500000
out_ymax <-  1500000

tile_size  <- 256L
pixel_res  <- 2000     # metres — ~2km regional overview

## --- Open source dataset, query its internal tiling ---

ds <- new(GDALRaster, cog_url)
src_crs    <- ds$getProjectionRef()
src_gt     <- ds$getGeoTransform()
src_dim    <- c(ds$getRasterXSize(), ds$getRasterYSize())
nodata     <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- -32768L

## COG internal tile structure
src_block  <- ds$getBlockSize(band = 1)  # c(block_width, block_height)
src_tw     <- src_block[1]
src_th     <- src_block[2]
src_ntx    <- ceiling(src_dim[1] / src_tw)
src_nty    <- ceiling(src_dim[2] / src_th)

cat(sprintf("Source: %d x %d, block %d x %d (%d x %d blocks)\n",
            src_dim[1], src_dim[2], src_tw, src_th, src_ntx, src_nty))
cat(sprintf("Source gt: %s\n", paste(src_gt, collapse = ", ")))

## --- Compute output grid ---

out_ncol <- as.integer(ceiling((out_xmax - out_xmin) / pixel_res))
out_nrow <- as.integer(ceiling((out_ymax - out_ymin) / pixel_res))

## Snap extent
out_xmax_snapped <- out_xmin + out_ncol * pixel_res
out_ymax_snapped <- out_ymin + out_nrow * pixel_res
out_extent <- c(out_xmin, out_xmax_snapped, out_ymin, out_ymax_snapped)
out_crs <- srs_to_wkt(fiji_lcc)

cat(sprintf("Output grid: %d x %d = %d pixels\n",
            out_ncol, out_nrow, out_ncol * out_nrow))

## --- Tile grid ---

n_tiles_x <- ceiling(out_ncol / tile_size)
n_tiles_y <- ceiling(out_nrow / tile_size)
cat(sprintf("Dest tile grid: %d x %d = %d tiles\n\n",
            n_tiles_x, n_tiles_y, n_tiles_x * n_tiles_y))

## --- Helper: source footprint for a destination tile ---
##
## Given a dest tile extent + dimension, transform a sampling of dest
## pixels back to source pixel space and return the range of source
## pixels needed plus which COG internal tiles that spans.
##
## Returns NULL if no valid source pixels (tile is entirely outside
## the source domain).

source_footprint <- function(tile_extent, tile_dim, out_crs, src_crs,
                             src_gt, src_dim, src_tw, src_th) {
  ## Sample dest tile edges + a cross through the middle
  ## (edges alone can miss interior features for conic projections)
  n_edge <- max(tile_dim)
  ex <- tile_extent

  ## Edge points
  top    <- cbind(seq(ex[1], ex[2], length.out = n_edge),
                  rep(ex[4], n_edge))
  bottom <- cbind(seq(ex[1], ex[2], length.out = n_edge),
                  rep(ex[3], n_edge))
  left   <- cbind(rep(ex[1], n_edge),
                  seq(ex[3], ex[4], length.out = n_edge))
  right  <- cbind(rep(ex[2], n_edge),
                  seq(ex[3], ex[4], length.out = n_edge))
  ## Cross through middle
  hmid   <- cbind(seq(ex[1], ex[2], length.out = n_edge),
                  rep((ex[3] + ex[4]) / 2, n_edge))
  vmid   <- cbind(rep((ex[1] + ex[2]) / 2, n_edge),
                  seq(ex[3], ex[4], length.out = n_edge))

  pts <- rbind(top, bottom, left, right, hmid, vmid)

  ## Transform to source CRS (geographic coords for GEBCO)
  src_pts <- transform_xy(pts, out_crs, src_crs)
  valid <- !is.na(src_pts[, 1]) & !is.na(src_pts[, 2])
  if (sum(valid) == 0) return(NULL)

  ## Convert to source pixel coords
  src_col <- (src_pts[valid, 1] - src_gt[1]) / src_gt[2]
  src_row <- (src_pts[valid, 2] - src_gt[4]) / src_gt[6]

  ## Source pixel bounds (1px padding, clamp to raster)
  col_min <- max(0L, as.integer(floor(min(src_col))) - 1L)
  col_max <- min(src_dim[1] - 1L, as.integer(ceiling(max(src_col))) + 1L)
  row_min <- max(0L, as.integer(floor(min(src_row))) - 1L)
  row_max <- min(src_dim[2] - 1L, as.integer(ceiling(max(src_row))) + 1L)

  ## Which COG internal tiles does this span?
  tile_col_min <- col_min %/% src_tw
  tile_col_max <- col_max %/% src_tw
  tile_row_min <- row_min %/% src_th
  tile_row_max <- row_max %/% src_th

  ## Read window aligned to COG tile boundaries (for efficient I/O)
  read_xoff  <- tile_col_min * src_tw
  read_yoff  <- tile_row_min * src_th
  read_xsize <- min((tile_col_max + 1L) * src_tw, src_dim[1]) - read_xoff
  read_ysize <- min((tile_row_max + 1L) * src_th, src_dim[2]) - read_yoff

  list(
    col_range = c(col_min, col_max),
    row_range = c(row_min, row_max),
    src_tile_cols = tile_col_min:tile_col_max,
    src_tile_rows = tile_row_min:tile_row_max,
    n_src_tiles = (tile_col_max - tile_col_min + 1L) *
      (tile_row_max - tile_row_min + 1L),
    read_xoff  = read_xoff,
    read_yoff  = read_yoff,
    read_xsize = read_xsize,
    read_ysize = read_ysize
  )
}

## --- Output directory ---

out_dir <- "tiles/tiles_fiji_lcc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## --- Phase 1: Build the dest→source tile mapping ---

cat("Phase 1: mapping dest tiles to source tiles...\n")

tile_map <- list()
total_src_tiles <- 0

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

    fp <- source_footprint(tile_extent, c(dst_ncol, dst_nrow),
                           out_crs, src_crs, src_gt, src_dim,
                           src_tw, src_th)

    key <- sprintf("%03d_%03d", tx, ty)
    tile_map[[key]] <- list(
      tx = tx, ty = ty,
      dst_col_off = dst_col_off, dst_row_off = dst_row_off,
      dst_ncol = dst_ncol, dst_nrow = dst_nrow,
      tile_extent = tile_extent,
      footprint = fp
    )
    if (!is.null(fp)) total_src_tiles <- total_src_tiles + fp$n_src_tiles
  }
}

active_tiles <- Filter(function(t) !is.null(t$footprint), tile_map)
cat(sprintf("  %d dest tiles have source data (of %d total)\n",
            length(active_tiles), length(tile_map)))
cat(sprintf("  Total source tile reads: %d (before dedup/caching)\n",
            total_src_tiles))

## Summary of source tile usage
all_src_tile_ids <- unlist(lapply(active_tiles, function(t) {
  fp <- t$footprint
  as.vector(outer(fp$src_tile_cols, fp$src_tile_rows,
                  function(c, r) paste(c, r, sep = ",")))
}))
unique_src_tiles <- unique(all_src_tile_ids)
cat(sprintf("  Unique source tiles needed: %d\n", length(unique_src_tiles)))

## Source tile read duplication factor
cat(sprintf("  Duplication factor: %.1fx (caching would reduce %d reads to %d)\n\n",
            length(all_src_tile_ids) / length(unique_src_tiles),
            length(all_src_tile_ids), length(unique_src_tiles)))

## --- Phase 2: Warp each dest tile from its source footprint ---

cat("Phase 2: warping tiles...\n")

tile_count <- 0
skip_count <- 0
t_total <- proc.time()
read_time_total <- 0
warp_time_total <- 0

for (tm in active_tiles) {
  fp <- tm$footprint

  ## Read source pixels for this tile's footprint
  t0 <- proc.time()
  src_pixels <- ds$read(
    band = 1,
    xoff = fp$read_xoff, yoff = fp$read_yoff,
    xsize = fp$read_xsize, ysize = fp$read_ysize,
    out_xsize = fp$read_xsize, out_ysize = fp$read_ysize
  )
  t_read <- (proc.time() - t0)[3]
  read_time_total <- read_time_total + t_read

  ## Tile geotransform
  tile_gt <- extent_dim_to_gt(tm$tile_extent, c(tm$dst_ncol, tm$dst_nrow))

  ## Warp using ApproxTransformer
  t0 <- proc.time()
  tile_pixels <- rust_warp_approx(
    src_crs     = src_crs,
    src_gt      = src_gt,
    dst_crs     = out_crs,
    dst_gt      = tile_gt,
    dst_dim     = as.integer(c(tm$dst_ncol, tm$dst_nrow)),
    src_pixels  = src_pixels,
    src_ncol    = as.integer(fp$read_xsize),
    src_nrow    = as.integer(fp$read_ysize),
    src_col_off = as.integer(fp$read_xoff),
    src_row_off = as.integer(fp$read_yoff),
    nodata      = as.integer(nodata),
    max_error   = 0.125
  )
  t_warp <- (proc.time() - t0)[3]
  warp_time_total <- warp_time_total + t_warp

  n_valid <- sum(tile_pixels != nodata)
  if (n_valid == 0) {
    skip_count <- skip_count + 1
    next
  }

  ## Write tile
  tile_file <- file.path(out_dir, sprintf("tile_%03d_%03d.tif", tm$tx, tm$ty))
  driver <- gdalraster::create(
    dst_filename = tile_file,
    format = "GTiff",
    xsize = tm$dst_ncol, ysize = tm$dst_nrow, nbands = 1,
    dataType = "Int16",
    options = c("COMPRESS=DEFLATE", "TILED=YES")
  )
  driver$setProjectionRef(out_crs)
  driver$setGeoTransform(tile_gt)
  driver$setNoDataValue(band = 1, as.double(nodata))
  driver$write(band = 1, xoff = 0, yoff = 0,
               xsize = tm$dst_ncol, ysize = tm$dst_nrow,
               as.integer(tile_pixels))
  driver$close()

  tile_count <- tile_count + 1
  if (tile_count %% 10 == 0 || tile_count <= 5) {
    cat(sprintf("  tile %d [%d,%d] src %.1fMpx (%d blks, %.3fs read) warp %.3fs, %d valid\n",
                tile_count, tm$tx, tm$ty,
                fp$read_xsize * fp$read_ysize / 1e6,
                fp$n_src_tiles, t_read,
                t_warp, n_valid))
  }
}

t_total_elapsed <- (proc.time() - t_total)[3]
cat(sprintf("\n%d tiles written, %d skipped (all nodata) in %.1fs\n",
            tile_count, skip_count, t_total_elapsed))
cat(sprintf("  read: %.1fs  warp: %.1fs  write+overhead: %.1fs\n",
            read_time_total, warp_time_total,
            t_total_elapsed - read_time_total - warp_time_total))
cat(sprintf("  avg per active tile: read %.3fs, warp %.3fs\n",
            read_time_total / max(length(active_tiles), 1),
            warp_time_total / max(length(active_tiles), 1)))

## --- Validation: compare centre tile vs GDAL ---

cat("\n=== Validation: centre tile vs gdalwarp (-ovr NONE) ===\n")

ctx <- ceiling(n_tiles_x / 2)
cty <- ceiling(n_tiles_y / 2)
centre_file <- file.path(out_dir, sprintf("tile_%03d_%03d.tif", ctx, cty))

if (file.exists(centre_file)) {
  cds <- new(GDALRaster, centre_file)
  tile_gt_v <- cds$getGeoTransform()
  centre_ext <- c(
    tile_gt_v[1],
    tile_gt_v[1] + cds$getRasterXSize() * tile_gt_v[2],
    tile_gt_v[4] + cds$getRasterYSize() * tile_gt_v[6],
    tile_gt_v[4]
  )
  centre_dim <- c(cds$getRasterXSize(), cds$getRasterYSize())
  our_pixels <- cds$read(band = 1, xoff = 0, yoff = 0,
                         xsize = centre_dim[1], ysize = centre_dim[2],
                         out_xsize = centre_dim[1], out_ysize = centre_dim[2])
  cds$close()

  ## GDAL warp for comparison (full-res source, no overview)
  gdal_file <- tempfile(fileext = ".tif")
  warp(ds, gdal_file, out_crs,
       cl_arg = c("-ts", centre_dim,
                  "-te", centre_ext[c(1, 3, 2, 4)],
                  "-r", "near", "-ovr", "NONE"))
  gds <- new(GDALRaster, gdal_file)
  gdal_pixels <- gds$read(band = 1, xoff = 0, yoff = 0,
                          xsize = centre_dim[1], ysize = centre_dim[2],
                          out_xsize = centre_dim[1], out_ysize = centre_dim[2])
  gds$close()

  n_total <- prod(centre_dim)
  n_match <- sum(our_pixels == gdal_pixels)
  cat(sprintf("  Approx vs GDAL: %d / %d match (%.1f%%)\n",
              n_match, n_total, 100 * n_match / n_total))

  diffs <- our_pixels - gdal_pixels
  cat(sprintf("  Diff range: [%d, %d], mean abs diff: %.2f\n",
              min(diffs), max(diffs), mean(abs(diffs))))
}

## --- Build VRT mosaic ---

tile_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(tile_files) > 1) {
  vrt_file <- file.path(out_dir, "fiji_lcc_mosaic.vrt")
  buildVRT(vrt_file, tile_files)
  cat(sprintf("\nVRT mosaic: %s (%d tiles)\n", vrt_file, length(tile_files)))
}

ds$close()
cat("Done.\n")
