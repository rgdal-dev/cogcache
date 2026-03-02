## IBCSO → Mawson web tiles
##
## Source: IBCSO v2 ice-surface COG (DEFLATE), EPSG:9354
## Output: 128×128 EPSG:3857 tiles at ~500m resolution (zoom 8, 128px tiles)
##
## Mawson extent in EPSG:3857:
##   [6553000, 7459000, -10840000, -9841000]
## (~906km × 999km, covering the Mawson coast region)

library(gdalraster)
library(reproj)
library(vaster)

## ---- Config ----
cog_url <- "/vsicurl/https://github.com/mdsumner/ibcso-cog/releases/download/latest/IBCSO_v2_ice-surface-cog_DEFLATE.tif"
tile_size <- 128L
out_dir <- "tiles/tiles_mawson"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Output extent around Mawson in EPSG:3857
dest_extent <- c(6553000, 7459000, -10840000, -9841000)

## Snap to web mercator zoom 8 grid (128px tiles)
## At zoom 8 with 128px tiles:
##   world extent in 3857 = c(-20037508.34, 20037508.34, -20037508.34, 20037508.34)
##   n_tiles = 2^8 = 256 tiles per axis
##   tile_width = 40075016.68 / 256 = 156543.03 m
##   pixel_res = 156543.03 / 128 = 1223.0 m  -- too coarse
##
## Zoom 9 with 128px:
##   tile_width = 40075016.68 / 512 = 78271.5 m
##   pixel_res = 78271.5 / 128 = 611.5 m  -- close to 500m
##
## Zoom 10 with 128px:
##   tile_width = 40075016.68 / 1024 = 39135.8 m
##   pixel_res = 39135.8 / 128 = 305.7 m  -- oversamples 500m source

zoom <- 10L  ## 611m pixel resolution, close to 500m source
world_half <- 20037508.342789244
tile_width <- 2 * world_half / 2^zoom

## Find tile indices covering our extent
tx_min <- as.integer(floor((dest_extent[1] + world_half) / tile_width))
tx_max <- as.integer(floor((dest_extent[2] + world_half) / tile_width))
ty_min <- as.integer(floor((world_half - dest_extent[4]) / tile_width))
ty_max <- as.integer(floor((world_half - dest_extent[3]) / tile_width))

tiles <- expand.grid(x = tx_min:tx_max, y = ty_min:ty_max)
cat(sprintf("Tile grid: x=%d..%d, y=%d..%d (%d tiles)\n",
            tx_min, tx_max, ty_min, ty_max, nrow(tiles)))
cat(sprintf("Zoom: %d, tile size: %dx%d, pixel resolution: %.1f m\n",
            zoom, tile_size, tile_size, tile_width / tile_size))

## ---- Source grid ----
cat("\n=== Opening source COG ===\n")

ds <- new(GDALRaster, cog_url)
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
cat(sprintf("  Source: %d x %d pixels, CRS: %s\n",
            src_dim[1], src_dim[2],
            substr(src_crs, 1, 60)))

nodata <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- -32768L
cat(sprintf("  NoData: %d\n", nodata))

dest_crs <- srs_to_wkt("EPSG:3857")

## ---- Tile extent helper ----
tile_extent_3857 <- function(z, tx, ty) {
  tw <- 2 * world_half / 2^z
  xmin <- -world_half + tx * tw
  xmax <- xmin + tw
  ymax <- world_half - ty * tw
  ymin <- ymax - tw
  c(xmin, xmax, ymin, ymax)
}

## ---- Warp tiles ----
cat(sprintf("\n=== Warping %d tiles ===\n\n", nrow(tiles)))

t_total <- proc.time()

for (i in seq_len(nrow(tiles))) {
  tx <- tiles$x[i]
  ty <- tiles$y[i]

  ## Destination grid for this tile
  dest_ext <- tile_extent_3857(zoom, tx, ty)
  dest_gt  <- extent_dim_to_gt(dest_ext, c(tile_size, tile_size))

  ## Compute source pixel coordinates for every destination pixel
  coords <- rust_gen_img_proj_transform(
    src_crs = src_crs, src_gt = src_gt,
    dst_crs = dest_crs, dst_gt = dest_gt,
    dst_dim = c(tile_size, tile_size)
  )

  valid <- !is.na(coords$src_x)
  if (sum(valid) == 0) {
    cat(sprintf("  [%d/%d] z=%d x=%d y=%d — no valid source pixels, skipping\n",
                i, nrow(tiles), zoom, tx, ty))
    next
  }

  src_col_range <- range(as.integer(coords$src_x[valid]))
  src_row_range <- range(as.integer(coords$src_y[valid]))

  ## Clamp to source raster bounds
  src_xoff  <- max(0L, src_col_range[1])
  src_yoff  <- max(0L, src_row_range[1])
  src_xend  <- min(src_dim[1] - 1L, src_col_range[2])
  src_yend  <- min(src_dim[2] - 1L, src_row_range[2])
  src_xsize <- src_xend - src_xoff + 1L
  src_ysize <- src_yend - src_yoff + 1L

  if (src_xsize <= 0 || src_ysize <= 0) {
    cat(sprintf("  [%d/%d] z=%d x=%d y=%d — empty source window, skipping\n",
                i, nrow(tiles), zoom, tx, ty))
    next
  }

  ## Read source pixels
  src_window <- ds$read(band = 1,
    xoff = src_xoff, yoff = src_yoff,
    xsize = src_xsize, ysize = src_ysize,
    out_xsize = src_xsize, out_ysize = src_ysize)

  ## Warp
  t0 <- proc.time()
  tile_pixels <- rust_warp_scanline(
    src_crs = src_crs, src_gt = src_gt,
    dst_crs = dest_crs, dst_gt = dest_gt,
    dst_dim = c(tile_size, tile_size),
    src_pixels = src_window,
    src_ncol = as.integer(src_xsize), src_nrow = as.integer(src_ysize),
    src_col_off = as.integer(src_xoff), src_row_off = as.integer(src_yoff),
    nodata = as.integer(nodata)
  )
  t_warp <- (proc.time() - t0)[3]

  ## Write output tile
  tile_path <- file.path(out_dir, sprintf("z%d_x%d_y%d.tif", zoom, tx, ty))

  create(
    format = "GTiff",
    dst_filename = tile_path,
    xsize = tile_size, ysize = tile_size,
    nbands = 1, dataType = "Int16",
    options = c("COMPRESS=DEFLATE", "TILED=YES")
  )
  tile_rds <- new(GDALRaster, tile_path, read_only = FALSE)
  tile_rds$setProjection(dest_crs)
  tile_rds$setGeoTransform(dest_gt)
  tile_rds$setNoDataValue(band = 1, as.double(nodata))
  tile_rds$write(band = 1, xoff = 0, yoff = 0,
                 xsize = tile_size, ysize = tile_size,
                 raster_data = tile_pixels)
  tile_rds$close()

  n_valid <- sum(tile_pixels != nodata)
  cat(sprintf("  [%d/%d] z=%d x=%d y=%d — %d/%d valid px, src %dx%d, warp %.3fs\n",
              i, nrow(tiles), zoom, tx, ty, n_valid, tile_size^2,
              src_xsize, src_ysize, t_warp))
}

t_total <- (proc.time() - t_total)[3]

cat(sprintf("\n=== Done: %d tiles in %.1fs (%.3fs per tile) ===\n",
            nrow(tiles), t_total, t_total / nrow(tiles)))

## ---- Summary ----
cat(sprintf("\nOutput: %s/\n", normalizePath(out_dir)))
cat(sprintf("Format: GeoTIFF, EPSG:3857, %dx%d, Int16, DEFLATE\n",
            tile_size, tile_size))
cat(sprintf("Pixel resolution: %.1f m\n", tile_width / tile_size))

ds$close()

## ---- Quick check: mosaic a few tiles with terra ----
cat("\n=== Quick visual check ===\n")
tifs <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(tifs) > 0) {
  cat(sprintf("  %d tiles written\n", length(tifs)))
  cat("  library(terra); r <- vrt(list.files('tiles/tiles_mawson', full=TRUE)); plot(r)\n")
}
