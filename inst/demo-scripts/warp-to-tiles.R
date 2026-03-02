## Warp-to-tiles workflow
##
## 1. Get footprint of source COG
## 2. Find all web tiles (z=10) covering the footprint
## 3. Warp each tile using our Rust pipeline
## 4. Save as individual GeoTIFFs in a tiles/ directory
##
## This demonstrates that the OUTPUT GRID defines the job,
## but the source footprint shapes which tiles we need.

library(gdalraster)
library(reproj)
library(vaster)
library(slippymath)
library(sf)

## ---- Config ----
cog_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"
zoom <- 10L
tile_size <- 256L
out_dir <- "tiles/tiles_tasmania"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Step 1: Source footprint ----
cat("=== Step 1: Compute source footprint ===\n")

ds <- new(GDALRaster, cog_url)
src_crs <- ds$getProjectionRef()
src_gt  <- ds$getGeoTransform()
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
cat(sprintf("  Source: %d x %d pixels\n", src_dim[1], src_dim[2]))

## Use gdal_footprint to get the data extent (handles nodata edges)
fp_file <- tempfile(fileext = ".geojson")
## gdal raster footprint is GDAL >= 3.8, fall back to bbox if not available
fp_ok <- tryCatch({
  system2("gdal", c("raster", "footprint", cog_url, fp_file),
          stdout = TRUE, stderr = TRUE)
  TRUE
}, error = function(e) FALSE)

if (!fp_ok || !file.exists(fp_file) || file.size(fp_file) == 0) {
  cat("  gdal raster footprint not available, using bbox\n")
  bbox <- ds$bbox()  # xmin, ymin, xmax, ymax in src CRS
  corners <- rbind(
    c(bbox[1], bbox[2]), c(bbox[3], bbox[2]),
    c(bbox[3], bbox[4]), c(bbox[1], bbox[4]),
    c(bbox[1], bbox[2])
  )
  ## Transform to 4326
  corners_ll <- transform_xy(corners, srs_from = src_crs,
                             srs_to = srs_to_wkt("EPSG:4326"))
  footprint <- st_sfc(st_polygon(list(corners_ll)), crs = 4326)
} else {
  footprint_sf <- st_read(fp_file, quiet = TRUE)
  footprint <- st_geometry(footprint_sf)
  if (st_crs(footprint)$epsg != 4326) {
    footprint <- st_transform(footprint, 4326)
  }
}

fp_bbox <- st_bbox(footprint)
cat(sprintf("  Footprint bbox (EPSG:4326): [%.4f, %.4f, %.4f, %.4f]\n",
            fp_bbox["xmin"], fp_bbox["ymin"], fp_bbox["xmax"], fp_bbox["ymax"]))

## ---- Step 2: Find covering tiles ----
cat("\n=== Step 2: Find web tiles at zoom %d ===\n", zoom)

## Get tile range from the footprint bbox
tile_min <- lonlat_to_tilenum(fp_bbox["xmin"], fp_bbox["ymax"], zoom)
tile_max <- lonlat_to_tilenum(fp_bbox["xmax"], fp_bbox["ymin"], zoom)

x_range <- tile_min$x:tile_max$x
y_range <- tile_min$y:tile_max$y

## Build tile grid and filter by footprint intersection
tile_extent_3857 <- function(z, x, y) {
  n <- 2^z
  lon <- c(x, x + 1) / n * 360 - 180
  lat <- atan(sinh(pi * (1 - 2 * c(y, y + 1) / n))) * 180 / pi
  reproj_extent(c(lon, rev(lat)), "EPSG:3857", source = "EPSG:4326")
}

tile_bbox_4326 <- function(z, x, y) {
  n <- 2^z
  lon <- c(x, x + 1) / n * 360 - 180
  lat <- atan(sinh(pi * (1 - 2 * c(y, y + 1) / n))) * 180 / pi
  c(xmin = lon[1], ymin = lat[2], xmax = lon[2], ymax = lat[1])
}

tiles <- expand.grid(x = x_range, y = y_range)
cat(sprintf("  Tile range: x=%d..%d, y=%d..%d (%d candidates)\n",
            min(x_range), max(x_range), min(y_range), max(y_range), nrow(tiles)))

## Filter tiles that intersect the footprint
keep <- logical(nrow(tiles))
for (i in seq_len(nrow(tiles))) {
  tb <- tile_bbox_4326(zoom, tiles$x[i], tiles$y[i])
  tile_poly <- st_sfc(st_polygon(list(rbind(
    c(tb["xmin"], tb["ymin"]), c(tb["xmax"], tb["ymin"]),
    c(tb["xmax"], tb["ymax"]), c(tb["xmin"], tb["ymax"]),
    c(tb["xmin"], tb["ymin"])
  ))), crs = 4326)
  keep[i] <- lengths(st_intersects(tile_poly, footprint)) > 0
}
tiles <- tiles[keep, ]
cat(sprintf("  Tiles intersecting footprint: %d\n", nrow(tiles)))

## ---- Step 3: Read source data (once) ----
cat("\n=== Step 3: Read source data ===\n")

## We need to know the full extent of all tiles in source pixel space
## to read a single source window that covers everything.
## For simplicity, read the relevant portion of the source for each tile.
## But first, let's see if there's a common source window.

nodata <- ds$getNoDataValue(band = 1)
if (is.na(nodata)) nodata <- 0L

dest_crs <- srs_to_wkt("EPSG:3857")

## ---- Step 4: Warp each tile ----
cat(sprintf("\n=== Step 4: Warp %d tiles ===\n\n", nrow(tiles)))

t_total <- proc.time()

for (i in seq_len(nrow(tiles))) {
  tx <- tiles$x[i]
  ty <- tiles$y[i]

  ## Destination grid for this tile
  dest_ext <- tile_extent_3857(zoom, tx, ty)
  dest_gt  <- extent_dim_to_gt(dest_ext, c(tile_size, tile_size))

  ## Compute source window for this tile
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

  tile_ds <- create(
    format = "GTiff",
    dst_filename = tile_path,
    xsize = tile_size, ysize = tile_size,
    nbands = 1, dataType = "UInt16",
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
  cat(sprintf("  [%d/%d] z=%d x=%d y=%d — %d valid px, src %dx%d, warp %.3fs → %s\n",
              i, nrow(tiles), zoom, tx, ty, n_valid,
              src_xsize, src_ysize, t_warp, tile_path))
}

t_total <- (proc.time() - t_total)[3]

cat(sprintf("\n=== Done: %d tiles in %.1fs (%.3fs per tile) ===\n",
            nrow(tiles), t_total, t_total / nrow(tiles)))

## ---- Summary ----
cat(sprintf("\nOutput directory: %s/\n", normalizePath(out_dir)))
cat(sprintf("Tile format: GeoTIFF, EPSG:3857, %dx%d, UInt16, DEFLATE\n",
            tile_size, tile_size))

ds$close()
