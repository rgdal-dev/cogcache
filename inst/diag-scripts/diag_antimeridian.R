## Walkthrough: what happens when a destination tile crosses the antimeridian?
##
## Fiji LCC is centred on lon=178. Tiles to the east of centre will have
## destination coordinates that, when inverse-projected to lon/lat, span
## across 180° (e.g. 179°E to 181°E = -179°E).
##
## GEBCO source is a standard -180 to +180 longitude grid.
## Source pixels at 179°E are near column 86160 (of 86400).
## Source pixels at 181°E = -179°E are near column 240.
##
## Question: does compute_source_window return a sensible window,
## or does it span the entire 86400-pixel width?

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

cat(sprintf("GEBCO: %d x %d, pixel size %.6f°\n", src_dim[1], src_dim[2], src_gt[2]))
cat(sprintf("  lon range: %.1f to %.1f\n", src_gt[1], src_gt[1] + src_dim[1] * src_gt[2]))
cat(sprintf("  Column at 178°E: %.0f\n", (178 - src_gt[1]) / src_gt[2]))
cat(sprintf("  Column at 180°E: %.0f\n", (180 - src_gt[1]) / src_gt[2]))
cat(sprintf("  Column at -179°E (=181°E): %.0f\n", (-179 - src_gt[1]) / src_gt[2]))

tile_size <- 256L
pixel_res <- 2000
out_xmin <- -2000000; out_ymin <- -1500000

## Survey a grid of tiles to find which ones cross the antimeridian
cat("\n=== Tile survey: source window for each tile column ===\n")
cat("(Row 3, columns 1-8)\n\n")

ty <- 3
for (tx in 1:8) {
  dst_col_off <- (tx - 1L) * tile_size
  dst_row_off <- (ty - 1L) * tile_size
  tile_xmin <- out_xmin + dst_col_off * pixel_res
  tile_xmax <- tile_xmin + tile_size * pixel_res
  tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) - dst_row_off * pixel_res
  tile_ymin <- tile_ymax - tile_size * pixel_res
  tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
  tile_gt   <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))

  sw <- rust_compute_source_window(
    src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
    dst_crs = out_crs, dst_gt = tile_gt,
    dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
    resample_padding = 1L
  )

  if (is.null(sw)) {
    cat(sprintf("  tile [%d,%d]: NULL (no valid transforms)\n", tx, ty))
    next
  }

  ## Convert source window to lon range
  sw_lon_min <- src_gt[1] + sw$xoff * src_gt[2]
  sw_lon_max <- src_gt[1] + (sw$xoff + sw$xsize) * src_gt[2]

  antimeridian_flag <- if (sw$xsize > 0.9 * src_dim[1]) " *** FULL WIDTH ***" else ""

  cat(sprintf("  tile [%d,%d]: sw xoff=%5d xsize=%5d  lon=[%7.1f, %7.1f]  fill=%.3f%s\n",
              tx, ty, sw$xoff, sw$xsize,
              sw_lon_min, sw_lon_max, sw$fill_ratio, antimeridian_flag))
}

## Now pick a tile that's clearly east of the antimeridian
## and trace the actual coordinates
cat("\n=== Detailed trace for tile [6,3] (should cross antimeridian) ===\n")
tx <- 6; ty <- 3
dst_col_off <- (tx - 1L) * tile_size
dst_row_off <- (ty - 1L) * tile_size
tile_xmin <- out_xmin + dst_col_off * pixel_res
tile_xmax <- tile_xmin + tile_size * pixel_res
tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) - dst_row_off * pixel_res
tile_ymin <- tile_ymax - tile_size * pixel_res
tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
tile_gt   <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))

cat(sprintf("  Destination extent: [%.0f, %.0f] x [%.0f, %.0f] metres\n",
            tile_xmin, tile_xmax, tile_ymin, tile_ymax))

## Transform the 4 corners to lon/lat
corners <- list(
  c(0.5, 0.5),           # top-left
  c(255.5, 0.5),         # top-right
  c(0.5, 255.5),         # bottom-left
  c(255.5, 255.5)        # bottom-right
)

cat("  Corner coordinates (dst pixel -> lon/lat):\n")
for (cn in seq_along(corners)) {
  px <- corners[[cn]]
  geo_x <- tile_gt[1] + px[1] * tile_gt[2] + px[2] * tile_gt[3]
  geo_y <- tile_gt[4] + px[1] * tile_gt[5] + px[2] * tile_gt[6]
  pt <- sf::st_sfc(sf::st_point(c(geo_x, geo_y)), crs = fiji_lcc)
  ll <- sf::st_coordinates(sf::st_transform(pt, "EPSG:4326"))
  cat(sprintf("    [%.1f, %.1f] -> (%.4f, %.4f)  %s\n",
              px[1], px[2], ll[1], ll[2],
              c("top-left", "top-right", "bottom-left", "bottom-right")[cn]))
}

## Transform all pixels and show the source coord distribution
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size)
)

valid <- !is.na(coords$src_x) & !is.na(coords$src_y)
cat(sprintf("\n  Valid transforms: %d/%d\n", sum(valid), length(valid)))
cat(sprintf("  Source X range: [%.1f, %.1f]  (GEBCO columns)\n",
            min(coords$src_x[valid]), max(coords$src_x[valid])))
cat(sprintf("  Source Y range: [%.1f, %.1f]  (GEBCO rows)\n",
            min(coords$src_y[valid]), max(coords$src_y[valid])))

## Convert source pixel coords back to lon/lat to check wrapping
src_lons <- src_gt[1] + coords$src_x[valid] * src_gt[2]
cat(sprintf("  Source lon range: [%.4f, %.4f]\n", min(src_lons), max(src_lons)))

## Is there a bimodal distribution (some near 180, some near -180)?
n_east <- sum(src_lons > 170)
n_west <- sum(src_lons < -170)
n_middle <- sum(src_lons >= -170 & src_lons <= 170)
cat(sprintf("  Lon distribution: %d east of 170°, %d west of -170°, %d in middle\n",
            n_east, n_west, n_middle))

## Source window for this tile
sw <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
if (!is.null(sw)) {
  sw_lon_min <- src_gt[1] + sw$xoff * src_gt[2]
  sw_lon_max <- src_gt[1] + (sw$xoff + sw$xsize) * src_gt[2]
  cat(sprintf("\n  Source window: xoff=%d xsize=%d  lon=[%.4f, %.4f]\n",
              sw$xoff, sw$xsize, sw_lon_min, sw_lon_max))
  cat(sprintf("  That's %.1f%% of GEBCO's width (%d / %d)\n",
              100 * sw$xsize / src_dim[1], sw$xsize, src_dim[1]))
  cat(sprintf("  fill_ratio=%.3f  n_failed=%d\n", sw$fill_ratio, sw$n_failed))

  if (sw$xsize > 0.5 * src_dim[1]) {
    cat("\n  *** This reads more than half the globe! ***\n")
    cat("  The antimeridian heuristic snapped to full width because\n")
    cat("  source coords span from near column 0 to near column 86400.\n")
    cat("  PROJ wraps 181°E to -179°, which maps to GEBCO column ~240.\n")
    cat("  The bounding box of [240, 86160] in source pixels spans\n")
    cat("  almost the entire raster, triggering the 90% heuristic.\n")
  }
}

ds$close()
cat("\nDone.\n")
