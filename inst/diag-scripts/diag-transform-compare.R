## Diagnostic: compare coordinate transforms directly
##
## Our Rust GenImgProjTransformer vs GDAL's actual transformer output
## at destination pixel [col=1, row=157]
##
## If these differ by ~2-3 source pixels, the bug is in the PROJ pipeline
## selection, not the bilinear kernel.

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

## --- Our Rust transformer ---
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size)
)

## --- GDAL's transformer via pixel_extract ---
## pixel_extract gives us GDAL's own coordinate transform result
## It returns the value at the source pixel, which tells us which pixel GDAL picked

## But we need the actual coordinates. Let's use sf/reproj to do the
## same three-step pipeline with system PROJ:

## Step 1: dst pixel -> dst geo (apply tile geotransform)
## For pixel [col=1, row=157], centre = (1.5, 157.5)
dst_px <- 1.5
dst_py <- 157.5

dst_geo_x <- tile_gt[1] + dst_px * tile_gt[2] + dst_py * tile_gt[3]
dst_geo_y <- tile_gt[4] + dst_px * tile_gt[5] + dst_py * tile_gt[6]
cat(sprintf("Dst pixel (1.5, 157.5) -> dst geo: (%.4f, %.4f)\n", dst_geo_x, dst_geo_y))

## Inverse source geotransform for geo -> pixel conversion
src_inv_gt <- gdalraster::inv_geotransform(src_gt)

## Step 2: dst geo -> src geo (PROJ transform)
## Use sf to do this with system PROJ
if (requireNamespace("sf", quietly = TRUE)) {
  pt <- sf::st_sfc(sf::st_point(c(dst_geo_x, dst_geo_y)), crs = fiji_lcc)
  pt_src <- sf::st_transform(pt, "EPSG:4326")  # GEBCO is WGS84 geographic
  src_geo <- sf::st_coordinates(pt_src)
  cat(sprintf("sf PROJ: dst geo -> src geo (lon, lat): (%.8f, %.8f)\n",
              src_geo[1], src_geo[2]))

  ## Step 3: src geo -> src pixel (apply inverse src geotransform)
  sf_src_px <- src_inv_gt[1] + src_geo[1] * src_inv_gt[2] + src_geo[2] * src_inv_gt[3]
  sf_src_py <- src_inv_gt[4] + src_geo[1] * src_inv_gt[5] + src_geo[2] * src_inv_gt[6]
  cat(sprintf("sf pipeline -> src pixel: (%.6f, %.6f)\n", sf_src_px, sf_src_py))
} else {
  cat("sf not available, skipping sf comparison\n")
}

## Also try reproj if available
if (requireNamespace("reproj", quietly = TRUE)) {
  rp <- reproj::reproj_xy(cbind(dst_geo_x, dst_geo_y),
                          target = "EPSG:4326",
                          source = fiji_lcc)
  cat(sprintf("reproj: dst geo -> src geo (lon, lat): (%.8f, %.8f)\n", rp[1], rp[2]))
  rp_src_px <- src_inv_gt[1] + rp[1] * src_inv_gt[2] + rp[2] * src_inv_gt[3]
  rp_src_py <- src_inv_gt[4] + rp[1] * src_inv_gt[5] + rp[2] * src_inv_gt[6]
  cat(sprintf("reproj pipeline -> src pixel: (%.6f, %.6f)\n", rp_src_px, rp_src_py))
}

## Our Rust result
idx <- 157 * tile_size + 1 + 1
rust_src_x <- coords$src_x[idx]
rust_src_y <- coords$src_y[idx]
cat(sprintf("\nRust GenImgProjTransformer -> src pixel: (%.6f, %.6f)\n",
            rust_src_x, rust_src_y))

## Difference
if (exists("sf_src_px")) {
  cat(sprintf("\nDifference (Rust - sf): dx=%.4f dy=%.4f pixels\n",
              rust_src_x - sf_src_px, rust_src_y - sf_src_py))
}
if (exists("rp_src_px")) {
  cat(sprintf("Difference (Rust - reproj): dx=%.4f dy=%.4f pixels\n",
              rust_src_x - rp_src_px, rust_src_y - rp_src_py))
}

## Now check ALL pixels - get a global picture of the offset
cat("\n=== Global coordinate comparison (Rust vs sf) at sampled pixels ===\n")
if (requireNamespace("sf", quietly = TRUE)) {
  ## Sample 100 pixels across the tile
  set.seed(42)
  sample_cols <- sample(0:255, 10)
  sample_rows <- sample(0:255, 10)

  diffs_x <- c()
  diffs_y <- c()
  for (r in sample_rows) {
    for (c in sample_cols) {
      idx_i <- r * tile_size + c + 1
      rust_x <- coords$src_x[idx_i]
      rust_y <- coords$src_y[idx_i]

      dpx <- c + 0.5
      dpy <- r + 0.5
      dgx <- tile_gt[1] + dpx * tile_gt[2] + dpy * tile_gt[3]
      dgy <- tile_gt[4] + dpx * tile_gt[5] + dpy * tile_gt[6]

      pt_i <- sf::st_sfc(sf::st_point(c(dgx, dgy)), crs = fiji_lcc)
      pt_src_i <- sf::st_transform(pt_i, "EPSG:4326")
      sg <- sf::st_coordinates(pt_src_i)

      sf_px <- src_inv_gt[1] + sg[1] * src_inv_gt[2] + sg[2] * src_inv_gt[3]
      sf_py <- src_inv_gt[4] + sg[1] * src_inv_gt[5] + sg[2] * src_inv_gt[6]

      diffs_x <- c(diffs_x, rust_x - sf_px)
      diffs_y <- c(diffs_y, rust_y - sf_py)
    }
  }
  cat(sprintf("  dx: mean=%.4f  sd=%.4f  min=%.4f  max=%.4f\n",
              mean(diffs_x), sd(diffs_x), min(diffs_x), max(diffs_x)))
  cat(sprintf("  dy: mean=%.4f  sd=%.4f  min=%.4f  max=%.4f\n",
              mean(diffs_y), sd(diffs_y), min(diffs_y), max(diffs_y)))
}

ds$close()
cat("\nDone.\n")
