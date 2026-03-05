## Test GCP polynomial transformer against GDAL
##
## Fabricate GCPs from a Sentinel-2 tile (55GDN, southwest Tasmania),
## compare our Rust polynomial transformer against ground truth
## (exact geotransform + PROJ).

library(gdalraster)

dsn <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/DN/2026/2/S2C_T55GDN_20260227T000650_L2A/B04.tif"

ds <- new(GDALRaster, dsn, TRUE)
src_gt  <- ds$getGeoTransform()
src_crs <- ds$getProjectionRef()
ncol    <- ds$getRasterXSize()
nrow    <- ds$getRasterYSize()
ds$close()

cat("Source tile: 55GDN (southwest Tasmania)\n")
cat("Dimensions:", ncol, "x", nrow, "\n")
cat("Pixel size:", src_gt[2], "x", src_gt[6], "\n\n")

## --- Fabricate GCPs ---
set.seed(42)
grid_cols <- seq(10, ncol - 10, length.out = 5)
grid_rows <- seq(10, nrow - 10, length.out = 5)
gcp_grid <- expand.grid(col = grid_cols, row = grid_rows)
n_gcp <- nrow(gcp_grid)
gcp_grid$col <- gcp_grid$col + runif(n_gcp, -5, 5)
gcp_grid$row <- gcp_grid$row + runif(n_gcp, -5, 5)

gcp_grid$x_proj <- src_gt[1] + gcp_grid$col * src_gt[2] + gcp_grid$row * src_gt[3]
gcp_grid$y_proj <- src_gt[4] + gcp_grid$col * src_gt[5] + gcp_grid$row * src_gt[6]

xy_proj <- cbind(gcp_grid$x_proj, gcp_grid$y_proj)
lonlat <- transform_xy(xy_proj, src_crs, srs_to_wkt("EPSG:4326"))
gcp_grid$lon <- lonlat[, 1]
gcp_grid$lat <- lonlat[, 2]

cat(sprintf("Fabricated %d GCPs (pixel/line -> lon/lat)\n", n_gcp))
cat(sprintf("  lon range: [%.4f, %.4f]\n", min(gcp_grid$lon), max(gcp_grid$lon)))
cat(sprintf("  lat range: [%.4f, %.4f]\n\n", min(gcp_grid$lat), max(gcp_grid$lat)))

## --- Test points (not in GCP set) ---
test_pts <- data.frame(
  name = c("centre", "TL_quarter", "BR_quarter", "random1", "random2"),
  col = c(ncol/2, ncol/4, 3*ncol/4, 1234, 8765),
  row = c(nrow/2, nrow/4, 3*nrow/4, 5678, 2345)
)

test_pts$x_true <- src_gt[1] + test_pts$col * src_gt[2] + test_pts$row * src_gt[3]
test_pts$y_true <- src_gt[4] + test_pts$col * src_gt[5] + test_pts$row * src_gt[6]
true_lonlat <- transform_xy(
  cbind(test_pts$x_true, test_pts$y_true), src_crs, srs_to_wkt("EPSG:4326")
)
test_pts$lon_true <- true_lonlat[, 1]
test_pts$lat_true <- true_lonlat[, 2]

## --- Test Rust GCP transformer ---
cat("=== Rust GCP transformer (forward: pixel/line -> lon/lat) ===\n\n")

for (ord in 1:3) {
  result <- rust_gcp_transform_fwd(
    gcp_pixel = gcp_grid$col,
    gcp_line  = gcp_grid$row,
    gcp_geo_x = gcp_grid$lon,
    gcp_geo_y = gcp_grid$lat,
    order     = ord,
    eval_pixel = test_pts$col,
    eval_line  = test_pts$row
  )

  cat(sprintf("Order %d (fitted as %d):\n", ord, result$order))
  for (i in seq_len(nrow(test_pts))) {
    lon_err <- result$x[i] - test_pts$lon_true[i]
    lat_err <- result$y[i] - test_pts$lat_true[i]
    ## Error in metres (approximate, at mid-latitudes)
    lon_err_m <- lon_err * 111320 * cos(test_pts$lat_true[i] * pi / 180)
    lat_err_m <- lat_err * 110540
    err_m <- sqrt(lon_err_m^2 + lat_err_m^2)
    cat(sprintf("  %-12s lon_err=%+.2e  lat_err=%+.2e  (~%.1f m)\n",
                test_pts$name[i], lon_err, lat_err, err_m))
  }
  cat("\n")
}

## --- Roundtrip test: forward then inverse ---
cat("=== Roundtrip: pixel -> lon/lat -> pixel (order 2) ===\n\n")

fwd <- rust_gcp_transform_fwd(
  gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat,
  order = 2L,
  eval_pixel = test_pts$col, eval_line = test_pts$row
)

inv <- rust_gcp_transform_inv(
  gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat,
  order = 2L,
  eval_x = fwd$x, eval_y = fwd$y
)

for (i in seq_len(nrow(test_pts))) {
  col_err <- inv$pixel[i] - test_pts$col[i]
  row_err <- inv$line[i] - test_pts$row[i]
  cat(sprintf("  %-12s col_err=%+.4f  row_err=%+.4f pixels\n",
              test_pts$name[i], col_err, row_err))
}

## --- Inspect coefficients ---
cat("\n=== Polynomial coefficients (order 2) ===\n\n")
coeffs <- rust_gcp_coefficients(
  gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat,
  order = 2L
)
cat("Forward (pixel/line -> lon):", round(coeffs$to_geo_x, 10), "\n")
cat("Forward (pixel/line -> lat):", round(coeffs$to_geo_y, 10), "\n")
cat("Reverse (lon/lat -> pixel):", round(coeffs$from_geo_pixel, 6), "\n")
cat("Reverse (lon/lat -> line): ", round(coeffs$from_geo_line, 6), "\n")

cat("\nDone.\n")
