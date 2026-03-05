## Test GCP polynomial transformer against GDAL
##
## Fabricate GCPs from a Sentinel-2 tile's known geotransform,
## then compare our Rust polynomial transformer against GDAL's
## using the vrt:// GCP URI syntax.

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
## 5x5 grid with jitter, mapping pixel/line -> lon/lat
set.seed(42)
grid_cols <- seq(10, ncol - 10, length.out = 5)
grid_rows <- seq(10, nrow - 10, length.out = 5)
gcp_grid <- expand.grid(col = grid_cols, row = grid_rows)
n_gcp <- nrow(gcp_grid)
gcp_grid$col <- gcp_grid$col + runif(n_gcp, -5, 5)
gcp_grid$row <- gcp_grid$row + runif(n_gcp, -5, 5)

## Pixel/line -> projected coords via geotransform
gcp_grid$x_proj <- src_gt[1] + gcp_grid$col * src_gt[2] + gcp_grid$row * src_gt[3]
gcp_grid$y_proj <- src_gt[4] + gcp_grid$col * src_gt[5] + gcp_grid$row * src_gt[6]

## Projected -> WGS84 lon/lat
xy_proj <- cbind(gcp_grid$x_proj, gcp_grid$y_proj)
lonlat <- transform_xy(xy_proj, src_crs, srs_to_wkt("EPSG:4326"))
gcp_grid$lon <- lonlat[, 1]
gcp_grid$lat <- lonlat[, 2]

cat(sprintf("Fabricated %d GCPs\n", n_gcp))
cat(sprintf("  pixel range: [%.1f, %.1f]\n", min(gcp_grid$col), max(gcp_grid$col)))
cat(sprintf("  line range:  [%.1f, %.1f]\n", min(gcp_grid$row), max(gcp_grid$row)))
cat(sprintf("  lon range:   [%.4f, %.4f]\n", min(gcp_grid$lon), max(gcp_grid$lon)))
cat(sprintf("  lat range:   [%.4f, %.4f]\n\n", min(gcp_grid$lat), max(gcp_grid$lat)))

## --- Test points (not in GCP set) ---
test_pts <- data.frame(
  name = c("centre", "TL_quarter", "BR_quarter", "random1", "random2"),
  col = c(ncol/2, ncol/4, 3*ncol/4, 1234, 8765),
  row = c(nrow/2, nrow/4, 3*nrow/4, 5678, 2345)
)

## Ground truth: exact geotransform + PROJ
test_pts$x_true <- src_gt[1] + test_pts$col * src_gt[2] + test_pts$row * src_gt[3]
test_pts$y_true <- src_gt[4] + test_pts$col * src_gt[5] + test_pts$row * src_gt[6]
true_lonlat <- transform_xy(
  cbind(test_pts$x_true, test_pts$y_true), src_crs, srs_to_wkt("EPSG:4326")
)
test_pts$lon_true <- true_lonlat[, 1]
test_pts$lat_true <- true_lonlat[, 2]

## --- Test our Rust GCP transformer ---
cat("=== Rust GCP transformer ===\n\n")

for (ord in 1:3) {
  cat(sprintf("Order %d:\n", ord))

  ## Create transformer (this will be rust_create_gcp_transformer)
  ## For now, placeholder showing the interface:
  ## t <- rust_create_gcp_transformer(
  ##   gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat, ord
  ## )

  ## Forward transform test points
  ## result <- rust_gcp_transform(t, test_pts$col, test_pts$row, inverse = FALSE)

  ## For each test point, compute error in lon/lat
  ## for (i in seq_len(nrow(test_pts))) {
  ##   lon_err <- result$x[i] - test_pts$lon_true[i]
  ##   lat_err <- result$y[i] - test_pts$lat_true[i]
  ##   cat(sprintf("  %s: lon_err=%.2e lat_err=%.2e\n",
  ##               test_pts$name[i], lon_err, lat_err))
  ## }

  cat("  (awaiting Rust implementation)\n\n")
}

## --- Validate against GDAL via VRT URI ---
cat("=== GDAL GCP transformer (via VRT) ===\n\n")

## Build VRT URI with GCPs
gcp_params <- paste0(
  sprintf("gcp=%.6f,%.6f,%.10f,%.10f",
          gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat),
  collapse = "&"
)
vrt_uri <- paste0("vrt://", dsn, "?", gcp_params)

cat("VRT URI length:", nchar(vrt_uri), "chars\n")

## Check that GDAL can open it
tryCatch({
  ds2 <- new(GDALRaster, vrt_uri, TRUE)
  n_gcps <- length(vapour::vapour_raster_gcp(vrt_uri)[[1]])
  cat(sprintf("GDAL opened VRT with %d GCPs\n", n_gcps))
  gcp_crs <-  vapour::vapour_raster_gcp(vrt_uri)[["CRS"]]
  cat("GCP CRS:", substr(gcp_crs, 1, 60), "...\n")
  ds2$close()
}, error = function(e) {
  cat("Error opening VRT:", conditionMessage(e), "\n")
  cat("Trying with a_srs parameter...\n")
  vrt_uri2 <- paste0(vrt_uri, "&a_srs=EPSG:4326")
  ds2 <- new(GDALRaster, vrt_uri2, TRUE)
  n_gcps <- length(vapour::vapour_raster_gcp(vrt_uri2)[[1]])
  cat(sprintf("GDAL opened VRT with %d GCPs (with a_srs)\n", n_gcps))
  ds2$close()
})

cat("\nDone. Ready for Rust GCP transformer validation.\n")
