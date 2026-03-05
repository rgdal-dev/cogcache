## Validate Rust GCP transformer against GDAL's GCP transformer
##
## Uses the vrt:// GCP URI to create a dataset with GCPs, then
## compares GDAL's polynomial transform (via gdalwarp/gdaltransform)
## against our Rust implementation.

library(gdalraster)

dsn <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/DN/2026/2/S2C_T55GDN_20260227T000650_L2A/B04.tif"

ds <- new(GDALRaster, dsn, TRUE)
src_gt  <- ds$getGeoTransform()
src_crs <- ds$getProjectionRef()
ncol    <- ds$getRasterXSize()
nrow    <- ds$getRasterYSize()
ds$close()

## --- Fabricate GCPs (same as test script, same seed) ---
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

## --- Build VRT URI with GCPs ---
gcp_params <- paste0(
  sprintf("gcp=%.6f,%.6f,%.10f,%.10f",
          gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat),
  collapse = "&"
)
vrt_uri <- paste0("vrt://", dsn, "?a_srs=EPSG:4326&", gcp_params)

cat("VRT URI length:", nchar(vrt_uri), "chars\n\n")

## --- Test points ---
test_pts <- data.frame(
  name = c("centre", "TL_quarter", "BR_quarter", "random1", "random2"),
  col = c(ncol/2, ncol/4, 3*ncol/4, 1234, 8765),
  row = c(nrow/2, nrow/4, 3*nrow/4, 5678, 2345)
)

## --- GDAL's GCP transform via gdaltransform CLI ---
## gdaltransform reads pixel/line from stdin and outputs geo coords
cat("=== Comparing Rust vs GDAL GCP polynomial (order 1, 2, 3) ===\n\n")

for (ord in 1:3) {
  cat(sprintf("--- Order %d ---\n", ord))

  ## Rust forward transform
  rust_fwd <- rust_gcp_transform_fwd(
    gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat,
    order = ord,
    eval_pixel = test_pts$col, eval_line = test_pts$row
  )

  ## GDAL forward transform via gdaltransform
  ## Build input: "col row" per line
  input_text <- paste(test_pts$col, test_pts$row, collapse = "\n")

  gdal_cmd <- sprintf(
    'echo "%s" | gdaltransform -order %d "%s"',
    input_text, ord, vrt_uri
  )

  gdal_output <- system(gdal_cmd, intern = TRUE)

  ## Parse GDAL output: "lon lat elevation" per line
  gdal_parsed <- do.call(rbind, strsplit(trimws(gdal_output), "\\s+"))
  gdal_lon <- as.numeric(gdal_parsed[, 1])
  gdal_lat <- as.numeric(gdal_parsed[, 2])

  ## Compare
  for (i in seq_len(nrow(test_pts))) {
    lon_diff <- rust_fwd$x[i] - gdal_lon[i]
    lat_diff <- rust_fwd$y[i] - gdal_lat[i]
    ## Convert to approximate metres
    lon_diff_m <- lon_diff * 111320 * cos(gdal_lat[i] * pi / 180)
    lat_diff_m <- lat_diff * 110540
    diff_m <- sqrt(lon_diff_m^2 + lat_diff_m^2)
    cat(sprintf("  %-12s rust=(%.8f,%.8f) gdal=(%.8f,%.8f) diff=%.4e m\n",
                test_pts$name[i],
                rust_fwd$x[i], rust_fwd$y[i],
                gdal_lon[i], gdal_lat[i],
                diff_m))
  }
  cat("\n")
}

## --- Also compare inverse (geo -> pixel/line) ---
cat("=== Inverse comparison (order 2): geo -> pixel/line ===\n\n")

## Use the order-2 forward results as input for inverse
rust_fwd2 <- rust_gcp_transform_fwd(
  gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat,
  order = 2L,
  eval_pixel = test_pts$col, eval_line = test_pts$row
)

## Rust inverse
rust_inv <- rust_gcp_transform_inv(
  gcp_grid$col, gcp_grid$row, gcp_grid$lon, gcp_grid$lat,
  order = 2L,
  eval_x = rust_fwd2$x, eval_y = rust_fwd2$y
)

## GDAL inverse via gdaltransform -i
input_geo <- paste(rust_fwd2$x, rust_fwd2$y, collapse = "\n")
gdal_inv_cmd <- sprintf(
  'echo "%s" | gdaltransform -i -order 2 "%s"',
  input_geo, vrt_uri
)
gdal_inv_output <- system(gdal_inv_cmd, intern = TRUE)
gdal_inv_parsed <- do.call(rbind, strsplit(trimws(gdal_inv_output), "\\s+"))
gdal_pixel <- as.numeric(gdal_inv_parsed[, 1])
gdal_line  <- as.numeric(gdal_inv_parsed[, 2])

for (i in seq_len(nrow(test_pts))) {
  pixel_diff <- rust_inv$pixel[i] - gdal_pixel[i]
  line_diff  <- rust_inv$line[i] - gdal_line[i]
  cat(sprintf("  %-12s rust=(%.4f,%.4f) gdal=(%.4f,%.4f) diff=(%.4e,%.4e) px\n",
              test_pts$name[i],
              rust_inv$pixel[i], rust_inv$line[i],
              gdal_pixel[i], gdal_line[i],
              pixel_diff, line_diff))
}

cat("\nDone.\n")
