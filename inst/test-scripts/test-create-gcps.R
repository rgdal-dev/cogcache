## Fabricate GCPs from Sentinel-2 geotransform
##
## Take a known Sentinel-2 tile (UTM projected), sample a grid of
## pixel/line positions, transform to lon/lat via PROJ, and produce
## a GCP table. This gives us ground truth to validate a GCP polynomial
## transformer against the exact geotransform + CRS path.

library(gdalraster)

## Macquarie Island tile from our test suite
s2_url <- "/vsicurl/https://e84-earth-search-sentinel-data.s3.us-west-2.amazonaws.com/sentinel-2-c1-l2a/55/G/CM/2026/2/S2C_T55GCM_20260227T000650_L2A/B04.tif"
ds <- new(GDALRaster, s2_url, TRUE)
src_gt  <- ds$getGeoTransform()
src_crs <- ds$getProjectionRef()
ncol    <- ds$getRasterXSize()
nrow    <- ds$getRasterYSize()
ds$close()

cat("Source CRS:", substr(src_crs, 1, 80), "...\n")
cat("Dimensions:", ncol, "x", nrow, "\n")
cat("GeoTransform:", paste(round(src_gt, 4), collapse = ", "), "\n")
cat("Pixel size:", src_gt[2], "x", src_gt[6], "\n\n")

## Sample a grid of pixel/line positions
## Use irregular spacing to make the polynomial fit non-trivial
set.seed(42)
n_gcp <- 25  # 5x5 grid (enough for 3rd order polynomial which needs 10)
grid_cols <- seq(10, ncol - 10, length.out = 5)
grid_rows <- seq(10, nrow - 10, length.out = 5)
gcp_grid <- expand.grid(col = grid_cols, row = grid_rows)

## Add some jitter to avoid perfectly regular grid
gcp_grid$col <- gcp_grid$col + runif(n_gcp, -5, 5)
gcp_grid$row <- gcp_grid$row + runif(n_gcp, -5, 5)

## Convert pixel/line to projected coordinates via geotransform
gcp_grid$x_proj <- src_gt[1] + gcp_grid$col * src_gt[2] + gcp_grid$row * src_gt[3]
gcp_grid$y_proj <- src_gt[4] + gcp_grid$col * src_gt[5] + gcp_grid$row * src_gt[6]

## Transform projected coords to WGS84 lon/lat
xy_proj <- cbind(gcp_grid$x_proj, gcp_grid$y_proj)
lonlat <- transform_xy(xy_proj, src_crs, srs_to_wkt("EPSG:4326"))

gcp_grid$lon <- lonlat[, 1]
gcp_grid$lat <- lonlat[, 2]

cat("GCP table (pixel, line, lon, lat):\n")
print(gcp_grid[, c("col", "row", "lon", "lat")], digits = 8)

cat("\nGCP extent:\n")
cat("  pixel:", range(gcp_grid$col), "\n")
cat("  line: ", range(gcp_grid$row), "\n")
cat("  lon:  ", range(gcp_grid$lon), "\n")
cat("  lat:  ", range(gcp_grid$lat), "\n")

## Save for use in GCP transformer tests
saveRDS(gcp_grid, "inst/test-data/s2_macquarie_gcps.rds")

## Also save as a simple text format
write.csv(gcp_grid[, c("col", "row", "lon", "lat")],
          "inst/test-data/s2_macquarie_gcps.csv",
          row.names = FALSE)

cat("\nSaved to inst/test-data/s2_macquarie_gcps.{rds,csv}\n")

## Validation: what does GDAL's GCP transformer produce?
## Create a VRT with GCPs instead of a geotransform
cat("\n=== Validation via GDAL GCP warp ===\n")

## Build GCP list for gdal_translate
gcp_args <- character(0)
for (i in seq_len(nrow(gcp_grid))) {
  gcp_args <- c(gcp_args, sprintf("-gcp %.6f %.6f %.10f %.10f",
                                  gcp_grid$col[i], gcp_grid$row[i],
                                  gcp_grid$lon[i], gcp_grid$lat[i]))
}

cat("Number of GCPs:", nrow(gcp_grid), "\n")
cat("Sample GCP arg:", gcp_args[1], "\n")

## Quick test: transform a few known points through GDAL's GCP path
## using gdaltransform
## We'll test the centre pixel
test_col <- ncol / 2
test_row <- nrow / 2

## Expected lon/lat from the exact geotransform + PROJ path
test_x <- src_gt[1] + test_col * src_gt[2] + test_row * src_gt[3]
test_y <- src_gt[4] + test_col * src_gt[5] + test_row * src_gt[6]
test_lonlat <- transform_xy(cbind(test_x, test_y), src_crs, srs_to_wkt("EPSG:4326"))

cat(sprintf("\nCentre pixel (%d, %d):\n", test_col, test_row))
cat(sprintf("  Expected (from GT+PROJ): lon=%.10f lat=%.10f\n",
            test_lonlat[1, 1], test_lonlat[1, 2]))
cat(sprintf("  Projected coords: x=%.4f y=%.4f\n", test_x, test_y))

## Also test corners
corners <- data.frame(
  name = c("TL", "TR", "BL", "BR", "centre"),
  col = c(0, ncol, 0, ncol, ncol/2),
  row = c(0, 0, nrow, nrow, nrow/2)
)
corners$x <- src_gt[1] + corners$col * src_gt[2] + corners$row * src_gt[3]
corners$y <- src_gt[4] + corners$col * src_gt[5] + corners$row * src_gt[6]
corner_lonlat <- transform_xy(cbind(corners$x, corners$y), src_crs,
                              srs_to_wkt("EPSG:4326"))
corners$lon <- corner_lonlat[, 1]
corners$lat <- corner_lonlat[, 2]

cat("\nCorner coordinates (ground truth for validation):\n")
print(corners, digits = 10)

saveRDS(corners, "inst/test-data/s2_macquarie_corners.rds")

cat("\nDone. Files saved to inst/test-data/\n")
