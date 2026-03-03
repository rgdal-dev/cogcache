## Source X distribution for tile [5,3] — the antimeridian-crosser
##
## Questions:
## 1. Is the distribution cleanly bimodal?
## 2. Where's the gap?
## 3. How many pixels are on each side?
## 4. What do the two source windows look like?

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
tx <- 5; ty <- 3
dst_col_off <- (tx - 1L) * tile_size
dst_row_off <- (ty - 1L) * tile_size
tile_xmin <- out_xmin + dst_col_off * pixel_res
tile_xmax <- tile_xmin + tile_size * pixel_res
tile_ymax <- (out_ymin + ceiling(3000000 / pixel_res) * pixel_res) - dst_row_off * pixel_res
tile_ymin <- tile_ymax - tile_size * pixel_res
tile_ext  <- c(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
tile_gt   <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))

## Transform all destination pixels to source
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size)
)

valid <- !is.na(coords$src_x) & !is.na(coords$src_y)
src_x <- coords$src_x[valid]
src_y <- coords$src_y[valid]

cat(sprintf("Tile [5,3]: %d valid of %d pixels\n", sum(valid), length(valid)))
cat(sprintf("Source X range: [%.1f, %.1f]  span=%.0f  (of %d total)\n",
            min(src_x), max(src_x), max(src_x) - min(src_x), src_dim[1]))
cat(sprintf("Source Y range: [%.1f, %.1f]\n", min(src_y), max(src_y)))

## Histogram of source X
cat("\n=== Source X distribution ===\n")
breaks <- seq(0, src_dim[1], length.out = 101)
h <- hist(src_x, breaks = breaks, plot = FALSE)
## Find the gap: largest empty run of bins
empty_runs <- rle(h$counts == 0)
if (any(empty_runs$values)) {
  run_lengths <- empty_runs$lengths[empty_runs$values]
  run_starts <- cumsum(c(1, empty_runs$lengths))[which(empty_runs$values)]
  biggest_gap_idx <- which.max(run_lengths)
  gap_start_bin <- run_starts[biggest_gap_idx]
  gap_end_bin <- gap_start_bin + run_lengths[biggest_gap_idx] - 1
  gap_start_col <- breaks[gap_start_bin]
  gap_end_col <- breaks[gap_end_bin + 1]
  gap_start_lon <- src_gt[1] + gap_start_col * src_gt[2]
  gap_end_lon <- src_gt[1] + gap_end_col * src_gt[2]

  cat(sprintf("Largest gap: columns %.0f to %.0f  (%.0f columns)\n",
              gap_start_col, gap_end_col, gap_end_col - gap_start_col))
  cat(sprintf("  = lon %.2f to %.2f\n", gap_start_lon, gap_end_lon))
  cat(sprintf("  = %.1f%% of source width\n",
              100 * (gap_end_col - gap_start_col) / src_dim[1]))

  ## Split at the gap midpoint
  split_col <- (gap_start_col + gap_end_col) / 2

  west_mask <- src_x < split_col
  east_mask <- src_x >= split_col

  cat(sprintf("\nSplit at column %.0f:\n", split_col))
  cat(sprintf("  West strip (near col 0):     %d pixels, X range [%.1f, %.1f]\n",
              sum(west_mask), min(src_x[west_mask]), max(src_x[west_mask])))
  cat(sprintf("  East strip (near col 86400): %d pixels, X range [%.1f, %.1f]\n",
              sum(east_mask), min(src_x[east_mask]), max(src_x[east_mask])))

  ## What source windows would we need?
  west_xmin <- floor(min(src_x[west_mask]))
  west_xmax <- ceiling(max(src_x[west_mask]))
  east_xmin <- floor(min(src_x[east_mask]))
  east_xmax <- ceiling(max(src_x[east_mask]))
  ymin <- floor(min(src_y))
  ymax <- ceiling(max(src_y))

  cat(sprintf("\n  West source window: xoff=%d xsize=%d  (lon %.1f to %.1f)\n",
              west_xmin, west_xmax - west_xmin,
              src_gt[1] + west_xmin * src_gt[2],
              src_gt[1] + west_xmax * src_gt[2]))
  cat(sprintf("  East source window: xoff=%d xsize=%d  (lon %.1f to %.1f)\n",
              east_xmin, east_xmax - east_xmin,
              src_gt[1] + east_xmin * src_gt[2],
              src_gt[1] + east_xmax * src_gt[2]))
  cat(sprintf("  Y window: yoff=%d ysize=%d\n", ymin, ymax - ymin))

  total_pixels <- (west_xmax - west_xmin + east_xmax - east_xmin) * (ymax - ymin)
  full_width_pixels <- src_dim[1] * (ymax - ymin)
  cat(sprintf("\n  Split read: %d pixels (%.1f%% of full-width read of %d)\n",
              total_pixels, 100 * total_pixels / full_width_pixels, full_width_pixels))
} else {
  cat("No gap found - source X is continuous\n")
}

## Plot
png("./antimeridian_src_x.png", width = 900, height = 600)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

## Top: histogram of source X (column indices)
hist(src_x, breaks = 200, col = "steelblue", border = NA,
     main = "Source X distribution for tile [5,3] (antimeridian crosser)",
     xlab = "GEBCO source column index", ylab = "Destination pixel count")
if (exists("split_col")) {
  abline(v = split_col, col = "red", lwd = 2, lty = 2)
  text(split_col, par("usr")[4] * 0.9, "split", col = "red", pos = 4)
}

## Bottom: source X vs destination column (shows the wrap structure)
## Reshape to matrix
src_x_mat <- matrix(NA_real_, nrow = tile_size, ncol = tile_size)
src_x_mat[valid] <- coords$src_x[valid]

## Plot a few rows to show the wrap
plot(NA, xlim = c(0, 256), ylim = c(0, src_dim[1]),
     xlab = "Destination column", ylab = "Source X (column)",
     main = "Source X vs destination column (selected rows)")
row_colors <- hcl.colors(8, "Set2")
for (i in 1:8) {
  row_idx <- (i - 1) * 32 + 1
  row_data <- src_x_mat[row_idx, ]
  points(1:256, row_data, pch = ".", col = row_colors[i], cex = 2)
}
if (exists("split_col")) {
  abline(h = split_col, col = "red", lwd = 2, lty = 2)
}
legend("right", paste("row", (0:7) * 32), col = row_colors,
       pch = 15, cex = 0.8, bg = "white")

dev.off()
cat("\nPlot saved to antimeridian_src_x.png\n")

ds$close()
cat("Done.\n")
