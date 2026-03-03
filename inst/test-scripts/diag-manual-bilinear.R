## Diagnostic: manually compute bilinear at [1,157] and compare with GDAL
##
## We know: buf_x=394.9785, buf_y=980.3255 (from exact transform with read_xoff offset)
## Bilinear: ix = floor(394.9785 - 0.5) = 394, iy = floor(980.3255 - 0.5) = 979
## rx = 1.5 - (394.9785 - 394) = 1.5 - 0.9785 = 0.5215
## ry = 1.5 - (980.3255 - 979) = 1.5 - 1.3255 = 0.1745

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
nodata <- -32768L

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

## Source window + read
sw <- rust_compute_source_window(
  src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 1L
)
src_block <- ds$getBlockSize(band = 1)
src_tw <- src_block[1]; src_th <- src_block[2]
read_xoff  <- (sw$xoff %/% src_tw) * src_tw
read_yoff  <- (sw$yoff %/% src_th) * src_th
read_xend  <- min(ceiling((sw$xoff + sw$xsize) / src_tw) * src_tw, src_dim[1])
read_yend  <- min(ceiling((sw$yoff + sw$ysize) / src_th) * src_th, src_dim[2])
read_xsize <- as.integer(read_xend - read_xoff)
read_ysize <- as.integer(read_yend - read_yoff)
src_pixels <- ds$read(band = 1, xoff = read_xoff, yoff = read_yoff,
                      xsize = read_xsize, ysize = read_ysize,
                      out_xsize = read_xsize, out_ysize = read_ysize)

## Get exact coordinates at [1,157]
coords <- rust_gen_img_proj_transform(
  src_crs = src_crs, src_gt = src_gt,
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size)
)
idx <- 157 * tile_size + 1 + 1  # row=157, col=1, 1-based
src_x_full <- coords$src_x[idx]
src_y_full <- coords$src_y[idx]

buf_x <- src_x_full - read_xoff
buf_y <- src_y_full - read_yoff

cat(sprintf("Full-image: src_x=%.6f src_y=%.6f\n", src_x_full, src_y_full))
cat(sprintf("Buffer-relative: buf_x=%.6f buf_y=%.6f\n", buf_x, buf_y))

## Bilinear kernel computation (matching warp.rs exactly)
ix <- floor(buf_x - 0.5)
iy <- floor(buf_y - 0.5)
rx <- 1.5 - (buf_x - ix)
ry <- 1.5 - (buf_y - iy)

cat(sprintf("\nBilinear kernel:\n"))
cat(sprintf("  ix=%d iy=%d\n", ix, iy))
cat(sprintf("  rx=%.6f ry=%.6f  (1-rx=%.6f 1-ry=%.6f)\n", rx, ry, 1-rx, 1-ry))

## 1-based R indexing into the row-major src_pixels vector
## src_pixels is row-major: index = row * ncol + col + 1 (1-based)
p00 <- src_pixels[iy     * read_xsize + ix     + 1]  # top-left
p10 <- src_pixels[iy     * read_xsize + (ix+1) + 1]  # top-right
p01 <- src_pixels[(iy+1) * read_xsize + ix     + 1]  # bottom-left
p11 <- src_pixels[(iy+1) * read_xsize + (ix+1) + 1]  # bottom-right

cat(sprintf("\nSource pixels (2x2 neighbourhood):\n"))
cat(sprintf("  [%d,%d]=%d  [%d,%d]=%d\n", ix, iy, p00, ix+1, iy, p10))
cat(sprintf("  [%d,%d]=%d  [%d,%d]=%d\n", ix, iy+1, p01, ix+1, iy+1, p11))

## Manual bilinear (matching our Rust formula line 231-233)
val_rust <- (p00 * rx + p10 * (1 - rx)) * ry +
  (p01 * rx + p11 * (1 - rx)) * (1 - ry)

## GDAL-style round
rust_round <- function(v) ifelse(v >= 0, as.integer(v + 0.5), as.integer(v - 0.5))
our_result <- rust_round(val_rust)

cat(sprintf("\nManual bilinear: %.4f -> rounded: %d\n", val_rust, our_result))

## What does our Rust warp actually produce?
our <- rust_warp_resample(
  src_crs = src_crs, src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
  src_ncol = read_xsize, src_nrow = read_ysize,
  src_col_off = as.integer(read_xoff), src_row_off = as.integer(read_yoff),
  nodata = nodata, max_error = 0.0, resample = "bilinear"
)

## What does GDAL produce?
gdal_file <- tempfile(fileext = ".tif")
warp(cog_url, gdal_file, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "bilinear", "-ovr", "NONE", "-et", "0"))
gds <- new(GDALRaster, gdal_file)
gdal_px <- gds$read(band = 1, xoff = 0, yoff = 0,
                    xsize = tile_size, ysize = tile_size,
                    out_xsize = tile_size, out_ysize = tile_size)
gds$close()

cat(sprintf("\nResults at [1,157]:\n"))
cat(sprintf("  Manual bilinear in R:  %d\n", our_result))
cat(sprintf("  Rust warp_resample:    %d\n", our[idx]))
cat(sprintf("  GDAL warp:             %d\n", gdal_px[idx]))
cat(sprintf("  diff (rust - gdal):    %d\n", our[idx] - gdal_px[idx]))

## What gradient exists here?
cat(sprintf("\nLocal gradient (pixels around the neighbourhood):\n"))
for (dy in -2:3) {
  vals <- sapply(-2:3, function(dx) {
    src_pixels[(iy + dy) * read_xsize + (ix + dx) + 1]
  })
  cat(sprintf("  row %d: %s\n", iy + dy, paste(vals, collapse = "  ")))
}

## Now check: what source coordinate does GDAL think this pixel maps to?
## Use gdalwarp -tap and check the inverse
cat(sprintf("\nFor GDAL to produce %d, what source coord would give that?\n", gdal_px[idx]))
cat("If GDAL bilinear at SAME coords gives different result, the formula differs.\n")
cat("If GDAL bilinear at DIFFERENT coords gives -4030, the transform differs.\n")

## Try manually: what buf_x would produce -4030 from these 4 pixels?
## val = (p00*rx + p10*(1-rx))*ry + (p01*rx + p11*(1-rx))*(1-ry)
## With p00, p10, p01, p11 known, solve for rx, ry that give -4030
## This is underdetermined (2 unknowns, 1 equation) but we can check
## if a small perturbation in buf_x explains it
cat(sprintf("\nSweep buf_x around %.4f to find what gives GDAL's %d:\n", buf_x, gdal_px[idx]))
for (delta in seq(-2, 2, by = 0.1)) {
  test_bx <- buf_x + delta
  test_ix <- floor(test_bx - 0.5)
  test_rx <- 1.5 - (test_bx - test_ix)
  ## recalculate with potentially different pixels if ix changed
  t_ix <- as.integer(test_ix)
  t_p00 <- src_pixels[iy * read_xsize + t_ix + 1]
  t_p10 <- src_pixels[iy * read_xsize + (t_ix+1) + 1]
  t_p01 <- src_pixels[(iy+1) * read_xsize + t_ix + 1]
  t_p11 <- src_pixels[(iy+1) * read_xsize + (t_ix+1) + 1]
  test_val <- (t_p00 * test_rx + t_p10 * (1 - test_rx)) * ry +
    (t_p01 * test_rx + t_p11 * (1 - test_rx)) * (1 - ry)
  test_rounded <- rust_round(test_val)
  if (abs(test_rounded - gdal_px[idx]) < 5) {
    cat(sprintf("  delta=%.1f  buf_x=%.4f  ix=%d  val=%.2f  rounded=%d  ***MATCH***\n",
                delta, test_bx, t_ix, test_val, test_rounded))
  }
}

ds$close()
cat("\nDone.\n")
