## Plot the source data neighbourhood around the problem pixel
## Show contours, our landing point, and where GDAL's -4030 would come from

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

## Read a 20x20 neighbourhood in full-image coords
## Our problem pixel: full-image col=84874, row=25555 (the bilinear ix,iy)
cx <- 84874L
cy <- 25555L
hw <- 10L  # half-width
region <- ds$read(band = 1,
                  xoff = cx - hw, yoff = cy - hw,
                  xsize = 2L * hw, ysize = 2L * hw,
                  out_xsize = 2L * hw, out_ysize = 2L * hw)
ds$close()

mat <- matrix(region, nrow = 2 * hw, ncol = 2 * hw, byrow = TRUE)

## Pixel centres in full-image coordinates
cols <- seq(cx - hw, cx + hw - 1)  # full-image column indices
rows <- seq(cy - hw, cy + hw - 1)  # full-image row indices

## Our exact landing point (from rust_gen_img_proj_transform)
our_src_x <- 84874.978478
our_src_y <- 25556.325546

## The bilinear 2x2 neighbourhood corners (ix=394+read_xoff=84874, iy=979+read_yoff=24576+979=25555)
## Actually: ix = floor(buf_x - 0.5) where buf_x = src_x - read_xoff
## buf_x = 84874.978478 - 84480 = 394.978478
## ix = floor(394.978478 - 0.5) = floor(394.478478) = 394
## full-image: 84480 + 394 = 84874. Correct.
## iy = floor(980.325546 - 0.5) = floor(979.825546) = 979
## full-image: 24576 + 979 = 25555. Correct.

png("bilinear_neighbourhood.png", width = 900, height = 800)
par(mar = c(5, 5, 4, 2))

## Image plot with pixel centres
image(cols + 0.5, rows + 0.5, t(mat),
      col = hcl.colors(50, "viridis"),
      xlab = "Full-image column", ylab = "Full-image row",
      main = "Source bathymetry around bilinear kernel at dst [1,157]",
      useRaster = TRUE, asp = 1)

## Contours
contour(cols + 0.5, rows + 0.5, t(mat),
        levels = seq(-6500, 0, by = 250),
        add = TRUE, col = "white", labcex = 0.7)

## Our landing point
points(our_src_x, our_src_y, pch = 4, cex = 3, col = "red", lwd = 3)
text(our_src_x + 0.3, our_src_y + 0.5,
     sprintf("Our point\n(%.2f, %.2f)\nbilinear=-5699", our_src_x, our_src_y),
     col = "red", adj = 0, cex = 0.9, font = 2)

## The 2x2 bilinear kernel box
rect(84874, 25555, 84876, 25557, border = "red", lwd = 2, lty = 2)

## Mark the 4 bilinear pixels with their values
text(84874.5, 25555.5, "-6306", col = "yellow", cex = 0.8, font = 2)
text(84875.5, 25555.5, "-5982", col = "yellow", cex = 0.8, font = 2)
text(84874.5, 25556.5, "-5726", col = "yellow", cex = 0.8, font = 2)
text(84875.5, 25556.5, "-5469", col = "yellow", cex = 0.8, font = 2)

## Where would GDAL's answer of -4030 come from?
## -4030 appears at full-image (84877, 25556) in our 6x6 dump
## and also nearby. Let's find the -4030 contour
## Actually, let's mark where -4028 and -4030 appear
for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    if (abs(mat[i, j] - (-4030)) < 50) {
      ## Mark candidate pixels
      points(cols[j] + 0.5, rows[i] + 0.5, pch = 1, cex = 1.5, col = "cyan", lwd = 2)
    }
  }
}

## Add the -4030 contour specifically
contour(cols + 0.5, rows + 0.5, t(mat),
        levels = -4030,
        add = TRUE, col = "cyan", lwd = 2, labcex = 0.9)

legend("topright",
       legend = c("Our bilinear point", "2x2 kernel box",
                  "-4030 contour (GDAL answer)", "Pixels near -4030"),
       pch = c(4, NA, NA, 1),
       lty = c(NA, 2, 1, NA),
       col = c("red", "red", "cyan", "cyan"),
       lwd = c(3, 2, 2, 2),
       bg = "grey20", text.col = "white", cex = 0.9)

dev.off()
cat("Plot saved to bilinear_neighbourhood.png\n")

## Also: compute a bilinear surface over the region to see what
## coordinates would give -4030
cat("\nBilinear sweep: where in the region would bilinear = -4030?\n")
cat("(Sweeping 0.1-pixel steps over a 6x6 area)\n")
found <- FALSE
for (test_x in seq(84872, 84878, by = 0.05)) {
  for (test_y in seq(25553, 25559, by = 0.05)) {
    tix <- floor(test_x - 0.5)
    tiy <- floor(test_y - 0.5)
    # check in range of our matrix
    ci <- tix - (cx - hw)  # 0-based col in mat
    ri <- tiy - (cy - hw)  # 0-based row in mat
    if (ci >= 0 && ci + 1 < 2*hw && ri >= 0 && ri + 1 < 2*hw) {
      trx <- 1.5 - (test_x - tix)
      try_ <- 1.5 - (test_y - tiy)
      p00 <- mat[ri + 1, ci + 1]      # R is 1-based
      p10 <- mat[ri + 1, ci + 2]
      p01 <- mat[ri + 2, ci + 1]
      p11 <- mat[ri + 2, ci + 2]
      val <- (p00 * trx + p10 * (1 - trx)) * try_ +
        (p01 * trx + p11 * (1 - trx)) * (1 - try_)
      rounded <- ifelse(val >= 0, as.integer(val + 0.5), as.integer(val - 0.5))
      if (abs(rounded - (-4030)) <= 1) {
        cat(sprintf("  src=(%.2f, %.2f) -> bilinear=%.1f -> rounded=%d\n",
                    test_x, test_y, val, rounded))
        found <- TRUE
      }
    }
  }
}
if (!found) cat("  No bilinear location in this region produces -4030!\n")

cat("\nDone.\n")
