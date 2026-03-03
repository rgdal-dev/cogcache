## Diagnostic: is GDAL reading from the same source data as us?
##
## Everything matches: coordinates, geotransform, kernel formula, PROJ.
## Yet results differ by 1669 at a steep gradient.
##
## Remaining hypothesis: GDAL reads different source pixels.
## - Different overview level?
## - Different source window → different tile boundaries?
## - COG internal tiling causing different reads?

library(gdalraster)
library(vaster)

cog_url <- "/vsicurl/https://data.source.coop/alexgleith/gebco-2024/GEBCO_2024.tif"

## Check overview structure
ds <- new(GDALRaster, cog_url)
cat(sprintf("Full resolution: %d x %d\n", ds$getRasterXSize(), ds$getRasterYSize()))
cat(sprintf("Block size: %s\n", paste(ds$getBlockSize(1), collapse="x")))
n_ovr <- ds$getOverviewCount(1)
cat(sprintf("Overview count: %d\n", n_ovr))
## gdalraster doesn't expose overview sizes directly, let's use gdalinfo
ds$close()

# ## Use gdalinfo to check overviews
# info <- gdalraster::gdalinfo(cog_url)
# cat("\n--- gdalinfo overviews ---\n")
# ovr_lines <- grep("Overview", info, value = TRUE, ignore.case = TRUE)
# cat(paste(ovr_lines, collapse = "\n"), "\n")

## Now the key test: extract source pixel values that GDAL actually uses.
## We know our pixel [1,157] maps to full-image src (84874.98, 25556.33).
## Read a small window around that location and compare.
ds <- new(GDALRaster, cog_url)
src_gt <- ds$getGeoTransform()

## Read 10x10 around the bilinear neighbourhood at the exact source location
## ix=394 iy=979 in our buffer, which is full-image col=84874, row=25555
## (read_xoff=84480, so 84480+394=84874)
read_xoff <- 84480L
read_yoff <- 24576L
full_col <- read_xoff + 394L  # = 84874
full_row <- read_yoff + 979L  # = 25555

cat(sprintf("\nFull-image location of bilinear kernel: col=%d row=%d\n", full_col, full_row))

## Read directly from full-res at this location
direct <- ds$read(band = 1, xoff = full_col - 2L, yoff = full_row - 2L,
                  xsize = 6L, ysize = 6L,
                  out_xsize = 6L, out_ysize = 6L)
cat("\nDirect read from full-res (6x6 around kernel):\n")
m <- matrix(direct, nrow = 6, ncol = 6, byrow = TRUE)
for (i in 1:6) {
  cat(sprintf("  row %d: %s\n", full_row - 2 + i - 1,
              paste(sprintf("%6d", m[i,]), collapse = " ")))
}

## Now: what does GDAL warp actually do with this pixel?
## Use nearest-neighbour to see WHERE GDAL thinks [1,157] maps
## (nearest picks a single source pixel, revealing the coordinate)
fiji_lcc <- paste(
  "+proj=lcc +lat_0=-18 +lon_0=178",
  "+lat_1=-10 +lat_2=-25",
  "+datum=WGS84 +units=m +no_defs"
)
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

## GDAL nearest-neighbour warp
gdal_nn <- tempfile(fileext = ".tif")
warp(cog_url, gdal_nn, out_crs,
     cl_arg = c("-ts", tile_size, tile_size,
                "-te", tile_ext[c(1, 3, 2, 4)],
                "-r", "near", "-ovr", "NONE", "-et", "0"))
gds_nn <- new(GDALRaster, gdal_nn)
gdal_nn_px <- gds_nn$read(band = 1, xoff = 0, yoff = 0,
                          xsize = tile_size, ysize = tile_size,
                          out_xsize = tile_size, out_ysize = tile_size)
gds_nn$close()

## Our nearest-neighbour
tile_gt <- extent_dim_to_gt(tile_ext, c(tile_size, tile_size))
src_dim <- c(ds$getRasterXSize(), ds$getRasterYSize())
sw <- rust_compute_source_window(
  src_crs = ds$getProjectionRef(), src_gt = src_gt, src_dim = as.integer(src_dim),
  dst_crs = out_crs, dst_gt = tile_gt,
  dst_off = c(0L, 0L), dst_size = c(tile_size, tile_size),
  resample_padding = 0L
)
src_block <- ds$getBlockSize(band = 1)
src_tw <- src_block[1]; src_th <- src_block[2]
r_xoff  <- (sw$xoff %/% src_tw) * src_tw
r_yoff  <- (sw$yoff %/% src_th) * src_th
r_xend  <- min(ceiling((sw$xoff + sw$xsize) / src_tw) * src_tw, src_dim[1])
r_yend  <- min(ceiling((sw$yoff + sw$ysize) / src_th) * src_th, src_dim[2])
r_xsize <- as.integer(r_xend - r_xoff)
r_ysize <- as.integer(r_yend - r_yoff)
src_pixels <- ds$read(band = 1, xoff = r_xoff, yoff = r_yoff,
                      xsize = r_xsize, ysize = r_ysize,
                      out_xsize = r_xsize, out_ysize = r_ysize)

our_nn <- rust_warp_resample(
  src_crs = ds$getProjectionRef(), src_gt = src_gt, dst_crs = out_crs, dst_gt = tile_gt,
  dst_dim = c(tile_size, tile_size), src_pixels = src_pixels,
  src_ncol = r_xsize, src_nrow = r_ysize,
  src_col_off = as.integer(r_xoff), src_row_off = as.integer(r_yoff),
  nodata = -32768L, max_error = 0.0, resample = "near"
)

idx <- 157 * tile_size + 1 + 1
cat(sprintf("\nNearest-neighbour at [1,157]:\n"))
cat(sprintf("  Ours: %d   GDAL: %d\n", our_nn[idx], gdal_nn_px[idx]))

## Check overall nearest match
nn_match <- sum(our_nn == gdal_nn_px)
cat(sprintf("  Overall NN match: %d/%d\n", nn_match, tile_size^2))

## Check the source pixel value that GDAL NN picks
## Our coordinates: src=(84874.98, 25556.33), truncated = (84874, 25556)
## = full-image col 84874, row 25556
direct_at_trunc <- ds$read(band = 1,
                           xoff = 84874L, yoff = 25556L,
                           xsize = 1L, ysize = 1L,
                           out_xsize = 1L, out_ysize = 1L)
cat(sprintf("  Direct read at (84874, 25556): %d\n", direct_at_trunc))

## The GDAL NN value tells us what source pixel GDAL landed on
## Search for that value nearby
gdal_nn_val <- gdal_nn_px[idx]
cat(sprintf("\n  GDAL NN picked value: %d\n", gdal_nn_val))
cat("  Searching for that value in 10x10 around our source location:\n")
search <- ds$read(band = 1, xoff = full_col - 5L, yoff = full_row - 5L,
                  xsize = 10L, ysize = 10L,
                  out_xsize = 10L, out_ysize = 10L)
search_m <- matrix(search, nrow = 10, ncol = 10, byrow = TRUE)
for (i in 1:10) {
  for (j in 1:10) {
    if (search_m[i,j] == gdal_nn_val) {
      cat(sprintf("    Found %d at full-image (%d, %d)\n",
                  gdal_nn_val, full_col - 5 + j - 1, full_row - 5 + i - 1))
    }
  }
}

ds$close()
cat("\nDone.\n")
