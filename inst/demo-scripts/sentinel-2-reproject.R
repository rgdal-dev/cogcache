
rr <- function(x) {
  library(terra)
  rr <- terra::rast(sprintf("/vsicurl/%s", x))
  x <- sprintf("vrt:///vsicurl/%s?ovr=3", x)
  r <- terra::rast(x)
  ## set index extent of full resolution to this overview
  terra::set.ext(r, c(0, ncol(rr), 0, nrow(rr)))
  r
}

vd <- function(x) {
  d <- vapour::gdal_raster_data(sprintf("/vsicurl/%s", x), target_dim = c(1024, 0))
  #attr(d, "extent") <- c(0, 2 * 10980, 0, 2* 10980)
  d
}
dst_crs = "EPSG:5070"
## resolution is 10m
#dst_ex <- c(380000,  420000, 4928000, 4984000)
dst_ex <- c(380000,  390000, 4928000, 4938000)
dst_gt <- vaster::extent_dim_to_gt(dst_ex, diff(dst_ex)[c(1, 3)]/10)
js <- jsonlite::fromJSON(sds::stacit(reproj::reproj_extent(dst_ex, "EPSG:4326", source = dst_crs),
                                     date = c("2023-06-01", "2023-08-31")))

href <- js$features$assets$red$href[255]
dd <- vd(href)
#ximage::ximage(dd)
## info for one of the tifs
info <- vapour::vapour_raster_info(sprintf("/vsicurl/%s", href))
src_gt <- info$geotransform
src_crs <- info$projection


src_dim <- as.integer(info$dimXY)
dst_ncol <- as.integer(diff(dst_ex)[1] / 10)  # 4000
dst_nrow <- as.integer(diff(dst_ex)[3] / 10)  # 5600

chunks <- cogcache:::rust_collect_chunk_list(
  src_crs, src_gt, src_dim,
  dst_crs, dst_gt,
  c(0L, 0L), c(dst_ncol, dst_nrow),
  1L, 0.5, 8L
)
length(chunks)
str(chunks[[1]])



ds <- new(gdalraster::GDALRaster, sprintf("/vsicurl/%s", href), TRUE)

chunk <- chunks[[1]]

src_buf <- ds$read(
  band = 1L,
  xoff = chunk$src_xoff,
  yoff = chunk$src_yoff,
  xsize = chunk$src_xsize,
  ysize = chunk$src_ysize,
  out_xsize = chunk$src_xsize,
  out_ysize = chunk$src_ysize
)

nodata <- as.integer(ds$getNoDataValue(1L))
if (is.na(nodata)) nodata <- 0L

ds$close()

# warp
result_vec <- cogcache:::rust_warp_resample(
  src_crs, src_gt,
  dst_crs, dst_gt,
  c(chunk$dst_xsize, chunk$dst_ysize),
  as.integer(src_buf),
  chunk$src_xsize, chunk$src_ysize,
  chunk$src_xoff,  chunk$src_yoff,
  nodata, 0.125, "bilinear"
)

result <- matrix(result_vec, chunk$dst_ysize, chunk$dst_xsize, byrow = TRUE)
result[result <=0L ] <- NA
ximage::ximage(result, extent = dst_ex)

par(mfrow = c(1, 2))
ximage::ximage(result, extent = dst_ex, main = "Sentinel-2 red EPSG:5070", asp = 1)
ximage::ximage(dd, asp = 1)
vaster::plot_extent(reproj::reproj_extent(dst_ex, info$projection, source = dst_crs), add = T)
