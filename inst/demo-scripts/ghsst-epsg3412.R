library(rustycogs)
library(cogcache)
library(dplyr)

# -- config ------------------------------------------------------------------
date   <- as.POSIXct("2025-11-26", tz = "UTC")
url    <- sprintf(
  "https://data.source.coop/ausantarctic/ghrsst-mur-v2/%s/%s090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1_analysed_sst.tif",
  format(date, "%Y/%m/%d"), format(date, "%Y%m%d")
)
nodata <- -32768L

# -- source grid (exact GHRSST geotransform) ---------------------------------
src_crs  <- "EPSG:4326"
src_gt   <- c(-179.995, 0.010, 0.0, 89.995, 0.0, -0.010)
src_full <- c(36000L, 17999L)

# -- target grid: EPSG:3412 (Antarctic), 25km, 40x40 centred lon=0 lat=-60 --
dst_crs  <- "EPSG:3412"
dst_res  <- 25000

pt <- reproj::reproj_xy(cbind(0, -60), target = dst_crs, source = "EPSG:4326")
cat(sprintf("lon=0, lat=-60 in EPSG:3412: x=%.0f y=%.0f\n", pt[1], pt[2]))

dst_ncol <- 40L
dst_nrow <- 40L
dst_gt   <- c(pt[1] - 20*dst_res, dst_res, 0, pt[2] + 20*dst_res, 0, -dst_res)
dst_off  <- c(0L, 0L)
dst_size <- c(dst_ncol, dst_nrow)

# -- byte-ref table ----------------------------------------------------------
message("fetching byte refs...")
refs <- rustycogs::tiff_refs(url, "", TRUE, concurrency = 180)
message(sprintf("%d tile refs", nrow(refs)))

# -- plan --------------------------------------------------------------------
message("planning...")
chunks <- cogcache:::rust_collect_chunk_list(
  src_crs, src_gt, src_full,
  dst_crs, dst_gt,
  dst_off, dst_size,
  1L, 0.5, 8L
)
message(sprintf("%d chunk(s)", length(chunks)))

# -- helpers -----------------------------------------------------------------
lookup_src_windows <- function(refs, chunk) {
  refs |>
    filter(
      ifd == 0L,
      tile_col >= chunk$src_xoff %/% 512L,
      tile_col <= (chunk$src_xoff + chunk$src_xsize - 1L) %/% 512L,
      tile_row >= chunk$src_yoff %/% 512L,
      tile_row <= (chunk$src_yoff + chunk$src_ysize - 1L) %/% 512L
    ) |>
    mutate(
      px_xoff  = tile_col * 512L,
      px_yoff  = tile_row * 512L,
      px_xsize = pmin(512L, 36000L - px_xoff),
      px_ysize = pmin(512L, 17999L - px_yoff)
    )
}

assemble_src_buffer <- function(hits, chunk) {
  buf <- matrix(nodata, nrow = chunk$src_ysize, ncol = chunk$src_xsize)
  for (i in seq_len(nrow(hits))) {
    tile <- tiff_tile(hits$path[i], ifd_index = 0L,
                      col = hits$tile_col[i], row = hits$tile_row[i],
                      region = "", anon = TRUE)
    m  <- matrix(tile$data,
                 nrow = hits$px_ysize[i],
                 ncol = hits$px_xsize[i],
                 byrow = TRUE)
    bx <- hits$px_xoff[i] - chunk$src_xoff
    by <- hits$px_yoff[i] - chunk$src_yoff

    x1 <- max(1L, bx + 1L);  x2 <- min(chunk$src_xsize, bx + hits$px_xsize[i])
    y1 <- max(1L, by + 1L);  y2 <- min(chunk$src_ysize, by + hits$px_ysize[i])
    if (x1 > x2 || y1 > y2) next

    mx1 <- x1 - bx;  mx2 <- x2 - bx
    my1 <- y1 - by;  my2 <- y2 - by

    buf[y1:y2, x1:x2] <- m[my1:my2, mx1:mx2]
  }
  as.integer(as.vector(t(buf)))
}

# -- fetch + warp ------------------------------------------------------------
dst_m <- matrix(nodata, dst_nrow, dst_ncol)

for (i in seq_along(chunks)) {
  chunk <- chunks[[i]]

  message(sprintf(
    "chunk %d/%d  src[%d,%d %dx%d]  dst[%d,%d %dx%d]  fill=%.3f",
    i, length(chunks),
    chunk$src_xoff, chunk$src_yoff, chunk$src_xsize, chunk$src_ysize,
    chunk$dst_xoff, chunk$dst_yoff, chunk$dst_xsize, chunk$dst_ysize,
    chunk$fill_ratio
  ))

  if (chunk$fill_ratio < 0.1) { message("  skipping"); next }

  hits <- lookup_src_windows(refs, chunk)
  message(sprintf("  n_tiles=%d", nrow(hits)))
  if (nrow(hits) == 0) next

  src_buf <- assemble_src_buffer(hits, chunk)
  message(sprintf("  valid=%.1f%%", mean(src_buf != nodata) * 100))

  dst_chunk <- cogcache:::rust_warp_resample(
    src_crs, src_gt,
    dst_crs, dst_gt,
    c(chunk$dst_xsize, chunk$dst_ysize),
    as.integer(src_buf),
    chunk$src_xsize, chunk$src_ysize,
    chunk$src_xoff,  chunk$src_yoff,
    nodata, 0.125, "bilinear"
  )

  dx <- chunk$dst_xoff - dst_off[1]
  dy <- chunk$dst_yoff - dst_off[2]
  chunk_m <- matrix(dst_chunk, chunk$dst_ysize, chunk$dst_xsize, byrow = TRUE)
  dst_m[(dy+1L):(dy+chunk$dst_ysize), (dx+1L):(dx+chunk$dst_xsize)] <- chunk_m
}

# -- view --------------------------------------------------------------------
result <- dst_m
result[result == nodata] <- NA
cat(sprintf("valid: %d / %d  (%.1f%%)\n",
            sum(!is.na(result)), length(result),
            mean(!is.na(result)) * 100))

ximage::ximage(result, main = "GHRSST SST - EPSG:3412 40x40 lon=0 lat=-60 (2025-11-26)")
