#' Plan efficient source reads for a warp operation
#'
#' Given a destination tile and source/destination CRS and geotransform,
#' compute the optimal set of source reads. For tiles that straddle the
#' antimeridian, this produces a split-read plan that avoids reading the
#' entire source latitude band.
#'
#' @param src_crs Character CRS string for source
#' @param src_gt Numeric length-6 source geotransform
#' @param src_dim Integer c(ncol, nrow) of source raster
#' @param dst_crs Character CRS string for destination
#' @param dst_gt Numeric length-6 destination geotransform
#' @param dst_off Integer c(col_off, row_off) destination window offset
#' @param dst_size Integer c(ncol, nrow) destination window size
#' @param resample_padding Integer kernel padding
#' @param min_fill_ratio Numeric threshold for subdivision (default 0.5)
#' @param min_dst_size Integer minimum chunk dimension (default 8)
#' @return A list of chunk plans
#' @export
plan_warp_reads <- function(src_crs, src_gt, src_dim, dst_crs, dst_gt,
                            dst_off, dst_size, resample_padding = 0L,
                            min_fill_ratio = 0.5,
                            min_dst_size = 8L) {

  chunks <- rust_collect_chunk_list(
    src_crs = src_crs, src_gt = src_gt, src_dim = as.integer(src_dim),
    dst_crs = dst_crs, dst_gt = dst_gt,
    dst_off = as.integer(dst_off), dst_size = as.integer(dst_size),
    resample_padding = as.integer(resample_padding),
    min_fill_ratio = min_fill_ratio,
    min_dst_size = as.integer(min_dst_size)
  )

  lapply(chunks, function(ch) {
    sw <- list(
      xoff = ch$src_xoff, yoff = ch$src_yoff,
      xsize = ch$src_xsize, ysize = ch$src_ysize
    )

    if (sw$xsize > src_dim[1] * 0.5 && ch$fill_ratio < min_fill_ratio) {
      strips <- detect_split_strips(
        src_crs = src_crs, src_gt = src_gt, src_dim = src_dim,
        dst_crs = dst_crs, dst_gt = dst_gt,
        dst_off = c(ch$dst_xoff, ch$dst_yoff),
        dst_size = c(ch$dst_xsize, ch$dst_ysize)
      )
      if (!is.null(strips)) {
        return(list(
          dst_off = c(ch$dst_xoff, ch$dst_yoff),
          dst_size = c(ch$dst_xsize, ch$dst_ysize),
          src_reads = strips,
          fill_ratio = ch$fill_ratio,
          is_split = TRUE
        ))
      }
    }

    list(
      dst_off = c(ch$dst_xoff, ch$dst_yoff),
      dst_size = c(ch$dst_xsize, ch$dst_ysize),
      src_reads = list(sw),
      fill_ratio = ch$fill_ratio,
      is_split = FALSE
    )
  })
}


#' Detect bimodal source X distribution and compute two strip windows
#' @keywords internal
detect_split_strips <- function(src_crs, src_gt, src_dim,
                                dst_crs, dst_gt,
                                dst_off, dst_size) {

  # Adjust geotransform origin for the chunk's position within the tile.
  # rust_gen_img_proj_transform transforms pixels (0..ncol, 0..nrow)
  # relative to the GT origin, so we shift the origin to the chunk corner.
  chunk_gt <- dst_gt
  chunk_gt[1] <- dst_gt[1] + dst_off[1] * dst_gt[2] + dst_off[2] * dst_gt[3]
  chunk_gt[4] <- dst_gt[4] + dst_off[1] * dst_gt[5] + dst_off[2] * dst_gt[6]

  coords <- rust_gen_img_proj_transform(
    src_crs = src_crs, src_gt = src_gt,
    dst_crs = dst_crs, dst_gt = chunk_gt,
    dst_dim = as.integer(dst_size)
  )

  src_x <- coords$src_x
  src_y <- coords$src_y
  valid <- !is.nan(src_x)
  if (sum(valid) < 10) return(NULL)

  sx <- src_x[valid]
  sy <- src_y[valid]

  sx_sorted <- sort(sx)
  gaps <- diff(sx_sorted)
  max_gap_idx <- which.max(gaps)
  max_gap <- gaps[max_gap_idx]

  if (max_gap < src_dim[1] * 0.5) return(NULL)

  gap_left <- sx_sorted[max_gap_idx]
  gap_right <- sx_sorted[max_gap_idx + 1]

  s1_mask <- sx <= gap_left
  s2_mask <- sx >= gap_right
  if (sum(s1_mask) < 2 || sum(s2_mask) < 2) return(NULL)

  y_min <- max(0L, as.integer(floor(min(sy))) - 2L)
  y_max <- min(as.integer(src_dim[2]), as.integer(ceiling(max(sy))) + 2L)

  s1_xmin <- max(0L, as.integer(floor(min(sx[s1_mask]))) - 2L)
  s1_xmax <- min(as.integer(src_dim[1]),
                 as.integer(ceiling(max(sx[s1_mask]))) + 2L)

  s2_xmin <- max(0L, as.integer(floor(min(sx[s2_mask]))) - 2L)
  s2_xmax <- min(as.integer(src_dim[1]),
                 as.integer(ceiling(max(sx[s2_mask]))) + 2L)

  list(
    list(xoff = s1_xmin, yoff = y_min,
         xsize = s1_xmax - s1_xmin, ysize = y_max - y_min),
    list(xoff = s2_xmin, yoff = y_min,
         xsize = s2_xmax - s2_xmin, ysize = y_max - y_min)
  )
}
