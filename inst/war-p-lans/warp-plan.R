## warp_plan.R
## Analyse the distortion field between a destination and source raster grid
## to map destination tiles to their required source tiles, with annotations
## for subdivision density and discontinuity detection.
##
## Requires: PROJ (>= 0.5.0)
## No other dependencies beyond base R.

#' Plan a raster warp: map destination tiles to source tiles
#'
#' Given a destination grid and source grid (each defined by extent, dimension,
#' and CRS), compute which source tiles are needed for each destination tile,
#' annotated with distortion metrics from the composite projection.
#'
#' The composite Jacobian of the dest->source mapping is computed analytically
#' via PROJ::proj_factors() evaluated at sample points within each destination
#' tile. Three cases are handled:
#'
#' - Source is geographic: composite = J_dest_inv (only dest factors needed)
#' - Dest is geographic: composite = J_source (only source factors needed)
#' - Both projected: composite = J_source %*% J_dest_inv (both needed)
#'
#' Scale factors are converted from CRS-unit ratios to pixel-to-pixel ratios
#' using the grid resolutions, so scale_a = 2.0 means one dest pixel spans
#' 2 source pixels in the direction of maximum stretch.
#'
#' A grid is defined by extent (c(xmin, xmax, ymin, ymax)), dimension
#' (c(ncol, nrow)), and CRS string. Tile layout is derived from tile_dim
#' over the same extent.
#'
#' @param dest_extent numeric length 4: c(xmin, xmax, ymin, ymax) in dest CRS
#' @param dest_dim integer length 2: c(ncol, nrow) pixel dimensions
#' @param dest_crs CRS string for the destination grid
#' @param source_extent numeric length 4: c(xmin, xmax, ymin, ymax) in source CRS
#' @param source_dim integer length 2: c(ncol, nrow) pixel dimensions
#' @param source_crs CRS string for the source grid
#' @param dest_tile_dim integer length 2: tile size in pixels c(ncol, nrow)
#' @param source_tile_dim integer length 2: tile size in pixels
#' @param sample_n integer; sample points per tile edge (default 5 -> 25 per tile)
#' @return A warp_plan data.frame with one row per dest_tile x source_tile pair
warp_plan <- function(dest_extent, dest_dim, dest_crs,
                      source_extent, source_dim, source_crs,
                      dest_tile_dim = c(256L, 256L),
                      source_tile_dim = c(256L, 256L),
                      sample_n = 5L) {

  ## ---- Validate ----
  stopifnot(length(dest_extent) == 4L, length(source_extent) == 4L)
  dest_dim <- as.integer(dest_dim)
  source_dim <- as.integer(source_dim)
  dest_tile_dim <- as.integer(dest_tile_dim)
  source_tile_dim <- as.integer(source_tile_dim)
  sample_n <- as.integer(sample_n)
  stopifnot(length(dest_dim) == 2L, length(source_dim) == 2L)
  stopifnot(all(dest_dim > 0L), all(source_dim > 0L))

  ## Detect geographic CRSs
  dest_is_geo   <- .is_geographic(dest_crs)
  source_is_geo <- .is_geographic(source_crs)
  if (dest_is_geo && source_is_geo) {
    stop("At least one of dest_crs / source_crs must be projected")
  }

  ## ---- Grid geometry ----
  dest_res <- c((dest_extent[2L] - dest_extent[1L]) / dest_dim[1L],
                (dest_extent[4L] - dest_extent[3L]) / dest_dim[2L])
  source_res <- c((source_extent[2L] - source_extent[1L]) / source_dim[1L],
                  (source_extent[4L] - source_extent[3L]) / source_dim[2L])

  dest_ntiles <- c(ceiling(dest_dim[1L] / dest_tile_dim[1L]),
                   ceiling(dest_dim[2L] / dest_tile_dim[2L]))
  source_ntiles <- c(ceiling(source_dim[1L] / source_tile_dim[1L]),
                     ceiling(source_dim[2L] / source_tile_dim[2L]))

  source_tile_width  <- source_tile_dim[1L] * source_res[1L]
  source_tile_height <- source_tile_dim[2L] * source_res[2L]

  ## Pixel-ratio scale factor: converts CRS-unit Jacobian to pixel-to-pixel
  ## For the composite J mapping (dest CRS units -> source CRS units),
  ## multiply by (dest pixel size / source pixel size) to get pixel ratios.
  ## We compute a mean pixel scale for x and y as a single representative value.
  pixel_ratio_x <- abs(dest_res[1L] / source_res[1L])
  pixel_ratio_y <- abs(dest_res[2L] / source_res[2L])
  pixel_ratio <- sqrt(pixel_ratio_x * pixel_ratio_y)  # geometric mean

  ## ---- Sample points for all dest tiles ----
  s <- seq(0.5 / sample_n, 1 - 0.5 / sample_n, length.out = sample_n)
  s_grid <- expand.grid(sx = s, sy = s)
  n_samples <- nrow(s_grid)

  tile_ids <- expand.grid(dtc = seq_len(dest_ntiles[1L]),
                          dtr = seq_len(dest_ntiles[2L]))
  n_tiles <- nrow(tile_ids)

  ## Dest tile bboxes
  tile_xmin <- dest_extent[1L] + (tile_ids$dtc - 1L) * dest_tile_dim[1L] * dest_res[1L]
  tile_ymin <- dest_extent[3L] + (tile_ids$dtr - 1L) * dest_tile_dim[2L] * dest_res[2L]
  tile_xmax <- pmin(tile_xmin + dest_tile_dim[1L] * dest_res[1L], dest_extent[2L])
  tile_ymax <- pmin(tile_ymin + dest_tile_dim[2L] * dest_res[2L], dest_extent[4L])

  ## Expand: tile x sample
  tile_idx <- rep(seq_len(n_tiles), each = n_samples)
  sx <- rep(s_grid$sx, times = n_tiles)
  sy <- rep(s_grid$sy, times = n_tiles)

  pts_x <- tile_xmin[tile_idx] + sx * (tile_xmax[tile_idx] - tile_xmin[tile_idx])
  pts_y <- tile_ymin[tile_idx] + sy * (tile_ymax[tile_idx] - tile_ymin[tile_idx])
  pts_dest <- cbind(pts_x, pts_y)

  ## ---- Transform: dest -> lonlat -> source ----
  pts_ll <- PROJ::proj_trans(pts_dest, source = dest_crs, target = "EPSG:4326")
  pts_source <- PROJ::proj_trans(pts_ll, source = "EPSG:4326", target = source_crs)

  ok <- is.finite(pts_ll[, 1L]) & is.finite(pts_ll[, 2L]) &
    is.finite(pts_source[, 1L]) & is.finite(pts_source[, 2L])

  ## ---- Composite Jacobian via proj_factors ----
  np <- length(ok)
  scale_a <- rep(NA_real_, np)
  scale_b <- rep(NA_real_, np)
  comp_area <- rep(NA_real_, np)
  rotation <- rep(NA_real_, np)

  if (any(ok)) {
    ll_ok <- pts_ll[ok, , drop = FALSE]

    composite <- .composite_jacobian(ll_ok, dest_crs, source_crs,
                                     dest_is_geo, source_is_geo)

    if (!is.null(composite)) {
      idx_full <- which(ok)
      idx_ok2 <- idx_full[composite$ok2]

      ## Apply pixel ratio scaling
      scale_a[idx_ok2] <- composite$sa * pixel_ratio
      scale_b[idx_ok2] <- composite$sb * pixel_ratio
      comp_area[idx_ok2] <- composite$area * pixel_ratio^2
      rotation[idx_ok2] <- composite$rot
    }
  }

  ## ---- Map sample points to source tiles ----
  stc <- floor((pts_source[, 1L] - source_extent[1L]) / source_tile_width) + 1L
  str <- floor((pts_source[, 2L] - source_extent[3L]) / source_tile_height) + 1L

  outside <- !ok |
    stc < 1L | stc > source_ntiles[1L] |
    str < 1L | str > source_ntiles[2L]
  stc[outside] <- NA_integer_
  str[outside] <- NA_integer_

  ## ---- Discontinuity detection per dest tile ----
  source_width  <- source_extent[2L] - source_extent[1L]
  source_height <- source_extent[4L] - source_extent[3L]

  tile_discon <- .detect_discon(tile_idx, ok, pts_source,
                                source_width, source_height,
                                scale_a, n_tiles)

  ## ---- Subdivision level per dest tile ----
  tile_subdiv <- .compute_subdiv(tile_idx, comp_area, n_samples,
                                 sample_n, n_tiles)

  ## ---- Aggregate per dest_tile x source_tile pair ----
  pt <- data.frame(tile_idx = tile_idx, stc = stc, str = str,
                   scale_a = scale_a, scale_b = scale_b,
                   comp_area = comp_area, rotation = rotation)
  pt_valid <- pt[!is.na(pt$stc), ]

  if (nrow(pt_valid) == 0L) {
    warning("No valid dest->source mappings found")
    return(.empty_plan(dest_extent, dest_dim, dest_crs, dest_tile_dim,
                       source_extent, source_dim, source_crs, source_tile_dim,
                       dest_ntiles, source_ntiles))
  }

  ## Unique pairs
  pair_id <- paste(pt_valid$tile_idx, pt_valid$stc, pt_valid$str, sep = ":")
  pairs <- split(pt_valid, pair_id, drop = TRUE)

  rows <- lapply(pairs, function(p) {
    ti <- p$tile_idx[1L]
    sb_safe <- pmax(p$scale_b, 1e-12)
    data.frame(
      dest_tile_col   = tile_ids$dtc[ti],
      dest_tile_row   = tile_ids$dtr[ti],
      source_tile_col = p$stc[1L],
      source_tile_row = p$str[1L],
      n_samples       = nrow(p),
      scale_area_min  = min(p$comp_area, na.rm = TRUE),
      scale_area_max  = max(p$comp_area, na.rm = TRUE),
      scale_a_max     = max(p$scale_a, na.rm = TRUE),
      anisotropy_max  = max(p$scale_a / sb_safe, na.rm = TRUE),
      convergence_sd  = if (sum(is.finite(p$rotation)) > 1L)
        sd(p$rotation, na.rm = TRUE) else NA_real_,
      subdiv          = tile_subdiv[ti],
      discontinuous   = tile_discon[ti],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  ## Clean up Inf/-Inf
  num_cols <- c("scale_area_min", "scale_area_max", "scale_a_max",
                "anisotropy_max", "convergence_sd")
  for (col in num_cols) {
    bad <- !is.finite(out[[col]])
    if (any(bad)) out[[col]][bad] <- NA_real_
  }

  out <- out[order(out$dest_tile_col, out$dest_tile_row,
                   out$source_tile_col, out$source_tile_row), ]

  attr(out, "dest_grid") <- list(
    extent = dest_extent, dim = dest_dim, crs = dest_crs,
    tile_dim = dest_tile_dim, ntiles = dest_ntiles
  )
  attr(out, "source_grid") <- list(
    extent = source_extent, dim = source_dim, crs = source_crs,
    tile_dim = source_tile_dim, ntiles = source_ntiles
  )
  class(out) <- c("warp_plan", "data.frame")
  out
}


## ==== Internal helpers ====

## Detect geographic CRS via PROJ
.is_geographic <- function(crs) {
  ## Try PROJ's own check first, fall back to string matching
  tryCatch({
    PROJ::proj_is_latlong(crs)
  }, error = function(e) {
    grepl("\\+proj=longlat|\\+proj=lonlat|\\+proj=latlong|EPSG:4326|OGC:CRS84",
          crs, ignore.case = TRUE)
  })
}


## Compute composite Jacobian entries for all valid lonlat points
## Returns list(sa, sb, area, rot, ok2) or NULL if nothing worked
.composite_jacobian <- function(ll_ok, dest_crs, source_crs,
                                dest_is_geo, source_is_geo) {

  if (source_is_geo) {
    ## Source is lonlat: dest->source = inverse of dest projection
    ## Composite = J_dest_inv
    fd <- .pf_safe(ll_ok, dest_crs)
    ok2 <- !is.na(fd[, "dx_dlam"])
    if (!any(ok2)) return(NULL)

    det_d <- fd[ok2, "dx_dlam"] * fd[ok2, "dy_dphi"] -
      fd[ok2, "dx_dphi"] * fd[ok2, "dy_dlam"]
    det_d[abs(det_d) < .Machine$double.eps * 100] <- NA_real_

    c_a <-  fd[ok2, "dy_dphi"] / det_d
    c_b <- -fd[ok2, "dx_dphi"] / det_d
    c_c <- -fd[ok2, "dy_dlam"] / det_d
    c_d <-  fd[ok2, "dx_dlam"] / det_d

  } else if (dest_is_geo) {
    ## Dest is lonlat: dest->source = source projection directly
    ## Composite = J_source
    fs <- .pf_safe(ll_ok, source_crs)
    ok2 <- !is.na(fs[, "dx_dlam"])
    if (!any(ok2)) return(NULL)

    c_a <- fs[ok2, "dx_dlam"]
    c_b <- fs[ok2, "dx_dphi"]
    c_c <- fs[ok2, "dy_dlam"]
    c_d <- fs[ok2, "dy_dphi"]

  } else {
    ## Both projected: full composite J_source %*% J_dest_inv
    fd <- .pf_safe(ll_ok, dest_crs)
    fs <- .pf_safe(ll_ok, source_crs)
    ok2 <- !is.na(fd[, "dx_dlam"]) & !is.na(fs[, "dx_dlam"])
    if (!any(ok2)) return(NULL)

    det_d <- fd[ok2, "dx_dlam"] * fd[ok2, "dy_dphi"] -
      fd[ok2, "dx_dphi"] * fd[ok2, "dy_dlam"]
    det_d[abs(det_d) < .Machine$double.eps * 100] <- NA_real_

    inv_a <-  fd[ok2, "dy_dphi"] / det_d
    inv_b <- -fd[ok2, "dx_dphi"] / det_d
    inv_c <- -fd[ok2, "dy_dlam"] / det_d
    inv_d <-  fd[ok2, "dx_dlam"] / det_d

    c_a <- fs[ok2, "dx_dlam"] * inv_a + fs[ok2, "dx_dphi"] * inv_c
    c_b <- fs[ok2, "dx_dlam"] * inv_b + fs[ok2, "dx_dphi"] * inv_d
    c_c <- fs[ok2, "dy_dlam"] * inv_a + fs[ok2, "dy_dphi"] * inv_c
    c_d <- fs[ok2, "dy_dlam"] * inv_b + fs[ok2, "dy_dphi"] * inv_d
  }

  ## 2x2 SVD closed form on composite entries
  tr  <- c_a^2 + c_b^2 + c_c^2 + c_d^2
  det <- c_a * c_d - c_b * c_c
  disc <- sqrt(pmax(tr^2 - 4 * det^2, 0))

  sa <- sqrt(pmax((tr + disc) / 2, 0))
  sb <- sqrt(pmax((tr - disc) / 2, 0))
  area <- abs(det)
  rot <- atan2(c_c - c_b, c_a + c_d) * 180 / pi

  list(sa = sa, sb = sb, area = area, rot = rot, ok2 = ok2)
}


## Safe wrapper around PROJ::proj_factors
.pf_safe <- function(ll, crs) {
  tryCatch(
    PROJ::proj_factors(ll, crs),
    error = function(e) {
      nms <- c("meridional_scale", "parallel_scale", "areal_scale",
               "angular_distortion", "meridian_parallel_angle",
               "meridian_convergence", "tissot_semimajor", "tissot_semiminor",
               "dx_dlam", "dx_dphi", "dy_dlam", "dy_dphi")
      matrix(NA_real_, nrow = nrow(ll), ncol = length(nms),
             dimnames = list(NULL, nms))
    }
  )
}


## Discontinuity detection per dest tile
.detect_discon <- function(tile_idx, ok, pts_source,
                           source_width, source_height,
                           scale_a, n_tiles) {
  discon <- logical(n_tiles)
  for (ti in seq_len(n_tiles)) {
    sel <- tile_idx == ti
    ok_i <- ok[sel]

    ## High NA rate
    if (mean(!ok_i) > 0.3) {
      discon[ti] <- TRUE
      next
    }

    src <- pts_source[sel, , drop = FALSE]
    src_ok <- src[ok_i, , drop = FALSE]
    if (nrow(src_ok) < 2L) next

    ## Coordinate range spans more than half the source extent -> wraparound
    x_range <- diff(range(src_ok[, 1L]))
    y_range <- diff(range(src_ok[, 2L]))
    if (x_range > source_width * 0.5 || y_range > source_height * 0.5) {
      discon[ti] <- TRUE
      next
    }

    ## Scale factor blow-up (in pixel units now)
    sa_i <- scale_a[sel][ok_i]
    if (any(is.finite(sa_i) & sa_i > 100, na.rm = TRUE)) {
      discon[ti] <- TRUE
    }
  }
  discon
}


## Subdivision level per dest tile
.compute_subdiv <- function(tile_idx, comp_area, n_samples,
                            sample_n, n_tiles) {
  subdiv <- integer(n_tiles)
  for (ti in seq_len(n_tiles)) {
    sel <- tile_idx == ti
    vals <- comp_area[sel]
    n_valid <- sum(is.finite(vals))

    if (n_valid < 4L) {
      subdiv[ti] <- 4L
      next
    }

    ## Corner indices in expand.grid(sx, sy) order — sx fastest
    corners <- c(1L, sample_n,
                 (sample_n - 1L) * sample_n + 1L,
                 sample_n * sample_n)
    corner_vals <- vals[corners]

    if (any(!is.finite(corner_vals))) {
      subdiv[ti] <- 4L
      next
    }

    ## Bilinear prediction at center vs actual
    center_pred <- mean(corner_vals)
    center_idx <- ceiling(n_samples / 2)
    center_actual <- vals[center_idx]

    if (!is.finite(center_actual) || abs(center_pred) < 1e-12) {
      subdiv[ti] <- 4L
      next
    }

    rel_err <- abs(center_actual - center_pred) / abs(center_pred)
    valid_vals <- vals[is.finite(vals)]
    cv <- sd(valid_vals) / abs(mean(valid_vals))

    subdiv[ti] <- if (rel_err < 0.02 && cv < 0.05) {
      1L
    } else if (rel_err < 0.10 && cv < 0.20) {
      2L
    } else if (rel_err < 0.25 && cv < 0.50) {
      3L
    } else {
      4L
    }
  }
  subdiv
}


## Empty plan skeleton
.empty_plan <- function(dest_extent, dest_dim, dest_crs, dest_tile_dim,
                        source_extent, source_dim, source_crs, source_tile_dim,
                        dest_ntiles, source_ntiles) {
  out <- data.frame(
    dest_tile_col = integer(0), dest_tile_row = integer(0),
    source_tile_col = integer(0), source_tile_row = integer(0),
    n_samples = integer(0),
    scale_area_min = numeric(0), scale_area_max = numeric(0),
    scale_a_max = numeric(0), anisotropy_max = numeric(0),
    convergence_sd = numeric(0),
    subdiv = integer(0), discontinuous = logical(0),
    stringsAsFactors = FALSE
  )
  attr(out, "dest_grid") <- list(
    extent = dest_extent, dim = dest_dim, crs = dest_crs,
    tile_dim = dest_tile_dim, ntiles = dest_ntiles
  )
  attr(out, "source_grid") <- list(
    extent = source_extent, dim = source_dim, crs = source_crs,
    tile_dim = source_tile_dim, ntiles = source_ntiles
  )
  class(out) <- c("warp_plan", "data.frame")
  out
}


## Print method
print.warp_plan <- function(x, ...) {
  dg <- attr(x, "dest_grid")
  sg <- attr(x, "source_grid")
  cat(sprintf("Warp plan: %d dest tiles x %d source tiles -> %d mappings\n",
              prod(dg$ntiles), prod(sg$ntiles), nrow(x)))
  cat(sprintf("  Dest:   %dx%d px, %dx%d tiles, %s\n",
              dg$dim[1L], dg$dim[2L], dg$ntiles[1L], dg$ntiles[2L],
              dg$crs))
  cat(sprintf("  Source: %dx%d px, %dx%d tiles, %s\n",
              sg$dim[1L], sg$dim[2L], sg$ntiles[1L], sg$ntiles[2L],
              sg$crs))
  if (nrow(x) > 0L) {
    n_discon <- sum(x$discontinuous, na.rm = TRUE)
    if (n_discon > 0L) {
      cat(sprintf("  Discontinuous: %d tile pairs\n", n_discon))
    }
    cat(sprintf("  Areal scale: [%.3f, %.3f] (source px per dest px)\n",
                min(x$scale_area_min, na.rm = TRUE),
                max(x$scale_area_max, na.rm = TRUE)))
    cat(sprintf("  Max anisotropy: %.2f\n",
                max(x$anisotropy_max, na.rm = TRUE)))
    subdiv_tab <- table(x$subdiv[!duplicated(paste(x$dest_tile_col, x$dest_tile_row))])
    cat(sprintf("  Subdivision levels: %s\n",
                paste(names(subdiv_tab), subdiv_tab, sep = ":", collapse = "  ")))
  }
  invisible(x)
}
