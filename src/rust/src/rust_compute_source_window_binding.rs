// --- Add to lib.rs: module declaration ---
// pub mod source_window;

// --- Add to lib.rs: extendr wrapper ---

/// Compute the source pixel window needed for a destination window.
///
/// This is the GDAL ComputeSourceWindow equivalent: given a destination
/// pixel window (offset + size) and CRS/geotransform specs for both source
/// and destination, determine which source pixels are needed.
///
/// Uses the same algorithm as GDAL: edge sampling with grid fallback,
/// antimeridian heuristic, resampling padding.
///
/// @param src_crs Character CRS string
/// @param src_gt Numeric vector length 6 (source geotransform)
/// @param src_dim Integer vector c(ncol, nrow) — full source raster size
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector length 6 (dest geotransform)
/// @param dst_off Integer vector c(col_off, row_off) — dest window offset
/// @param dst_size Integer vector c(ncol, nrow) — dest window size
/// @param resample_padding Integer — extra pixels for resampling kernel
///   (0 = nearest, 1 = bilinear, 2 = cubic, 3 = lanczos)
/// @return List with xoff, yoff, xsize, ysize, fill_ratio, n_failed, n_samples
///   or NULL if computation fails (too few valid transform points)
/// @export
#[extendr]
fn rust_compute_source_window(
    src_crs: &str,
    src_gt: &[f64],
    src_dim: &[i32],
    dst_crs: &str,
    dst_gt: &[f64],
    dst_off: &[i32],
    dst_size: &[i32],
    resample_padding: i32,
) -> extendr_api::Result<Robj> {
    if src_gt.len() != 6 || dst_gt.len() != 6 {
        return Err(Error::Other("geotransform must be length 6".to_string()));
    }
    if src_dim.len() != 2 || dst_off.len() != 2 || dst_size.len() != 2 {
        return Err(Error::Other(
            "src_dim, dst_off, dst_size must be length 2".to_string(),
        ));
    }

    let src_gt_arr: [f64; 6] = src_gt
        .try_into()
        .map_err(|_| Error::Other("src_gt conversion failed".to_string()))?;
    let dst_gt_arr: [f64; 6] = dst_gt
        .try_into()
        .map_err(|_| Error::Other("dst_gt conversion failed".to_string()))?;

    let transformer = transform::GenImgProjTransformer::new(
        src_crs, src_gt_arr, dst_crs, dst_gt_arr,
    )
    .map_err(Error::Other)?;

    let result = source_window::compute_source_window(
        &transformer,
        [dst_off[0], dst_off[1]],
        [dst_size[0], dst_size[1]],
        [src_dim[0], src_dim[1]],
        resample_padding,
    );

    match result {
        Some(sw) => {
            let out = list!(
                xoff = sw.xoff,
                yoff = sw.yoff,
                xsize = sw.xsize,
                ysize = sw.ysize,
                fill_ratio = sw.fill_ratio,
                n_failed = sw.n_failed,
                n_samples = sw.n_samples
            );
            Ok(out.into())
        }
        None => Ok(r!(NULL)),
    }
}

// --- Add to extendr_module! ---
// fn rust_compute_source_window;
