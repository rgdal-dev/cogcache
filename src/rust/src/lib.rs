// cogcache: Rust-native COG warp pipeline
//
// Module structure:
//   transform      — Transformer trait, GenImgProjTransformer (mirrors GDAL gdaltransformer.cpp)
//   approx         — ApproxTransformer (mirrors GDAL GDALApproxTransform)
//   source_window  — ComputeSourceWindow (mirrors GDAL gdalwarpoperation.cpp)
//   warp           — Warp kernel (mirrors GDAL gdalwarpkernel.cpp)
//
// The extendr wrappers below expose these to R.
// The Rust modules are designed to be extractable as standalone crates.

use extendr_api::prelude::*;
use flate2::read::ZlibDecoder;
use std::io::Read;

pub mod transform;
pub mod approx;
pub mod source_window;
pub mod warp;

// ===========================================================================
// Tile decode (existing, unchanged)
// ===========================================================================

fn decode_deflate_pixels(
    compressed: &[u8],
    tile_width: usize,
    tile_height: usize,
    predictor: i32,
) -> std::result::Result<Vec<i32>, String> {
    let n_pixels = tile_width * tile_height;

    let mut decoder = ZlibDecoder::new(compressed);
    let mut decompressed = Vec::with_capacity(n_pixels * 2);
    decoder
        .read_to_end(&mut decompressed)
        .map_err(|e| format!("DEFLATE decompression failed: {}", e))?;

    if decompressed.len() != n_pixels * 2 {
        return Err(format!(
            "Decompressed size {} != expected {} ({}x{}x2)",
            decompressed.len(),
            n_pixels * 2,
            tile_width,
            tile_height
        ));
    }

    let mut pixels: Vec<i32> = decompressed
        .chunks_exact(2)
        .map(|b| u16::from_le_bytes([b[0], b[1]]) as i32)
        .collect();

    if predictor == 2 {
        for row in 0..tile_height {
            let row_start = row * tile_width;
            for col in 1..tile_width {
                let idx = row_start + col;
                pixels[idx] = (pixels[idx] + pixels[idx - 1]) % 65536;
            }
        }
    }

    Ok(pixels)
}

// ===========================================================================
// extendr wrappers — tile decode
// ===========================================================================

/// Decode a DEFLATE-compressed tile with predictor=2 undo.
/// @param raw_bytes Raw vector of compressed tile bytes
/// @param tile_width Integer tile width in pixels
/// @param tile_height Integer tile height in pixels
/// @param predictor Integer predictor type (1 = none, 2 = horizontal diff)
/// @return Integer vector of decoded UInt16 pixel values
/// @export
#[extendr]
fn rust_decode_tile(
    raw_bytes: Raw,
    tile_width: i32,
    tile_height: i32,
    predictor: i32,
) -> extendr_api::Result<Robj> {
    let pixels = decode_deflate_pixels(
        raw_bytes.as_slice(),
        tile_width as usize,
        tile_height as usize,
        predictor,
    )
    .map_err(Error::Other)?;
    Ok(r!(pixels))
}

/// Fetch a tile via HTTP range request and decode it.
/// @param url Character URL
/// @param byte_offset Numeric byte offset
/// @param byte_length Integer byte count
/// @param tile_width Integer tile width
/// @param tile_height Integer tile height
/// @param predictor Integer predictor type
/// @return Integer vector of decoded pixel values
/// @export
#[extendr]
fn rust_fetch_decode_tile(
    url: &str,
    byte_offset: f64,
    byte_length: i32,
    tile_width: i32,
    tile_height: i32,
    predictor: i32,
) -> extendr_api::Result<Robj> {
    let offset = byte_offset as u64;
    let length = byte_length as u64;
    let range_end = offset + length - 1;
    let range_header = format!("bytes={}-{}", offset, range_end);

    let resp = ureq::get(url)
        .set("Range", &range_header)
        .call()
        .map_err(|e| Error::Other(format!("HTTP request failed: {}", e)))?;

    let mut compressed = Vec::with_capacity(byte_length as usize);
    resp.into_reader()
        .read_to_end(&mut compressed)
        .map_err(|e| Error::Other(format!("Failed to read response body: {}", e)))?;

    let pixels = decode_deflate_pixels(
        &compressed,
        tile_width as usize,
        tile_height as usize,
        predictor,
    )
    .map_err(Error::Other)?;
    Ok(r!(pixels))
}

// ===========================================================================
// extendr wrappers — transformer
// ===========================================================================

/// Create a GenImgProjTransformer and transform destination pixel coordinates
/// to source pixel coordinates.
///
/// This is the GDAL GenImgProjTransform equivalent: for each destination pixel
/// centre, compute the corresponding source pixel coordinate via:
///   dst pixel → dst geo → src geo → src pixel
///
/// Returns a list with `src_x` and `src_y` (fractional source pixel coords,
/// 0-based, NaN for failed transforms).
///
/// @param src_crs Character CRS string
/// @param src_gt Numeric vector length 6 (source geotransform)
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector length 6 (dest geotransform)
/// @param dst_dim Integer vector c(ncol, nrow)
/// @return List with `src_x`, `src_y` numeric vectors (fractional pixel coords)
/// @export
#[extendr]
fn rust_gen_img_proj_transform(
    src_crs: &str,
    src_gt: &[f64],
    dst_crs: &str,
    dst_gt: &[f64],
    dst_dim: &[i32],
) -> extendr_api::Result<Robj> {
    if src_gt.len() != 6 || dst_gt.len() != 6 {
        return Err(Error::Other("geotransform must be length 6".to_string()));
    }
    if dst_dim.len() != 2 {
        return Err(Error::Other("dim must be length 2 (ncol, nrow)".to_string()));
    }

    let src_gt_arr: [f64; 6] = src_gt.try_into()
        .map_err(|_| Error::Other("src_gt conversion failed".to_string()))?;
    let dst_gt_arr: [f64; 6] = dst_gt.try_into()
        .map_err(|_| Error::Other("dst_gt conversion failed".to_string()))?;
    let dst_dim_arr = [dst_dim[0] as usize, dst_dim[1] as usize];

    let transformer = transform::GenImgProjTransformer::new(
        src_crs, src_gt_arr, dst_crs, dst_gt_arr,
    )
    .map_err(Error::Other)?;

    let (src_x, src_y) = transform::transform_grid(&transformer, &dst_dim_arr);

    let result = list!(src_x = src_x, src_y = src_y);
    Ok(result.into())
}

/// Warp source pixels to destination grid using GenImgProjTransformer.
///
/// This is the full per-scanline warp kernel: for each destination scanline,
/// transform pixel coordinates to source, then sample source pixels.
/// Nearest neighbour resampling.
///
/// @param src_crs Character CRS string
/// @param src_gt Numeric vector length 6
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector length 6
/// @param dst_dim Integer vector c(ncol, nrow)
/// @param src_pixels Integer vector (source pixel buffer, row-major)
/// @param src_ncol Integer source buffer width
/// @param src_nrow Integer source buffer height
/// @param src_col_off Integer source buffer column offset in full image
/// @param src_row_off Integer source buffer row offset in full image
/// @param nodata Integer nodata value
/// @return Integer vector of warped pixels (row-major)
/// @export
#[extendr]
fn rust_warp_scanline(
    src_crs: &str,
    src_gt: &[f64],
    dst_crs: &str,
    dst_gt: &[f64],
    dst_dim: &[i32],
    src_pixels: &[i32],
    src_ncol: i32,
    src_nrow: i32,
    src_col_off: i32,
    src_row_off: i32,
    nodata: i32,
) -> extendr_api::Result<Robj> {
    if src_gt.len() != 6 || dst_gt.len() != 6 {
        return Err(Error::Other("geotransform must be length 6".to_string()));
    }
    if dst_dim.len() != 2 {
        return Err(Error::Other("dim must be length 2 (ncol, nrow)".to_string()));
    }

    let src_gt_arr: [f64; 6] = src_gt.try_into()
        .map_err(|_| Error::Other("src_gt conversion failed".to_string()))?;
    let dst_gt_arr: [f64; 6] = dst_gt.try_into()
        .map_err(|_| Error::Other("dst_gt conversion failed".to_string()))?;

    let transformer = transform::GenImgProjTransformer::new(
        src_crs, src_gt_arr, dst_crs, dst_gt_arr,
    )
    .map_err(Error::Other)?;

    let output = warp::warp_nearest(
        &transformer,
        src_pixels,
        src_ncol as usize,
        src_nrow as usize,
        src_col_off as usize,
        src_row_off as usize,
        dst_dim[0] as usize,
        dst_dim[1] as usize,
        nodata,
    );

    Ok(r!(output))
}

/// Warp source pixels using ApproxTransformer (GDAL-style adaptive interpolation).
///
/// Matches GDAL's actual default architecture (gdalwarp_lib.cpp line 3260):
/// the ApproxTransformer wraps the ENTIRE GenImgProjTransform. The max_error
/// parameter is in source pixel units (output of GenImgProjTransform when
/// called with bDstToSrc=TRUE). GDAL default: 0.125 pixels.
///
/// @param src_crs Character CRS string
/// @param src_gt Numeric vector length 6
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector length 6
/// @param dst_dim Integer vector c(ncol, nrow)
/// @param src_pixels Integer vector (source pixel buffer, row-major)
/// @param src_ncol Integer source buffer width
/// @param src_nrow Integer source buffer height
/// @param src_col_off Integer source buffer column offset in full image
/// @param src_row_off Integer source buffer row offset in full image
/// @param nodata Integer nodata value
/// @param max_error Numeric max interpolation error in pixels (default 0.125)
/// @return Integer vector of warped pixels (row-major)
/// @export
#[extendr]
fn rust_warp_approx(
    src_crs: &str,
    src_gt: &[f64],
    dst_crs: &str,
    dst_gt: &[f64],
    dst_dim: &[i32],
    src_pixels: &[i32],
    src_ncol: i32,
    src_nrow: i32,
    src_col_off: i32,
    src_row_off: i32,
    nodata: i32,
    max_error: f64,
) -> extendr_api::Result<Robj> {
    if src_gt.len() != 6 || dst_gt.len() != 6 {
        return Err(Error::Other("geotransform must be length 6".to_string()));
    }
    if dst_dim.len() != 2 {
        return Err(Error::Other("dim must be length 2 (ncol, nrow)".to_string()));
    }

    let src_gt_arr: [f64; 6] = src_gt.try_into()
        .map_err(|_| Error::Other("src_gt conversion failed".to_string()))?;
    let dst_gt_arr: [f64; 6] = dst_gt.try_into()
        .map_err(|_| Error::Other("dst_gt conversion failed".to_string()))?;

    // gdalwarp default architecture (gdalwarp_lib.cpp line 3260):
    // ApproxTransformer wraps the ENTIRE GenImgProjTransform.
    // max_error is in source pixel units (output space of GenImgProjTransform
    // when called with bDstToSrc=TRUE). Default: 0.125 pixels.
    let exact = transform::GenImgProjTransformer::new(
        src_crs, src_gt_arr, dst_crs, dst_gt_arr,
    )
    .map_err(Error::Other)?;

    let transformer = approx::ApproxTransformer::new(exact, max_error);

    let output = warp::warp_resample(
        &transformer,
        src_pixels,
        src_ncol as usize,
        src_nrow as usize,
        src_col_off as usize,
        src_row_off as usize,
        dst_dim[0] as usize,
        dst_dim[1] as usize,
        nodata,
        warp::ResampleAlg::NearestNeighbour,
    );

    Ok(r!(output))
}

// ===========================================================================
// extendr wrappers — source window computation
// ===========================================================================

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

/// Warp source pixels with selectable resampling algorithm.
///
/// Same as rust_warp_approx but with a `resample` parameter:
///   "near" = nearest neighbour (default)
///   "bilinear" = bilinear (2×2)
///   "cubic" = cubic convolution (4×4, Catmull-Rom)
///   "lanczos" = Lanczos windowed sinc (6×6)
///
/// @param src_crs Character CRS string
/// @param src_gt Numeric vector length 6
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector length 6
/// @param dst_dim Integer vector c(ncol, nrow)
/// @param src_pixels Integer vector (source pixel buffer, row-major)
/// @param src_ncol Integer source buffer width
/// @param src_nrow Integer source buffer height
/// @param src_col_off Integer source buffer column offset in full image
/// @param src_row_off Integer source buffer row offset in full image
/// @param nodata Integer nodata value
/// @param max_error Numeric max interpolation error in pixels (default 0.125)
/// @param resample Character resampling method: "near", "bilinear", "cubic", "lanczos"
/// @return Integer vector of warped pixels (row-major)
/// @export
#[extendr]
fn rust_warp_resample(
    src_crs: &str,
    src_gt: &[f64],
    dst_crs: &str,
    dst_gt: &[f64],
    dst_dim: &[i32],
    src_pixels: &[i32],
    src_ncol: i32,
    src_nrow: i32,
    src_col_off: i32,
    src_row_off: i32,
    nodata: i32,
    max_error: f64,
    resample: &str,
) -> extendr_api::Result<Robj> {
    if src_gt.len() != 6 || dst_gt.len() != 6 {
        return Err(Error::Other("geotransform must be length 6".to_string()));
    }
    if dst_dim.len() != 2 {
        return Err(Error::Other("dim must be length 2 (ncol, nrow)".to_string()));
    }

    let alg = match resample {
        "near" | "nearest" => warp::ResampleAlg::NearestNeighbour,
        "bilinear" => warp::ResampleAlg::Bilinear,
        "cubic" => warp::ResampleAlg::Cubic,
        "lanczos" => warp::ResampleAlg::Lanczos,
        other => return Err(Error::Other(
            format!("unknown resample method '{}', use: near, bilinear, cubic, lanczos", other)
        )),
    };

    let src_gt_arr: [f64; 6] = src_gt.try_into()
        .map_err(|_| Error::Other("src_gt conversion failed".to_string()))?;
    let dst_gt_arr: [f64; 6] = dst_gt.try_into()
        .map_err(|_| Error::Other("dst_gt conversion failed".to_string()))?;

    let exact = transform::GenImgProjTransformer::new(
        src_crs, src_gt_arr, dst_crs, dst_gt_arr,
    )
    .map_err(Error::Other)?;

    let transformer = approx::ApproxTransformer::new(exact, max_error);

    let output = warp::warp_resample(
        &transformer,
        src_pixels,
        src_ncol as usize,
        src_nrow as usize,
        src_col_off as usize,
        src_row_off as usize,
        dst_dim[0] as usize,
        dst_dim[1] as usize,
        nodata,
        alg,
    );

    Ok(r!(output))
}

// ===========================================================================
// Legacy wrappers (kept for backward compatibility with test scripts)
// ===========================================================================

/// Compute a warp mapping (legacy interface, wraps GenImgProjTransformer).
/// @param src_crs Character CRS string
/// @param src_gt Numeric vector length 6
/// @param src_dim Integer vector c(ncol, nrow)
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector length 6
/// @param dst_dim Integer vector c(ncol, nrow)
/// @return List with `src_cols`, `src_rows` integer vectors (1-based, NA for OOB)
/// @export
#[extendr]
fn rust_warp_map(
    src_crs: &str,
    src_gt: &[f64],
    src_dim: &[i32],
    dst_crs: &str,
    dst_gt: &[f64],
    dst_dim: &[i32],
) -> extendr_api::Result<Robj> {
    if src_gt.len() != 6 || dst_gt.len() != 6 {
        return Err(Error::Other("geotransform must be length 6".to_string()));
    }
    if src_dim.len() != 2 || dst_dim.len() != 2 {
        return Err(Error::Other("dim must be length 2 (ncol, nrow)".to_string()));
    }

    let src_gt_arr: [f64; 6] = src_gt.try_into()
        .map_err(|_| Error::Other("src_gt conversion failed".to_string()))?;
    let dst_gt_arr: [f64; 6] = dst_gt.try_into()
        .map_err(|_| Error::Other("dst_gt conversion failed".to_string()))?;

    let src_dim_arr = [src_dim[0] as usize, src_dim[1] as usize];
    let dst_dim_arr = [dst_dim[0] as usize, dst_dim[1] as usize];

    // Use GenImgProjTransformer + transform_grid
    let transformer = transform::GenImgProjTransformer::new(
        src_crs, src_gt_arr, dst_crs, dst_gt_arr,
    )
    .map_err(Error::Other)?;

    let (src_x, src_y) = transform::transform_grid(&transformer, &dst_dim_arr);

    // Convert fractional pixel coords to 1-based integer indices for R
    // NaN → NA, round for nearest neighbour, bounds check
    let src_cols_r: Vec<Option<i32>> = src_x.iter().enumerate().map(|(_, &sx)| {
        if sx.is_nan() { return None; }
        let col = sx as i32; // truncate (GDAL convention)
        if col >= 0 && col < src_dim_arr[0] as i32 {
            Some(col + 1)  // 1-based for R
        } else {
            None
        }
    }).collect();

    let src_rows_r: Vec<Option<i32>> = src_y.iter().enumerate().map(|(_, &sy)| {
        if sy.is_nan() { return None; }
        let row = sy as i32;
        if row >= 0 && row < src_dim_arr[1] as i32 {
            Some(row + 1)
        } else {
            None
        }
    }).collect();

    let result = list!(src_cols = src_cols_r, src_rows = src_rows_r);
    Ok(result.into())
}

/// Apply a warp map to source pixels (legacy interface).
/// @param src_pixels Integer vector of source pixel values (row-major)
/// @param src_ncol Integer source width
/// @param src_nrow Integer source height
/// @param src_cols Integer vector of source column indices (1-based, NA for OOB)
/// @param src_rows Integer vector of source row indices (1-based, NA for OOB)
/// @param nodata Integer nodata value
/// @return Integer vector of destination pixel values
/// @export
#[extendr]
fn rust_apply_warp(
    src_pixels: &[i32],
    src_ncol: i32,
    src_nrow: i32,
    src_cols: &[i32],
    src_rows: &[i32],
    nodata: i32,
) -> extendr_api::Result<Robj> {
    let n = src_cols.len();
    if src_rows.len() != n {
        return Err(Error::Other(
            "src_cols and src_rows must be same length".to_string(),
        ));
    }

    let src_dim = [src_ncol as usize, src_nrow as usize];
    let mut dst_pixels = vec![nodata; n];
    let na_int = i32::MIN;

    for i in 0..n {
        let c = src_cols[i];
        let r = src_rows[i];
        if c != na_int && r != na_int && c > 0 && r > 0 {
            let sc = (c - 1) as usize;
            let sr = (r - 1) as usize;
            let src_idx = sr * src_dim[0] + sc;
            if src_idx < src_pixels.len() {
                dst_pixels[i] = src_pixels[src_idx];
            }
        }
    }

    Ok(r!(dst_pixels))
}

// ===========================================================================
// Module registration
// ===========================================================================

extendr_module! {
    mod cogcache;
    fn rust_decode_tile;
    fn rust_fetch_decode_tile;
    fn rust_gen_img_proj_transform;
    fn rust_warp_scanline;
    fn rust_warp_approx;
    fn rust_warp_resample;
    fn rust_compute_source_window;
    fn rust_warp_map;
    fn rust_apply_warp;
}



#[cfg(test)]
mod tests {
    #[test]
    fn test_antarctic_crs() {
        let proj = proj::Proj::new_known_crs("EPSG:3857", "EPSG:3031", None).unwrap();
        let result = proj.convert((7006000.0, -10340500.0)).unwrap();
        println!("3031: {:?}", result);

        let proj2 = proj::Proj::new_known_crs("EPSG:3857", "EPSG:3412", None).unwrap();
        let result2 = proj2.convert((7006000.0, -10340500.0)).unwrap();
        println!("3412: {:?}", result2);
    }
}
