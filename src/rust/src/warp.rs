//! Warp kernel: the actual pixel-shuffling loop.
//!
//! Maps GDAL's GWKNearestThread / GWKGeneralCaseThread
//! (gdalwarpkernel.cpp ~line 5510).
//!
//! The pattern:
//!   for each output scanline:
//!     build array of destination pixel coordinates
//!     call transformer to map dst→src
//!     for each output pixel:
//!       resample source pixels at the transformed coordinate
//!
//! We separate the "resample" step as a trait for future bilinear/cubic/lanczos.
//! For now, nearest neighbour only.

use crate::transform::{transform_scanline, Transformer};

/// Nearest-neighbour warp: for each destination pixel, find the closest source pixel.
///
/// # Arguments
/// * `transformer` - maps dst pixel coords → src pixel coords
/// * `src_pixels` - source pixel buffer (row-major, `src_ncol` wide)
/// * `src_ncol` - source buffer width
/// * `src_nrow` - source buffer height
/// * `src_col_off` - source buffer column offset in the full source image
/// * `src_row_off` - source buffer row offset in the full source image
/// * `dst_ncol` - destination width
/// * `dst_nrow` - destination height
/// * `nodata` - fill value for unmapped pixels
///
/// # Returns
/// Row-major destination pixel buffer.
///
/// # GDAL equivalent
/// GWKNearestThread() in gdalwarpkernel.cpp, simplified for:
/// - No mask/density handling
/// - No multi-band (single band)
/// - Source buffer is a window into the full image (offsets matter)
pub fn warp_nearest(
    transformer: &impl Transformer,
    src_pixels: &[i32],
    src_ncol: usize,
    src_nrow: usize,
    src_col_off: usize,
    src_row_off: usize,
    dst_ncol: usize,
    dst_nrow: usize,
    nodata: i32,
) -> Vec<i32> {
    let mut output = vec![nodata; dst_ncol * dst_nrow];

    for dst_row in 0..dst_nrow {
        // Build scanline coords and transform dst→src
        // x[i] = col + 0.5, y[i] = row + 0.5 (cell centres)
        // After transform: x[i], y[i] = fractional src pixel coords
        let (src_x, src_y, ok) = transform_scanline(
            transformer, dst_row, dst_ncol, 0, 0,
        );

        for dst_col in 0..dst_ncol {
            if !ok[dst_col] {
                continue;
            }

            // Nearest neighbour: GDAL uses truncation (floor for positive values)
            // static_cast<int>(padfX[iDstX]) — this truncates toward zero
            //
            // The transformer output is in full-image pixel coords.
            // Adjust for the source buffer window offset.
            let src_col_f = src_x[dst_col];
            let src_row_f = src_y[dst_col];

            // Truncate to integer pixel index (GDAL convention)
            let src_col_img = src_col_f as i64;
            let src_row_img = src_row_f as i64;

            // Convert from full-image coords to buffer-relative coords
            let buf_col = src_col_img - src_col_off as i64;
            let buf_row = src_row_img - src_row_off as i64;

            // Bounds check against the source buffer
            if buf_col < 0
                || buf_col >= src_ncol as i64
                || buf_row < 0
                || buf_row >= src_nrow as i64
            {
                continue;
            }

            let src_idx = buf_row as usize * src_ncol + buf_col as usize;
            let dst_idx = dst_row * dst_ncol + dst_col;

            let val = src_pixels[src_idx];
            if val != nodata {
                output[dst_idx] = val;
            }
        }
    }

    output
}
