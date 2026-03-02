use extendr_api::prelude::*;
use flate2::read::ZlibDecoder;
use std::io::Read;

// ---------------------------------------------------------------------------
// Grid arithmetic — mirrors vaster semantics (0-based internally)
// ---------------------------------------------------------------------------

/// Geotransform: [origin_x, pixel_width, 0, origin_y, 0, pixel_height]
/// pixel_height is typically negative (north-up).

/// Cell centre x-coordinate from 0-based column index.
fn x_from_col(gt: &[f64; 6], col: f64) -> f64 {
    gt[0] + (col + 0.5) * gt[1]
}

/// Cell centre y-coordinate from 0-based row index.
fn y_from_row(gt: &[f64; 6], row: f64) -> f64 {
    gt[3] + (row + 0.5) * gt[5]
}

/// 0-based column from x-coordinate (floor = nearest left).
fn col_from_x(gt: &[f64; 6], x: f64) -> f64 {
    (x - gt[0]) / gt[1] - 0.5
}

/// 0-based row from y-coordinate (floor = nearest top).
fn row_from_y(gt: &[f64; 6], y: f64) -> f64 {
    (y - gt[3]) / gt[5] - 0.5
}

/// Build a geotransform from extent (xmin, xmax, ymin, ymax) and dimensions (ncol, nrow).
fn extent_dim_to_gt(extent: &[f64; 4], dim: &[usize; 2]) -> [f64; 6] {
    let pixel_width = (extent[1] - extent[0]) / dim[0] as f64;
    let pixel_height = -(extent[3] - extent[2]) / dim[1] as f64;
    [extent[0], pixel_width, 0.0, extent[3], 0.0, pixel_height]
}

// ---------------------------------------------------------------------------
// Warp map — the core extracted module
// ---------------------------------------------------------------------------

/// A warp mapping: for each destination pixel, the source pixel coordinates.
/// src_cols and src_rows are 0-based. -1 means out-of-bounds.
/// Values are stored in row-major order (row 0 col 0, row 0 col 1, ...).
pub struct WarpMap {
    pub dst_dim: [usize; 2],   // (ncol, nrow)
    pub src_cols: Vec<i32>,
    pub src_rows: Vec<i32>,
}

/// Compute the warp mapping from a destination grid to a source grid.
///
/// For each destination pixel centre:
///   1. Compute its (x, y) in destination CRS
///   2. Transform to source CRS via PROJ
///   3. Convert to source pixel index (0-based, nearest neighbour)
///   4. Mark out-of-bounds as -1
fn compute_warp_map(
    src_crs: &str,
    src_gt: &[f64; 6],
    src_dim: &[usize; 2],
    dst_crs: &str,
    dst_gt: &[f64; 6],
    dst_dim: &[usize; 2],
) -> std::result::Result<WarpMap, String> {
    // Create PROJ transformer: dst CRS → src CRS
    let transformer = proj::Proj::new_known_crs(dst_crs, src_crs, None)
        .map_err(|e| format!("PROJ error creating transformer: {}", e))?;

    let n_pixels = dst_dim[0] * dst_dim[1];
    let mut src_cols = vec![-1i32; n_pixels];
    let mut src_rows = vec![-1i32; n_pixels];

    // Row-major iteration: row 0 col 0, row 0 col 1, ...
    for row in 0..dst_dim[1] {
        for col in 0..dst_dim[0] {
            let idx = row * dst_dim[0] + col;

            // Destination pixel centre in destination CRS
            let dst_x = x_from_col(dst_gt, col as f64);
            let dst_y = y_from_row(dst_gt, row as f64);

            // Transform to source CRS
            // proj crate normalises to (lon, lat) / (easting, northing)
            let src_coord = match transformer.convert((dst_x, dst_y)) {
                Ok(c) => c,
                Err(_) => continue, // transform failed, leave as -1
            };

            // Source pixel index (0-based)
            let sc = col_from_x(src_gt, src_coord.0);
            let sr = row_from_y(src_gt, src_coord.1);

            // Nearest neighbour: round to nearest pixel centre
            let sc_nn = sc.round() as i32;
            let sr_nn = sr.round() as i32;

            // Bounds check
            if sc_nn >= 0 && sc_nn < src_dim[0] as i32
                && sr_nn >= 0 && sr_nn < src_dim[1] as i32
            {
                src_cols[idx] = sc_nn;
                src_rows[idx] = sr_nn;
            }
        }
    }

    Ok(WarpMap {
        dst_dim: *dst_dim,
        src_cols,
        src_rows,
    })
}

/// Apply a warp map to a source pixel buffer, producing destination pixels.
/// src_pixels is row-major, src_dim[0] wide, src_dim[1] tall.
/// Returns row-major destination pixels, with nodata for unmapped pixels.
fn apply_warp_map(
    warp_map: &WarpMap,
    src_pixels: &[i32],
    src_dim: &[usize; 2],
    nodata: i32,
) -> Vec<i32> {
    let n = warp_map.dst_dim[0] * warp_map.dst_dim[1];
    let mut dst_pixels = vec![nodata; n];

    for i in 0..n {
        let sc = warp_map.src_cols[i];
        let sr = warp_map.src_rows[i];
        if sc >= 0 && sr >= 0 {
            let src_idx = sr as usize * src_dim[0] + sc as usize;
            if src_idx < src_pixels.len() {
                dst_pixels[i] = src_pixels[src_idx];
            }
        }
    }

    dst_pixels
}

// ---------------------------------------------------------------------------
// Tile decode (existing)
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// extendr wrappers
// ---------------------------------------------------------------------------

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

/// Compute a warp mapping from destination grid to source grid.
///
/// Returns a list with components `src_cols` and `src_rows` (integer vectors,
/// 0-based, -1 for out-of-bounds), in row-major order matching
/// `xy_from_cell(dst_dim, dst_ext, 1:n)`.
///
/// @param src_crs Character CRS string (WKT, EPSG:XXXX, or PROJ string)
/// @param src_gt Numeric vector of length 6 (GDAL geotransform)
/// @param src_dim Integer vector c(ncol, nrow)
/// @param dst_crs Character CRS string
/// @param dst_gt Numeric vector of length 6 (GDAL geotransform)
/// @param dst_dim Integer vector c(ncol, nrow)
/// @return List with `src_cols` and `src_rows` integer vectors
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

    let wmap = compute_warp_map(
        src_crs, &src_gt_arr, &src_dim_arr,
        dst_crs, &dst_gt_arr, &dst_dim_arr,
    )
    .map_err(Error::Other)?;

    // Return as R list with src_cols and src_rows
    // Convert to 1-based for R: -1 stays as NA
    let src_cols_r: Vec<Option<i32>> = wmap.src_cols.iter()
        .map(|&c| if c >= 0 { Some(c + 1) } else { None })
        .collect();
    let src_rows_r: Vec<Option<i32>> = wmap.src_rows.iter()
        .map(|&r| if r >= 0 { Some(r + 1) } else { None })
        .collect();

    let result = list!(src_cols = src_cols_r, src_rows = src_rows_r);
    Ok(result.into())
}

/// Warp source pixels through a precomputed warp map.
///
/// @param src_pixels Integer vector of source pixel values (row-major)
/// @param src_ncol Integer number of columns in source
/// @param src_nrow Integer number of rows in source
/// @param src_cols Integer vector of source column indices (1-based, NA for OOB)
/// @param src_rows Integer vector of source row indices (1-based, NA for OOB)
/// @param nodata Integer nodata value
/// @return Integer vector of destination pixel values (row-major)
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
        return Err(Error::Other("src_cols and src_rows must be same length".to_string()));
    }

    let src_dim = [src_ncol as usize, src_nrow as usize];
    let mut dst_pixels = vec![nodata; n];
    let na_int = i32::MIN; // R's NA_integer_

    for i in 0..n {
        let c = src_cols[i];
        let r = src_rows[i];
        if c != na_int && r != na_int && c > 0 && r > 0 {
            // Convert from R 1-based to 0-based
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

extendr_module! {
    mod cogcache;
    fn rust_decode_tile;
    fn rust_fetch_decode_tile;
    fn rust_warp_map;
    fn rust_apply_warp;
}
