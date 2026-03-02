//! Grid specification and geotransform arithmetic.
//!
//! Mirrors the semantics of the R `vaster` package:
//! - Extent is (xmin, xmax, ymin, ymax)
//! - Dimensions are (ncol, nrow)
//! - Cell indices are 0-based internally (1-based at the R boundary)
//! - Geotransform follows GDAL convention: [originX, pixelW, 0, originY, 0, pixelH]
//!
//! Reference: GDAL geotransform tutorial
//! https://gdal.org/tutorials/geotransforms_tut.html

/// A regular raster grid: CRS + extent + dimensions.
/// Everything else (resolution, geotransform, cell centres) derives from these.
#[derive(Debug, Clone)]
pub struct GridSpec {
    pub crs: String,          // WKT, EPSG:XXXX, or PROJ string
    pub extent: [f64; 4],     // xmin, xmax, ymin, ymax (vaster convention)
    pub dim: [usize; 2],      // ncol, nrow
}

impl GridSpec {
    pub fn new(crs: &str, extent: [f64; 4], dim: [usize; 2]) -> Self {
        Self {
            crs: crs.to_string(),
            extent,
            dim,
        }
    }

    /// Build the GDAL geotransform from extent and dimensions.
    ///
    /// gt[0] = xmin (left edge)
    /// gt[1] = pixel width (positive)
    /// gt[2] = 0 (no rotation)
    /// gt[3] = ymax (top edge)
    /// gt[4] = 0 (no rotation)
    /// gt[5] = -pixel height (negative for north-up)
    pub fn geotransform(&self) -> [f64; 6] {
        extent_dim_to_gt(&self.extent, &self.dim)
    }

    /// Cell centre x-coordinate from 0-based column index.
    pub fn x_from_col(&self, col: f64) -> f64 {
        let gt = self.geotransform();
        x_from_col(&gt, col)
    }

    /// Cell centre y-coordinate from 0-based row index.
    pub fn y_from_row(&self, row: f64) -> f64 {
        let gt = self.geotransform();
        y_from_row(&gt, row)
    }

    /// 0-based fractional column from x-coordinate.
    pub fn col_from_x(&self, x: f64) -> f64 {
        let gt = self.geotransform();
        col_from_x(&gt, x)
    }

    /// 0-based fractional row from y-coordinate.
    pub fn row_from_y(&self, y: f64) -> f64 {
        let gt = self.geotransform();
        row_from_y(&gt, y)
    }

    /// Cell centre (x, y) from 0-based cell index (row-major).
    pub fn xy_from_cell(&self, cell: usize) -> (f64, f64) {
        let col = cell % self.dim[0];
        let row = cell / self.dim[0];
        (self.x_from_col(col as f64), self.y_from_row(row as f64))
    }

    /// Number of cells.
    pub fn ncells(&self) -> usize {
        self.dim[0] * self.dim[1]
    }
}

// ---------------------------------------------------------------------------
// Free functions operating on raw geotransforms
// ---------------------------------------------------------------------------

/// Cell centre x-coordinate from 0-based column index.
pub fn x_from_col(gt: &[f64; 6], col: f64) -> f64 {
    gt[0] + (col + 0.5) * gt[1]
}

/// Cell centre y-coordinate from 0-based row index.
pub fn y_from_row(gt: &[f64; 6], row: f64) -> f64 {
    gt[3] + (row + 0.5) * gt[5]
}

/// 0-based fractional column from x-coordinate.
/// Floor the result for nearest-left pixel index.
pub fn col_from_x(gt: &[f64; 6], x: f64) -> f64 {
    (x - gt[0]) / gt[1] - 0.5
}

/// 0-based fractional row from y-coordinate.
/// Floor the result for nearest-top pixel index.
pub fn row_from_y(gt: &[f64; 6], y: f64) -> f64 {
    (y - gt[3]) / gt[5] - 0.5
}

/// Build a geotransform from extent (xmin, xmax, ymin, ymax) and dimensions (ncol, nrow).
pub fn extent_dim_to_gt(extent: &[f64; 4], dim: &[usize; 2]) -> [f64; 6] {
    let pixel_width = (extent[1] - extent[0]) / dim[0] as f64;
    let pixel_height = -(extent[3] - extent[2]) / dim[1] as f64;
    [extent[0], pixel_width, 0.0, extent[3], 0.0, pixel_height]
}

/// Invert a geotransform (for non-rotated grids only: gt[2] == 0 && gt[4] == 0).
/// Returns the inverse geotransform that maps geo coords back to pixel coords.
///
/// GDAL equivalent: GDALInvGeoTransform()
pub fn inv_geotransform(gt: &[f64; 6]) -> Option<[f64; 6]> {
    // For non-rotated case (gt[2] == 0, gt[4] == 0):
    // inv[0] = -gt[0] / gt[1]
    // inv[1] = 1.0 / gt[1]
    // inv[3] = -gt[3] / gt[5]
    // inv[5] = 1.0 / gt[5]
    if gt[1].abs() < 1e-15 || gt[5].abs() < 1e-15 {
        return None;
    }
    Some([
        -gt[0] / gt[1],
        1.0 / gt[1],
        0.0,
        -gt[3] / gt[5],
        0.0,
        1.0 / gt[5],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_col() {
        let gt = [100.0, 10.0, 0.0, 200.0, 0.0, -10.0];
        let x = x_from_col(&gt, 5.0); // cell centre of column 5
        let c = col_from_x(&gt, x);
        assert!((c - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_inv_geotransform() {
        let gt = [100.0, 10.0, 0.0, 500.0, 0.0, -10.0];
        let inv = inv_geotransform(&gt).unwrap();
        // pixel (0,0) centre is at geo (105, 495)
        let x = 105.0;
        let y = 495.0;
        // inv maps geo → pixel
        let col = inv[0] + x * inv[1]; // should be 0.5 (centre of pixel 0)
        let row = inv[3] + y * inv[5]; // should be 0.5
        assert!((col - 0.5).abs() < 1e-10);
        assert!((row - 0.5).abs() < 1e-10);
    }
}
