//! Coordinate transformers for the warp pipeline.
//!
//! Maps GDAL's transformer hierarchy (gdaltransformer.cpp):
//!
//!   Transformer trait          ← GDALTransformerFunc callback signature
//!   GenImgProjTransformer      ← GDALGenImgProjTransform (~line 2800)
//!   ApproxTransformer          ← GDALApproxTransform (future Phase 2)
//!
//! The GenImgProjTransformer composes three steps:
//!   1. dst pixel → dst geo     (apply dst geotransform)
//!   2. dst geo → src geo       (PROJ coordinate transform)
//!   3. src geo → src pixel     (apply inverse src geotransform)
//!
//! This is the "back-transform coordinates to source" approach.
//! For dst_to_src=false, the steps are reversed.

use crate::grid::{inv_geotransform, x_from_col, y_from_row};

// ---------------------------------------------------------------------------
// Transformer trait
// ---------------------------------------------------------------------------

/// Trait for coordinate transformers.
///
/// Matches the semantics of GDAL's GDALTransformerFunc:
/// - Takes arrays of x, y coordinates (modified in place)
/// - Returns per-point success flags
/// - `dst_to_src=true`: input is dest pixel coords, output is src pixel coords
/// - `dst_to_src=false`: input is src pixel coords, output is dest pixel coords
///
/// Coordinates are in pixel/line space (0-based, +0.5 = cell centre).
pub trait Transformer {
    fn transform(
        &self,
        dst_to_src: bool,
        x: &mut [f64],
        y: &mut [f64],
    ) -> Vec<bool>;
}

// ---------------------------------------------------------------------------
// GenImgProjTransformer
// ---------------------------------------------------------------------------

/// The workhorse transformer: pixel ↔ georef ↔ CRS pipeline.
///
/// Equivalent to GDAL's GDALGenImgProjTransform (gdaltransformer.cpp).
///
/// Holds precomputed geotransforms and their inverses for both source
/// and destination grids, plus a PROJ transformer for the CRS step.
pub struct GenImgProjTransformer {
    // Source grid affine transforms
    src_gt: [f64; 6],
    src_inv_gt: [f64; 6],

    // Destination grid affine transforms
    dst_gt: [f64; 6],
    dst_inv_gt: [f64; 6],

    // PROJ: dst CRS → src CRS (for the dst_to_src direction)
    proj_dst_to_src: proj::Proj,

    // PROJ: src CRS → dst CRS (for the src_to_dst direction)
    proj_src_to_dst: proj::Proj,
}

impl GenImgProjTransformer {
    /// Create a new GenImgProjTransformer.
    ///
    /// # Arguments
    /// * `src_crs` - Source CRS (WKT, EPSG:XXXX, or PROJ string)
    /// * `src_gt` - Source geotransform (GDAL convention, 6 elements)
    /// * `dst_crs` - Destination CRS
    /// * `dst_gt` - Destination geotransform
    ///
    /// # GDAL equivalent
    /// GDALCreateGenImgProjTransformer2() with source and dest datasets.
    /// We skip the dataset handles and take the CRS + geotransform directly,
    /// which is what GDAL extracts from the datasets anyway.
    pub fn new(
        src_crs: &str,
        src_gt: [f64; 6],
        dst_crs: &str,
        dst_gt: [f64; 6],
    ) -> Result<Self, String> {
        let src_inv_gt = inv_geotransform(&src_gt)
            .ok_or_else(|| "Source geotransform is not invertible".to_string())?;
        let dst_inv_gt = inv_geotransform(&dst_gt)
            .ok_or_else(|| "Destination geotransform is not invertible".to_string())?;

        // PROJ transformers — one for each direction
        // new_known_crs normalises axis order to (easting, northing) / (lon, lat)
        let proj_dst_to_src = proj::Proj::new_known_crs(dst_crs, src_crs, None)
            .map_err(|e| format!("PROJ dst→src failed: {}", e))?;
        let proj_src_to_dst = proj::Proj::new_known_crs(src_crs, dst_crs, None)
            .map_err(|e| format!("PROJ src→dst failed: {}", e))?;

        Ok(Self {
            src_gt,
            src_inv_gt,
            dst_gt,
            dst_inv_gt,
            proj_dst_to_src,
            proj_src_to_dst,
        })
    }

    /// Create from GridSpec pairs (convenience).
    pub fn from_grids(
        src: &crate::grid::GridSpec,
        dst: &crate::grid::GridSpec,
    ) -> Result<Self, String> {
        Self::new(&src.crs, src.geotransform(), &dst.crs, dst.geotransform())
    }
}

impl Transformer for GenImgProjTransformer {
    /// Transform an array of coordinates.
    ///
    /// When `dst_to_src = true`:
    ///   input:  x[i], y[i] = destination pixel coords (col + 0.5, row + 0.5)
    ///   output: x[i], y[i] = source pixel coords (fractional)
    ///
    /// When `dst_to_src = false`:
    ///   input:  x[i], y[i] = source pixel coords
    ///   output: x[i], y[i] = destination pixel coords
    ///
    /// The three steps for dst_to_src:
    ///   1. Apply dst geotransform: pixel → geo
    ///   2. PROJ: dst CRS → src CRS
    ///   3. Apply inverse src geotransform: geo → pixel
    ///
    /// This matches GDAL's GDALGenImgProjTransform() exactly.
    fn transform(
        &self,
        dst_to_src: bool,
        x: &mut [f64],
        y: &mut [f64],
    ) -> Vec<bool> {
        let n = x.len();
        assert_eq!(n, y.len());
        let mut success = vec![true; n];

        if dst_to_src {
            // Step 1: dst pixel → dst geo (apply dst geotransform)
            //
            // GDAL (gdaltransformer.cpp ~line 2830):
            //   dfGeoX = gt[0] + dfPixel * gt[1] + dfLine * gt[2]
            //   dfGeoY = gt[3] + dfPixel * gt[4] + dfLine * gt[5]
            //
            // For non-rotated grids (gt[2]=0, gt[4]=0), this simplifies to
            // our x_from_col / y_from_row — but those add +0.5 for cell centre.
            // Here x[i] already IS pixel + 0.5, so we use the raw geotransform.
            for i in 0..n {
                let pixel = x[i];
                let line = y[i];
                x[i] = self.dst_gt[0] + pixel * self.dst_gt[1] + line * self.dst_gt[2];
                y[i] = self.dst_gt[3] + pixel * self.dst_gt[4] + line * self.dst_gt[5];
            }

            // Step 2: dst geo → src geo (PROJ)
            for i in 0..n {
                match self.proj_dst_to_src.convert((x[i], y[i])) {
                    Ok((sx, sy)) => {
                        x[i] = sx;
                        y[i] = sy;
                    }
                    Err(_) => {
                        success[i] = false;
                    }
                }
            }

            // Step 3: src geo → src pixel (apply inverse src geotransform)
            //
            // GDAL (gdaltransformer.cpp ~line 2870):
            //   dfPixel = inv_gt[0] + dfGeoX * inv_gt[1] + dfGeoY * inv_gt[2]
            //   dfLine  = inv_gt[3] + dfGeoX * inv_gt[4] + dfGeoY * inv_gt[5]
            for i in 0..n {
                if !success[i] {
                    continue;
                }
                let geo_x = x[i];
                let geo_y = y[i];
                x[i] = self.src_inv_gt[0]
                    + geo_x * self.src_inv_gt[1]
                    + geo_y * self.src_inv_gt[2];
                y[i] = self.src_inv_gt[3]
                    + geo_x * self.src_inv_gt[4]
                    + geo_y * self.src_inv_gt[5];
            }
        } else {
            // src_to_dst: reverse the three steps

            // Step 1: src pixel → src geo
            for i in 0..n {
                let pixel = x[i];
                let line = y[i];
                x[i] = self.src_gt[0] + pixel * self.src_gt[1] + line * self.src_gt[2];
                y[i] = self.src_gt[3] + pixel * self.src_gt[4] + line * self.src_gt[5];
            }

            // Step 2: src geo → dst geo (PROJ)
            for i in 0..n {
                match self.proj_src_to_dst.convert((x[i], y[i])) {
                    Ok((dx, dy)) => {
                        x[i] = dx;
                        y[i] = dy;
                    }
                    Err(_) => {
                        success[i] = false;
                    }
                }
            }

            // Step 3: dst geo → dst pixel
            for i in 0..n {
                if !success[i] {
                    continue;
                }
                let geo_x = x[i];
                let geo_y = y[i];
                x[i] = self.dst_inv_gt[0]
                    + geo_x * self.dst_inv_gt[1]
                    + geo_y * self.dst_inv_gt[2];
                y[i] = self.dst_inv_gt[3]
                    + geo_x * self.dst_inv_gt[4]
                    + geo_y * self.dst_inv_gt[5];
            }
        }

        success
    }
}

// ---------------------------------------------------------------------------
// Convenience: transform a whole scanline
// ---------------------------------------------------------------------------

/// Build destination pixel coordinates for a scanline and transform to source.
///
/// This is what GDAL does at the top of GWKNearestThread / GWKGeneralCaseThread:
///   for (iDstX = 0; iDstX < nDstXSize; iDstX++) {
///       padfX[iDstX] = iDstX + 0.5 + nDstXOff;
///       padfY[iDstX] = iDstY + 0.5 + nDstYOff;
///   }
///   pfnTransformer(..., TRUE, nDstXSize, padfX, padfY, padfZ, pabSuccess);
///
/// Returns (src_x, src_y, success) arrays.
pub fn transform_scanline(
    transformer: &impl Transformer,
    dst_y: usize,        // destination row (0-based)
    dst_ncol: usize,     // destination width
    dst_x_off: usize,    // destination column offset (usually 0)
    dst_y_off: usize,    // destination row offset (usually 0)
) -> (Vec<f64>, Vec<f64>, Vec<bool>) {
    let mut x: Vec<f64> = (0..dst_ncol)
        .map(|col| (col + dst_x_off) as f64 + 0.5)
        .collect();
    let mut y = vec![(dst_y + dst_y_off) as f64 + 0.5; dst_ncol];

    let success = transformer.transform(true, &mut x, &mut y);

    (x, y, success)
}

// ---------------------------------------------------------------------------
// Convenience: full grid transform (replaces old compute_warp_map)
// ---------------------------------------------------------------------------

/// Compute source pixel coordinates for every destination pixel.
///
/// Returns parallel vectors of (src_col, src_row) as fractional f64.
/// Failed transforms are NaN.
///
/// This is equivalent to calling transform_scanline for every row
/// and collecting the results.
pub fn transform_grid(
    transformer: &impl Transformer,
    dst_dim: &[usize; 2],
) -> (Vec<f64>, Vec<f64>) {
    let n = dst_dim[0] * dst_dim[1];
    let mut all_x = Vec::with_capacity(n);
    let mut all_y = Vec::with_capacity(n);

    for row in 0..dst_dim[1] {
        let (sx, sy, ok) = transform_scanline(transformer, row, dst_dim[0], 0, 0);
        for col in 0..dst_dim[0] {
            if ok[col] {
                all_x.push(sx[col]);
                all_y.push(sy[col]);
            } else {
                all_x.push(f64::NAN);
                all_y.push(f64::NAN);
            }
        }
    }

    (all_x, all_y)
}
