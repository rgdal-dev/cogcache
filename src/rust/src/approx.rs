//! Approximate transformer: adaptive interpolation of coordinate transforms.
//!
//! Maps GDAL's `GDALApproxTransform` (gdaltransformer.cpp ~line 3500).
//!
//! Instead of calling PROJ for every pixel on a scanline, transform only a few
//! points exactly and linearly interpolate the rest. For smooth transforms
//! (UTM ↔ WebMercator), this reduces PROJ calls from ~256 to ~3–9 per scanline.
//!
//! ## Algorithm
//!
//! Given N points to transform:
//!   1. Transform 3 points exactly: first, middle, last
//!   2. Check: does the true midpoint match the linear interpolation?
//!   3. If error < max_error → linearly interpolate all interior points
//!   4. If error ≥ max_error → split at midpoint, recurse on each half
//!
//! Error is measured as L1 (Manhattan) distance in the output coordinate space
//! (pixels for GenImgProjTransformer). GDAL's default threshold is 0.125 pixels.
//!
//! ## GDAL reference
//!
//! `GDALApproxTransform()` in gdaltransformer.cpp (commit a3b7b01d3e):
//!   - Creates wrapper: `GDALCreateApproxTransformer()` ~line 3440
//!   - Transform function: `GDALApproxTransform()` ~line 3500
//!   - Error check: L1 norm of (interpolated_mid - true_mid) ~line 3570
//!   - Recursion: splits at `nMiddle = (nPoints-1)/2` ~line 3600

use crate::transform::Transformer;

/// Approximate transformer: wraps an inner (exact) transformer with adaptive
/// linear interpolation.
///
/// `max_error` is in the output coordinate units of the inner transformer.
/// For GenImgProjTransformer, that's pixels. GDAL default: 0.125 pixels.
pub struct ApproxTransformer<T: Transformer> {
    inner: T,
    max_error: f64,
}

impl<T: Transformer> ApproxTransformer<T> {
    /// Create an approximate transformer wrapping an exact one.
    ///
    /// # Arguments
    /// * `inner` - The exact transformer (e.g. GenImgProjTransformer)
    /// * `max_error` - Maximum interpolation error. GDAL default: 0.125
    ///
    /// # GDAL equivalent
    /// `GDALCreateApproxTransformer(pfnBaseTransformer, pBaseCBData, dfMaxError)`
    pub fn new(inner: T, max_error: f64) -> Self {
        Self { inner, max_error }
    }

    /// Access the inner transformer (e.g. for direct exact transforms).
    pub fn inner(&self) -> &T {
        &self.inner
    }

    /// Recursive approximation on the range [start .. start+len) of the arrays.
    ///
    /// After this call, x[start..start+len] and y[start..start+len] contain
    /// transformed coordinates, and success[start..start+len] are set.
    fn recurse(
        &self,
        dst_to_src: bool,
        x: &mut [f64],
        y: &mut [f64],
        success: &mut [bool],
        start: usize,
        len: usize,
    ) {
        // Base case: very small segments get exact transforms.
        //
        // GDAL uses nPoints < 3 as its base case (2 or fewer points).
        // We match this to get identical recursion depth.
        if len < 3 {
            let ok = self.inner.transform(
                dst_to_src,
                &mut x[start..start + len],
                &mut y[start..start + len],
            );
            for i in 0..len {
                success[start + i] = ok[i];
            }
            return;
        }

        let mid = start + (len - 1) / 2;
        let end = start + len - 1;

        // Transform 3 control points exactly (copies — originals preserved)
        let mut tx = [x[start], x[mid], x[end]];
        let mut ty = [y[start], y[mid], y[end]];
        let ok3 = self.inner.transform(dst_to_src, &mut tx, &mut ty);

        // If any control point fails, fall back to exact transform of all points.
        // GDAL does the same — if the endpoints can't be transformed, interpolation
        // is meaningless.
        if !ok3[0] || !ok3[1] || !ok3[2] {
            let ok = self.inner.transform(
                dst_to_src,
                &mut x[start..start + len],
                &mut y[start..start + len],
            );
            for i in 0..len {
                success[start + i] = ok[i];
            }
            return;
        }

        // Check interpolation error at midpoint.
        //
        // GDAL (gdaltransformer.cpp ~line 3570):
        //   dfError = fabs((x2[0] + x2[2]) / 2.0 - x2[1])
        //           + fabs((y2[0] + y2[2]) / 2.0 - y2[1]);
        let interp_x = (tx[0] + tx[2]) * 0.5;
        let interp_y = (ty[0] + ty[2]) * 0.5;
        let error = (interp_x - tx[1]).abs() + (interp_y - ty[1]).abs();

        if error <= self.max_error {
            // Error acceptable: linearly interpolate all points.
            //
            // GDAL (gdaltransformer.cpp ~line 3585):
            //   for (i = nPoints - 1; i >= 0; i--) {
            //       dfRatio = i / (double)(nPoints - 1);
            //       x[i] = x2[0] * (1.0 - dfRatio) + x2[2] * dfRatio;
            //   }
            let n_f = (len - 1) as f64;
            for i in 0..len {
                let ratio = i as f64 / n_f;
                x[start + i] = tx[0] + (tx[2] - tx[0]) * ratio;
                y[start + i] = ty[0] + (ty[2] - ty[0]) * ratio;
                success[start + i] = true;
            }
        } else {
            // Error too large: split at midpoint and recurse.
            //
            // Both halves share the midpoint index. The left half will
            // transform x[mid]/y[mid] in place. The right half needs
            // the ORIGINAL (untransformed) value at mid as its starting
            // point. Save it before the left recursion overwrites it.
            //
            // GDAL avoids this because its recursive calls operate on
            // separate pointer offsets and re-transform their own
            // endpoints via the 3-point test. But at the base case
            // (len < 3), the shared point IS read directly from the array.
            let orig_mid_x = x[mid];
            let orig_mid_y = y[mid];

            let left_len = mid - start + 1;
            self.recurse(dst_to_src, x, y, success, start, left_len);

            // Restore the original midpoint for the right half
            x[mid] = orig_mid_x;
            y[mid] = orig_mid_y;

            let right_len = end - mid + 1;
            self.recurse(dst_to_src, x, y, success, mid, right_len);
        }
    }
}

impl<T: Transformer> Transformer for ApproxTransformer<T> {
    fn transform(
        &self,
        dst_to_src: bool,
        x: &mut [f64],
        y: &mut [f64],
    ) -> Vec<bool> {
        let n = x.len();

        // Bypass approximation if disabled
        if self.max_error <= 0.0 {
            return self.inner.transform(dst_to_src, x, y);
        }

        let mut success = vec![false; n];
        self.recurse(dst_to_src, x, y, &mut success, 0, n);
        success
    }
}
