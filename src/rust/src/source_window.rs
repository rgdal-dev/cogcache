//! Compute the source pixel window needed for a given destination window.
//!
//! Maps GDAL's `GDALWarpOperation::ComputeSourceWindow()` and
//! `ComputeSourceWindowTransformPoints()` (gdalwarpoperation.cpp,
//! commit a3b7b01d3e, GDAL 3.13.0dev, lines 2751–3367).
//!
//! ## Architecture
//!
//! This is a pure function of a `Transformer` — no dataset handles, no
//! warp options struct, no mutable state. It answers:
//!
//! > "Given a destination pixel window and a coordinate transformer,
//! >  which source pixels are needed?"
//!
//! The result is a `SourceWindow` (xoff, yoff, xsize, ysize) in source
//! pixel coordinates, plus diagnostic fields (fill_ratio, n_failed).
//!
//! ## Algorithm (matching GDAL)
//!
//! 1. **Edge sampling** (default): sample `step_count * 4` points around
//!    the perimeter of the destination window.
//! 2. **Corner check**: if any of the 4 dest corners fail to transform,
//!    switch to ALL mode (every pixel on every edge).
//! 3. **Grid fallback**: if any edge points fail, retry with a
//!    `(step_count+2)²` grid through the interior.
//! 4. **Bounds clamping**: clamp to source raster, add resampling padding.
//! 5. **Antimeridian heuristic**: if source window covers >90% of source
//!    width, snap to full width (GDAL line 3325).
//! 6. **Fill ratio**: ratio of clamped to unclamped window area — low
//!    values suggest the dest window should be split.
//!
//! ## What we skip (vs GDAL)
//!
//! - `CHECK_WITH_INVERT_PROJ` retry (PROJ-4 era workaround)
//! - `ComputeSourceWindowStartingFromSource` forward probe (refinement)
//! - `aDstXYSpecialPoints` pole detection (could add later)
//! - `SAMPLE_GRID` / `SAMPLE_STEPS` warp options (we use fixed defaults)

use crate::transform::Transformer;

/// Result of compute_source_window.
#[derive(Debug, Clone)]
pub struct SourceWindow {
    /// Column offset in source raster (0-based)
    pub xoff: i32,
    /// Row offset in source raster (0-based)
    pub yoff: i32,
    /// Width in source pixels
    pub xsize: i32,
    /// Height in source pixels
    pub ysize: i32,
    /// Ratio of clamped window area to unclamped — low values mean
    /// much of the source window falls outside the raster.
    pub fill_ratio: f64,
    /// Number of sample points that failed to transform.
    pub n_failed: i32,
    /// Total sample points attempted.
    pub n_samples: i32,
}

/// Default step count for edge/grid sampling (GDAL: DEFAULT_STEP_COUNT = 21).
const DEFAULT_STEP_COUNT: usize = 21;

/// Compute the source pixel window needed for a destination window.
///
/// # Arguments
///
/// * `transformer` — any `Transformer` (exact or approx)
/// * `dst_off` — destination window offset `[col_off, row_off]` in pixels
/// * `dst_size` — destination window size `[ncol, nrow]` in pixels
/// * `src_raster_size` — full source raster `[ncol, nrow]`
/// * `resample_padding` — extra pixels for resampling kernel
///   (0 = nearest, 1 = bilinear, 2 = cubic, 3 = lanczos)
///
/// Returns `None` if too few points transform successfully (<5 valid).
pub fn compute_source_window(
    transformer: &impl Transformer,
    dst_off: [i32; 2],
    dst_size: [i32; 2],
    src_raster_size: [i32; 2],
    resample_padding: i32,
) -> Option<SourceWindow> {
    let step_count = DEFAULT_STEP_COUNT;

    // --- Phase 1: Try edge-only sampling ---
    //
    // First check corners. If any fail, use ALL mode (every edge pixel).
    // (GDAL lines 3050–3073)
    let use_all = {
        let mut cx = vec![
            dst_off[0] as f64,
            (dst_off[0] + dst_size[0]) as f64,
            dst_off[0] as f64,
            (dst_off[0] + dst_size[0]) as f64,
        ];
        let mut cy = vec![
            dst_off[1] as f64,
            dst_off[1] as f64,
            (dst_off[1] + dst_size[1]) as f64,
            (dst_off[1] + dst_size[1]) as f64,
        ];
        let ok = transformer.transform(true, &mut cx, &mut cy);
        !ok[0] || !ok[1] || !ok[2] || !ok[3]
    };

    let (min_x, min_y, max_x, max_y, n_samples, n_failed) = if use_all {
        // ALL mode: every pixel on every edge (GDAL lines 2870–2896)
        sample_edges_all(transformer, dst_off, dst_size)
    } else {
        // Step mode: step_count points per edge (GDAL lines 2898–2922)
        sample_edges_step(transformer, dst_off, dst_size, step_count)
    };

    // --- Phase 2: Grid fallback if any edge points failed ---
    // (GDAL lines 3153–3164)
    let (min_x, min_y, max_x, max_y, n_samples, n_failed) = if n_failed > 0 {
        sample_grid(transformer, dst_off, dst_size, step_count)
    } else {
        (min_x, min_y, max_x, max_y, n_samples, n_failed)
    };

    // --- Phase 3: Check if too many failures ---
    // (GDAL line 3170: nFailedCount > nSamplePoints - 5)
    if n_failed > n_samples - 5 {
        return None;
    }

    // --- Phase 4: Early exit if completely outside source ---
    // (GDAL lines 3214–3228)
    let src_w = src_raster_size[0] as f64;
    let src_h = src_raster_size[1] as f64;
    if min_x > src_w || max_x < 0.0 || min_y > src_h || max_y < 0.0 {
        return Some(SourceWindow {
            xoff: 0, yoff: 0, xsize: 0, ysize: 0,
            fill_ratio: 0.0, n_failed, n_samples,
        });
    }

    // --- Phase 5: Round-if-close-enough ---
    // (GDAL lines 3232–3243)
    let snap = |v: f64| -> f64 {
        let r = v.round();
        if (r - v).abs() < 1e-6 { r } else { v }
    };
    let min_x = snap(min_x);
    let min_y = snap(min_y);
    let max_x = snap(max_x);
    let max_y = snap(max_y);

    // --- Phase 6: Resampling kernel padding ---
    // (GDAL lines 3262–3278)
    // Scale padding if downsampling significantly
    let dx_scale = (dst_size[0] as f64 / (max_x - min_x)).max(1e-3);
    let dy_scale = (dst_size[1] as f64 / (max_y - min_y)).max(1e-3);
    let mut x_radius = if dx_scale < 0.95 {
        (resample_padding as f64 / dx_scale).ceil() as i32
    } else {
        resample_padding
    };
    let mut y_radius = if dy_scale < 0.95 {
        (resample_padding as f64 / dy_scale).ceil() as i32
    } else {
        resample_padding
    };

    // Extra pixel beyond the kernel radius to ensure interpolation kernels
    // always have their full neighbourhood. Without this, transformed coords
    // near the buffer edge can put kernel pixels outside the buffer, forcing
    // the weight-renormalisation fallback path which produces different values
    // than GDAL. GDAL achieves this via nExtraSrcPixels (gdalwarpoperation.cpp
    // ~L3270). See: bilinear edge bug at pixel [1,157] in GEBCO diagnostics.
    if resample_padding > 0 {
        x_radius += 1;
        y_radius += 1;
    }

    // Extra padding if there were failures (GDAL lines 3293–3297)
    if n_failed > 0 {
        x_radius += 10;
        y_radius += 10;
    }

    // --- Phase 7: Clamp to source raster bounds ---
    // (GDAL lines 3309–3351)
    let min_x_clamped = min_x.max(0.0) as i32;
    let min_y_clamped = min_y.max(0.0) as i32;
    let max_x_clamped = max_x.ceil().min(src_w) as i32;
    let max_y_clamped = max_y.ceil().min(src_h) as i32;

    let src_x_size_raw = (max_x - min_x)
        .max(0.0)
        .min((src_raster_size[0] - min_x_clamped) as f64);
    let src_y_size_raw = (max_y - min_y)
        .max(0.0)
        .min((src_raster_size[1] - min_y_clamped) as f64);

    // Antimeridian heuristic: if >90% of source width, use full width
    // (GDAL lines 3325–3329)
    let (xoff, xsize) = if (max_x_clamped - min_x_clamped) as f64 > 0.9 * src_w {
        (0, src_raster_size[0])
    } else {
        let xoff = 0.max((min_x_clamped - x_radius).min(src_raster_size[0]));
        let xsize = 0.max(
            (max_x_clamped - xoff + x_radius).min(src_raster_size[0] - xoff),
        );
        (xoff, xsize)
    };

    let (yoff, ysize) = if (max_y_clamped - min_y_clamped) as f64 > 0.9 * src_h {
        (0, src_raster_size[1])
    } else {
        let yoff = 0.max((min_y_clamped - y_radius).min(src_raster_size[1]));
        let ysize = 0.max(
            (max_y_clamped - yoff + y_radius).min(src_raster_size[1] - yoff),
        );
        (yoff, ysize)
    };

    // Fill ratio (GDAL lines 3360–3364)
    let fill_ratio = (xsize as f64 * ysize as f64)
        / ((max_x - min_x + 2.0 * x_radius as f64)
            * (max_y - min_y + 2.0 * y_radius as f64))
            .max(1.0);

    Some(SourceWindow {
        xoff, yoff, xsize, ysize,
        fill_ratio, n_failed, n_samples,
    })
}

// =========================================================================
// Sampling strategies
// =========================================================================

/// Sample every pixel on every edge of the dest window.
/// Returns (min_x, min_y, max_x, max_y, n_samples, n_failed).
fn sample_edges_all(
    transformer: &impl Transformer,
    dst_off: [i32; 2],
    dst_size: [i32; 2],
) -> (f64, f64, f64, f64, i32, i32) {
    let n_max = 2 * (dst_size[0] + dst_size[1]) as usize;
    let mut xs = Vec::with_capacity(n_max);
    let mut ys = Vec::with_capacity(n_max);

    // Top and bottom edges
    for ix in 0..=dst_size[0] {
        xs.push((dst_off[0] + ix) as f64);
        ys.push(dst_off[1] as f64);
        xs.push((dst_off[0] + ix) as f64);
        ys.push((dst_off[1] + dst_size[1]) as f64);
    }
    // Left and right edges (excluding corners already done)
    for iy in 1..dst_size[1] {
        xs.push(dst_off[0] as f64);
        ys.push((dst_off[1] + iy) as f64);
        xs.push((dst_off[0] + dst_size[0]) as f64);
        ys.push((dst_off[1] + iy) as f64);
    }

    collect_bounds(transformer, &mut xs, &mut ys)
}

/// Sample step_count points along each edge of the dest window.
/// Returns (min_x, min_y, max_x, max_y, n_samples, n_failed).
fn sample_edges_step(
    transformer: &impl Transformer,
    dst_off: [i32; 2],
    dst_size: [i32; 2],
    step_count: usize,
) -> (f64, f64, f64, f64, i32, i32) {
    let n_max = step_count * 4;
    let mut xs = Vec::with_capacity(n_max);
    let mut ys = Vec::with_capacity(n_max);

    let step_size = 1.0 / (step_count - 1) as f64;
    let mut ratio = 0.0;
    while ratio <= 1.0 + step_size * 0.5 {
        let dx = ratio * dst_size[0] as f64;
        let dy = ratio * dst_size[1] as f64;

        // Top edge
        xs.push(dst_off[0] as f64 + dx);
        ys.push(dst_off[1] as f64);

        // Bottom edge
        xs.push(dst_off[0] as f64 + dx);
        ys.push((dst_off[1] + dst_size[1]) as f64);

        // Left edge
        xs.push(dst_off[0] as f64);
        ys.push(dst_off[1] as f64 + dy);

        // Right edge
        xs.push((dst_off[0] + dst_size[0]) as f64);
        ys.push(dst_off[1] as f64 + dy);

        ratio += step_size;
    }

    collect_bounds(transformer, &mut xs, &mut ys)
}

/// Sample on a (step_count+2)² grid through the dest window interior.
/// The +2 adds near-edge points at 0.5/nDstXSize offset (GDAL lines 2846–2863).
/// Returns (min_x, min_y, max_x, max_y, n_samples, n_failed).
fn sample_grid(
    transformer: &impl Transformer,
    dst_off: [i32; 2],
    dst_size: [i32; 2],
    step_count: usize,
) -> (f64, f64, f64, f64, i32, i32) {
    let n_grid = step_count + 2;
    let n_max = n_grid * n_grid;
    let mut xs = Vec::with_capacity(n_max);
    let mut ys = Vec::with_capacity(n_max);

    let step_size = 1.0 / (step_count - 1) as f64;
    let near_edge = 0.5 / dst_size[0].max(1) as f64;

    for iy in 0..n_grid {
        let ratio_y = if iy == 0 {
            near_edge
        } else if iy <= step_count {
            (iy - 1) as f64 * step_size
        } else {
            1.0 - near_edge
        };

        for ix in 0..n_grid {
            let ratio_x = if ix == 0 {
                near_edge
            } else if ix <= step_count {
                (ix - 1) as f64 * step_size
            } else {
                1.0 - near_edge
            };

            xs.push(ratio_x * dst_size[0] as f64 + dst_off[0] as f64);
            ys.push(ratio_y * dst_size[1] as f64 + dst_off[1] as f64);
        }
    }

    collect_bounds(transformer, &mut xs, &mut ys)
}

/// Transform points and collect bounding box of valid results.
fn collect_bounds(
    transformer: &impl Transformer,
    xs: &mut Vec<f64>,
    ys: &mut Vec<f64>,
) -> (f64, f64, f64, f64, i32, i32) {
    let n = xs.len();
    let ok = transformer.transform(true, xs, ys);

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut n_failed = 0i32;

    for i in 0..n {
        if !ok[i] || xs[i].is_nan() || ys[i].is_nan() {
            n_failed += 1;
            continue;
        }
        if xs[i] < min_x { min_x = xs[i]; }
        if xs[i] > max_x { max_x = xs[i]; }
        if ys[i] < min_y { min_y = ys[i]; }
        if ys[i] > max_y { max_y = ys[i]; }
    }

    (min_x, min_y, max_x, max_y, n as i32, n_failed)
}
