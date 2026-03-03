## Changes to lib.rs for resampling support

### 1. Update rust_warp_approx to use warp_resample internally

In `rust_warp_approx`, change:

```rust
    let output = warp::warp_nearest(
```

to:

```rust
    let output = warp::warp_resample(
```

and add `warp::ResampleAlg::NearestNeighbour,` as the last argument before the closing `)`.

### 2. Add rust_warp_resample after rust_warp_approx

Add this new function and register it in `extendr_module!`:

```rust
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
```

### 3. Register in extendr_module!

```rust
extendr_module! {
    mod cogcache;
    fn rust_decode_tile;
    fn rust_fetch_decode_tile;
    fn rust_gen_img_proj_transform;
    fn rust_warp_scanline;
    fn rust_warp_approx;
    fn rust_warp_resample;       // ← NEW
    fn rust_compute_source_window;
    fn rust_warp_map;
    fn rust_apply_warp;
}
```

### 4. Source window padding

When calling `rust_compute_source_window`, pass the resampling radius:

```r
sw <- rust_compute_source_window(
  ...,
  resample_padding = alg_radius  # 0=near, 1=bilinear, 2=cubic, 3=lanczos
)
```

This is already wired — the source_window adds the padding to the read window.
