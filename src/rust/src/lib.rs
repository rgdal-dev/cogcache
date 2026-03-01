use extendr_api::prelude::*;
use flate2::read::ZlibDecoder;
use std::io::Read;

/// Core decode logic — pure Rust, no extendr types.
/// Takes compressed bytes, tile dimensions, predictor.
/// Returns decoded UInt16 pixel values as Vec<i32>.
fn decode_deflate_pixels(
    compressed: &[u8],
    tile_width: usize,
    tile_height: usize,
    predictor: i32,
) -> extendr_api::Result<Robj> {
    let n_pixels = tile_width * tile_height;

    // Decompress with zlib (TIFF DEFLATE uses zlib-wrapped deflate)
    let mut decoder = ZlibDecoder::new(compressed);
    let mut decompressed = Vec::with_capacity(n_pixels * 2);
    decoder
        .read_to_end(&mut decompressed)
        .map_err(|e| Error::Other(format!("DEFLATE decompression failed: {}", e)))?;
       // .expect("DEFLATE decompression failed");

    assert_eq!(
        decompressed.len(),
        n_pixels * 2,
        "Decompressed size {} != expected {} ({}x{}x2)",
        decompressed.len(),
        n_pixels * 2,
        tile_width,
        tile_height
    );

    // Read as little-endian UInt16 into i32 (R integer range)
    let mut pixels: Vec<i32> = decompressed
    .chunks_exact(2)
    .map(|b| u16::from_le_bytes([b[0], b[1]]) as i32)
    .collect();

    // Undo predictor=2: cumsum each row, mod 65536
    if predictor == 2 {
        for row in 0..tile_height {
            let row_start = row * tile_width;
            for col in 1..tile_width {
                let idx = row_start + col;
                pixels[idx] = (pixels[idx] + pixels[idx - 1]) % 65536;
            }
        }
    }

    Ok(r!(pixels))
}

/// Decode a DEFLATE-compressed tile with predictor=2 undo.
///
/// Takes raw compressed bytes (as R raw vector), tile dimensions, and
/// predictor type. Returns decoded pixel values as integer vector.
///
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
) -> Robj {
    let pixels = decode_deflate_pixels(
        raw_bytes.as_slice(),
        tile_width as usize,
        tile_height as usize,
        predictor,
    );
    r!(pixels)
}

/// Fetch a tile via HTTP range request and decode it.
///
/// Performs the full pipeline in Rust: HTTP range fetch, DEFLATE decompress,
/// predictor undo. Returns decoded pixel values as integer vector.
///
/// @param url Character string URL (without /vsicurl/ prefix)
/// @param byte_offset Numeric byte offset in the file
/// @param byte_length Integer number of bytes to fetch
/// @param tile_width Integer tile width in pixels
/// @param tile_height Integer tile height in pixels
/// @param predictor Integer predictor type (1 = none, 2 = horizontal diff)
/// @return Integer vector of decoded UInt16 pixel values
/// @export
#[extendr]
fn rust_fetch_decode_tile(
    url: &str,
    byte_offset: f64,
    byte_length: i32,
    tile_width: i32,
    tile_height: i32,
    predictor: i32,
) -> extendr_api::Result<Robj>  {
    //this works because R doesn't have 64-bit integers and
    // the offsets fit in f64's 53-bit mantissa
    let offset = byte_offset as u64;
    let length = byte_length as u64;
    let range_end = offset + length - 1;
    let range_header = format!("bytes={}-{}", offset, range_end);

    // HTTP range request
    let resp = ureq::get(url)
        .set("Range", &range_header)
        .call()
        .map_err(|e| Error::Other(format!("HTTP request failed: {}", e)))?;
//        .expect("HTTP request failed");

    let mut compressed = Vec::with_capacity(byte_length as usize);
    resp.into_reader()
        .read_to_end(&mut compressed)
        .map_err(|e| Error::Other(format!("Failed to read response body: {}", e)))?;
        //.expect("Failed to read response body");

    let pixels = decode_deflate_pixels(
        &compressed,
        tile_width as usize,
        tile_height as usize,
        predictor,
    );
   Ok(r!(pixels))
}

extendr_module! {
    mod cogcache;
    fn rust_decode_tile;
    fn rust_fetch_decode_tile;
}
