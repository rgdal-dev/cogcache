#' Decode a DEFLATE-compressed COG tile (Rust backend)
#'
#' Takes raw compressed bytes, decompresses with zlib, undoes predictor=2
#' horizontal differencing, and returns UInt16 pixel values as an integer
#' vector. This is the Rust equivalent of the R pipeline:
#' `memDecompress() |> readBin() |> cumsum per row`.
#'
#' @param raw_bytes Raw vector of compressed tile bytes
#' @param tile_width Integer tile width in pixels
#' @param tile_height Integer tile height in pixels
#' @param predictor Integer predictor type (1 = none, 2 = horizontal diff)
#' @return Integer vector of decoded pixel values (row-major order)
#' @export
rust_decode_tile <- function(raw_bytes, tile_width, tile_height, predictor = 2L) {
  .Call("wrap__rust_decode_tile", raw_bytes, tile_width, tile_height, predictor)
}

#' Fetch and decode a COG tile via HTTP range request (Rust backend)
#'
#' Performs the full pipeline in Rust: HTTP range fetch, DEFLATE decompress,
#' predictor undo. Returns decoded pixel values.
#'
#' @param url Character string URL (without /vsicurl/ prefix)
#' @param byte_offset Numeric byte offset in the file
#' @param byte_length Integer number of bytes to fetch
#' @param tile_width Integer tile width in pixels
#' @param tile_height Integer tile height in pixels
#' @param predictor Integer predictor type (1 = none, 2 = horizontal diff)
#' @return Integer vector of decoded pixel values (row-major order)
#' @export
rust_fetch_decode_tile <- function(url, byte_offset, byte_length,
                                    tile_width, tile_height, predictor = 2L) {
  .Call("wrap__rust_fetch_decode_tile", url, byte_offset, byte_length,
        tile_width, tile_height, predictor)
}
