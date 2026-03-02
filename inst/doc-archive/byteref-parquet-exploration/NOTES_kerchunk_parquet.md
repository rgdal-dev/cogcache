## kerchunk parquet from rustycogs — notes

### You don't need VirtualiZarr for this

The `.zmetadata` is just JSON describing the zarr array structure. For a
single-variable COG collection where every file has identical structure
(same dimensions, tile size, dtype, compression), everything needed is in
the gdalinfo output. The R script constructs it directly.

VirtualiZarr's value is in handling heterogeneous collections and complex
multi-variable NetCDF/HDF5 structures where the zarr metadata isn't obvious.
For a uniform COG timeseries, it's overhead — and the v2/v3 codec mismatch
makes it actively painful.

### Actual COG specs (from gdalinfo)

```
36000 x 18000 pixels, 512x512 tiles
COMPRESSION=ZSTD, PREDICTOR=2
Band 1 Type=Int16, NoData=-32768
Offset: 25, Scale: 0.001, Unit: celsius
Origin: (-179.995, 89.995), Pixel: (0.01, -0.01)
Overviews: 18000x9000, 9000x4500, 4500x2250, 2250x1125, 1125x563
```

Tile grid: ceil(36000/512) × ceil(18000/512) = 71 × 36 = **2556 tiles/file**

### The compressor question

TIFF tile bytes are encoded as: raw pixels → predictor=2 → ZSTD compress.

Three possible zarr v2 translations (script supports all as modes A/B/C):

| Mode | compressor | filters | Hypothesis |
|------|-----------|---------|------------|
| A | `{"id":"zstd"}` | `[{"id":"delta","dtype":"<i2"}]` | Full translation: zarr decompresses ZSTD, then delta filter reverses predictor |
| B | `{"id":"zstd"}` | `null` | GDAL handles predictor internally for kerchunk TIFF refs |
| C | `null` | `null` | GDAL handles entire decompression chain |

**Start with B.** GDAL's kerchunk implementation is TIFF-aware (Even Rouault
specifically discussed predictor handling in the gdal-dev thread). Mode A
is the "pure zarr" answer but `delta` may not match TIFF predictor=2 exactly
(delta is element-wise subtraction, TIFF predictor=2 is row-wise horizontal
differencing on the raw byte stream including interleaved bytes of multi-byte
types). Mode C is the fallback.

### Testing workflow

```r
source("assemble_kerchunk_parquet.R")

## 1. single file, 2D — fastest iteration
test_one_file(compressor_mode = "B")

## 2. check GDAL can open it
gdalraster::gdalmdiminfo("test_one.parq")

## 3. compare actual pixels
url <- "/vsicurl/https://data.source.coop/ausantarctic/ghrsst-mur-v2/2002/06/01/20020601090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1_analysed_sst.tif"
direct <- vapour::vapour_read_raster(url, window = c(0, 0, 512, 512))
zarr <- vapour::vapour_read_raster(
  'ZARR:"test_one.parq":/analysed_sst',
  window = c(0, 0, 512, 512))
identical(direct[[1]], zarr[[1]])

## 4. if wrong, iterate
test_one_file(compressor_mode = "A")  # try with delta filter
test_one_file(compressor_mode = "C")  # try null
```

### Full catalogue numbers

- ~8400 days (2002-06 to 2025-02)
- 2556 tiles/file
- **~21.5M reference rows** total
- Single parquet: ~200-400 MB compressed (snappy)
- Generation time: dominated by the HTTP IFD reads, not the parquet write
  (rustycogs is fast, the write is one arrow::write_parquet call)

### Partition size

Default `record_size = nrow(refs)` → one file. The kerchunk spec's default
of 10,000 would create 2,150 partition files (slower to write, no read
benefit for GDAL which loads the whole partition on first access anyway).

If GDAL complains about a single huge partition, fall back to
`record_size = 1000000` (~22 files).

### Overlap with hypertidy

- **rustycogs**: IFD scanning → reference table
- **arrow**: parquet I/O
- **grout**: tile index ↔ pixel coordinate mapping
- **vaster**: raster grid logic (extent, dimension, resolution)
- **vapour/gdalraster**: pixel data access (consumes the virtual store)
- **ndr**: eventual lazy array consumer

The parquet reference table is the interchange format. Each package
does one thing. This is just tables and byte ranges — no special
spatial types needed.
