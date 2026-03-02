# Sys.setenv(LD_LIBRARY_PATH = "")
# reticulate::use_python("~/.venv/bin/python3")
# library(reticulate)
# repl_python()

# library(mirai)
# daemons(4)
# files <- sds::ghrsst(vsi = FALSE)
# library(rustycogs)
# library(purrr)
# fun <- in_parallel(function(.x) {
#   rustycogs::tiff_refs(.x, "", TRUE, 180)
# })
# 
# refs <- map_dfr(files$source, fun)
# arrow::write_parquet(refs, "~/wholecat/ghrsst-refs.parquet")

import os
os.environ["LD_LIBRARY_PATH"] = ""

import numpy as np
import xarray as xr
import pyarrow.parquet as pq
from virtualizarr.manifests import ChunkManifest, ManifestArray

# ---- Step 1: read refs from wherever ----
#table = pq.read_table("/perm_storage/home/mdsumner/wholecat/ghrsst-refs.parquet")

table = pq.read_table(
     "/perm_storage/home/mdsumner/wholecat/ghrsst-refs.parquet",
     filters=[("ifd", "==", 0)]
 )

# columns: path (string), offset (uint64), length (uint64)
# either with chunk_t/chunk_y/chunk_x columns, or in C-contiguous order
import pyarrow.compute as pc

unique_paths = table["path"].unique()   # pa.Array
# ---- Step 2: build shaped numpy arrays (the thing to benchmark) ----
n = len(table)

nn = len(unique_paths)



shape = (nn,
    #pc.max(table["chunk_t"]).as_py() + 1,
    pc.max(table["tile_row"]).as_py() + 1,
    pc.max(table["tile_col"]).as_py() + 1,
)

# Method A: just reshape (if C-contiguous and complete)
paths = table["path"].to_numpy(zero_copy_only=False).reshape(shape)
offsets = table["offset"].to_numpy().reshape(shape)
lengths = table["length"].to_numpy().reshape(shape)

# Method B: scatter via index columns
encoded = pc.dictionary_encode(table["path"]).combine_chunks()
chunk_t = encoded.indices.to_numpy()

idx = (
    chunk_t,
    table["tile_row"].to_numpy(),
    table["tile_col"].to_numpy(),
)
paths = np.empty(shape, dtype=np.dtypes.StringDType())
offsets = np.zeros(shape, dtype="uint64")
lengths = np.zeros(shape, dtype="uint64")
paths[idx] = table["path"].to_numpy(zero_copy_only=False)
offsets[idx] = table["offset"].to_numpy()
lengths[idx] = table["length"].to_numpy()

# ---- Step 3: construct the virtual dataset ----
manifest = ChunkManifest.from_arrays(
    paths=paths, offsets=offsets, lengths=lengths
)

# from zarr.codecs import BytesCodec
# from virtual_tiff.codecs import HorizontalDeltaCodec
# from virtual_tiff.imagecodecs import ZstdCodec
# 
# codecs = [
#     HorizontalDeltaCodec(),
#     BytesCodec(endian="little"),
#     ZstdCodec(level=9),
# ]

# Steal metadata from a reference vds, or build from scratch:
from virtualizarr.manifests.utils import create_v3_array_metadata
metadata = create_v3_array_metadata(
    shape=(nn, 17999, 36000),
    data_type=np.dtype("int16"),
    chunk_shape=(1, 1023, 2047),  # or whatever
    fill_value=-32768,
   # codecs= codecs,  # from vds, or from the enriched parquet zarray, or whatever
    dimension_names=("time", "lat", "lon"),
)

ma = ManifestArray(metadata=metadata, chunkmanifest=manifest)
print("about to create vds")
# Build the xarray Dataset
vds = xr.Dataset({"analysed_sst": xr.Variable(("time", "lat", "lon"), ma)})
print("about to pickle vds: vds.pkl")
import json, pickle
pickle.dump(vds, open("vds.pkl", "wb"))

# ---- Now you have a live virtual dataset ----
# You can save/restore the pieces:
#   - metadata is JSON-serializable (it's zarr v3 metadata)
#   - manifest is 3 numpy arrays you can save with np.save or as parquet
#   - vds itself can be written to kerchunk/icechunk
