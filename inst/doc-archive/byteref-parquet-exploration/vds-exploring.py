import xarray as xr
from obstore.store import from_url, LocalStore
from obspec_utils.registry import ObjectStoreRegistry

from virtualizarr import open_virtual_dataset, open_virtual_mfdataset
from virtualizarr.parsers import HDFParser

bucket = "https://www.ncei.noaa.gov"
path = "data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/198109/oisst-avhrr-v02r01.19810901.nc"
store = from_url(bucket)
bucket  = 'https://esgf-data.ucar.edu'
path = 'thredds/fileServer/esg_dataroot/CMIP6/CMIP/NCAR/CESM2/historical/r3i1p1f1/day/tas/gn/v20190308/tas_day_CESM2_historical_r3i1p1f1_gn_19200101-19291231.nc'
store = from_url(bucket)

url = "file:///rdsi/PUBLIC/raad/data/podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2002/152/20020601090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc"

registry = ObjectStoreRegistry()
registry.register(url, LocalStore())

vds = open_virtual_dataset(url=url, parser=HDFParser(), registry=registry)

m = vds["analysed_sst"].data.manifest
m.shape_chunk_grid
m._offsets.shape
