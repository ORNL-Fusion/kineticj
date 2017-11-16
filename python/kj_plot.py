from netCDF4 import Dataset

dataset = Dataset('output/jP2.nc')

print dataset.file_format
