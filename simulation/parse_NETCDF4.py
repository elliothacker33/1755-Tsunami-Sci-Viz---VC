# ============ Parse the NETCDF4 FILE to GeoClaw Format ================== #
# ==== Authors: Diogo Silva, Frederico Silva, Tomás Pereira ============== #
# ======================================================================== #

# File format link:
# https://www.earthdata.nasa.gov/about/edsis/esco/standards-practices/ascii-file-format-guidelines-earth-science-data/aerogeophysics #

import os
from netCDF4 import Dataset

# ============================================= #
# Coordinates class to process on netCDF4 files #
# ============================================= #

class Coordinates:

    def __init__(self, lat_min, lat_max, lon_min, lon_max):
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lon_min = lon_min
        self.lon_max = lon_max

# ===============================================#
# Convert the NetCDF4 to ASCII with AMR mesh     #
# ============================================== #
# The file generated will follow the rules of    #
# NASA Aerogeophysics ASCII File Format          #
# ============================================== #

def amr_ascii_convert(input_file, output_file, coordinates):

    # Get input file
    root_cdf = Dataset(input_file, "r", format="NETCDF4")
    print(root_cdf)

    lons = root_cdf.variables['lon'][:]
    lats = root_cdf.variables['lat'][:]
    elevation = root_cdf.variables['elevation'][:]

    lon_indices = [i for i, lon in enumerate(lons) if coordinates.lon_min <= lon <= coordinates.lon_max]
    lat_indices = [j for j, lat in enumerate(lats) if coordinates.lat_min <= lat <= coordinates.lat_max]

    if not lon_indices or not lat_indices:
        print("No points found in the specified coordinate range.")
        root_cdf.close()
        return

    elv_region = elevation[lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1]
    lons_region = lons[lon_indices[0]:lon_indices[-1]+1]
    lats_region = lats[lat_indices[0]:lat_indices[-1]+1]

    nlon = len(lons_region)
    nlat = len(lats_region)

    with open(output_file, "w") as f:

        # Header
        f.write(f"# NASA Aerogeophysics ASCII File Format Convention\n")
        f.write(f"# Dataset: NetCDF4 elevation subset\n")
        f.write(f"# Authors: Diogo Silva, Frederico Silva, Tomás Pereira\n")
        f.write(f"# Longitude range: {coordinates.lon_min} to {coordinates.lon_max}\n")
        f.write(f"# Latitude range: {coordinates.lat_min} to {coordinates.lat_max}\n")
        f.write(f"# Columns: {nlon}, Rows: {nlat}\n")
        f.write(f"# Data: Bathymetry / Elevation (meters)\n")
        f.write(f"# End of header\n")

        for i in range(nlat-1, -1, -1):
            row_data = ' '.join(f"{elv_region[i, j]:.2f}" for j in range(nlon))
            f.write(row_data + '\n')

    root_cdf.close()

# ============================================ #
# Tsunami in different regions                 #
# ============================================ #

# Using Adaptive Mesh Refinement for different regions
gebco_path = "../data/GEBCO_data/"
gebco_file = gebco_path + gebco_data.nc

coords_coarse = Coordinates(36.0, 40.0, -12.0, -6.0)
coords_medium = Coordinates(37.0, 39.0, -9.5, -8.5)
coords_fine = Coordinates(38.6, 38.8, -9.3, -9.0)

# 1. Coarse-grained
amr_ascii_convert(gebco_file, gebco_path + "gebco_coarse_data.asc", coords_coarse)

# 2. Medium-grained
amr_ascii_convert(gebco_file, gebco_path + "gebco_medium_data.asc", coords_medium)

# 3. Fine-grained
amr_ascii_convert(gebco_file, gebco_path + "gebco_fine_data.asc", coords_fine)

