# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Diogo Silva, Frederico Silva, Tomás Pereira

# ============ Parse the NETCDF4 FILE to GeoClaw Format ================== #
# ==== Authors: Diogo Silva, Frederico Silva, Tomás Pereira ============== #
# ======================================================================== #

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LightSource

import scienceplots
import cmocean

# ====================================== #
#  Elsevier plotting configurations      #
# ====================================== #

width_inch = 190 / 25.4
height_inch = width_inch * 0.5

plt.rcParams.update({
    "figure.figsize": (width_inch, height_inch),
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 10,
    "savefig.dpi": 600,
    "figure.autolayout": True
})

class Coordinates:
    def __init__(self, lat_min, lat_max, lon_min, lon_max):
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lon_min = lon_min
        self.lon_max = lon_max

def amr_ascii_convert(input_file, output_file, plot_name, coordinates):

    root_cdf = Dataset(input_file, "r", format="NETCDF4")
    print(root_cdf)

    lons = root_cdf.variables['lon'][:]
    lats = root_cdf.variables['lat'][:]
    elevation = root_cdf.variables['elevation'][:]

    lon_idx = np.where((lons >= coordinates.lon_min) & (lons <= coordinates.lon_max))[0]
    lat_idx = np.where((lats >= coordinates.lat_min) & (lats <= coordinates.lat_max))[0]

    if lon_idx.size == 0 or lat_idx.size == 0:
        print("No points found in coordinate range.")
        root_cdf.close()
        return

    elv_region = elevation[lat_idx[0]:lat_idx[-1]+1, lon_idx[0]:lon_idx[-1]+1]
    lons_region = lons[lon_idx[0]:lon_idx[-1]+1]
    lats_region = lats[lat_idx[0]:lat_idx[-1]+1]

    nlon = len(lons_region)
    nlat = len(lats_region)

    with open(output_file, "w") as f:
        f.write("# NASA Aerogeophysics ASCII File Format Convention\n")
        f.write("# Dataset: NetCDF4 elevation subset\n")
        f.write("# Authors: Diogo Silva, Frederico Silva, Tomás Pereira\n")
        f.write(f"# Longitude range: {coordinates.lon_min} to {coordinates.lon_max}\n")
        f.write(f"# Latitude range: {coordinates.lat_min} to {coordinates.lat_max}\n")
        f.write(f"# Columns: {nlon}, Rows: {nlat}\n")
        f.write("# Data: Bathymetry / Elevation (meters)\n")
        f.write("# End of header\n")

        for i in range(nlat-1, -1, -1):
            f.write(" ".join(f"{elv_region[i, j]:.2f}" for j in range(nlon)) + "\n")

    plotting_map_bathymetry(lons_region, lats_region, elv_region, plot_name)

    root_cdf.close()

def plotting_map_bathymetry(lons, lats, elevation, plot_name):

    plt.style.use(['science', 'no-latex'])

    lon2d, lat2d = np.meshgrid(lons, lats)

    ls = LightSource(azdeg=315, altdeg=45)
    rgb = ls.shade(
        elevation,
        cmap=cmocean.cm.deep,
        vert_exag=0.5,
        blend_mode='soft'
    )

    fig = plt.figure(figsize=(width_inch, height_inch))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent(
        [lons.min(), lons.max(), lats.min(), lats.max()],
        crs=ccrs.PlateCarree()
    )

    ax.imshow(
        rgb,
        extent=[lons.min(), lons.max(), lats.min(), lats.max()],
        origin='lower',
        transform=ccrs.PlateCarree(),
        zorder=1
    )

    levels = np.concatenate([
        np.arange(-6000, -1000, 500),
        np.arange(-1000, -100, 100),
        np.arange(-100, 1, 20)
    ])

    cf = ax.contourf(
        lon2d, lat2d, elevation,
        levels=levels,
        cmap=cmocean.cm.deep,
        alpha=0.75,
        transform=ccrs.PlateCarree(),
        extend='both',
        zorder=2
    )

    isobath_levels = [-5000, -4000, -3000, -2000, -1000, -500, -200, -100]
    cs = ax.contour(
        lon2d, lat2d, elevation,
        levels=isobath_levels,
        colors='black',
        linewidths=0.5,
        alpha=0.6,
        transform=ccrs.PlateCarree(),
        zorder=3
    )

    ax.clabel(cs, inline=True, fmt='%d', fontsize=6.5, inline_spacing=3)

    ax.add_feature(cfeature.LAND, facecolor='#e8e8e8', edgecolor='none', zorder=10)
    ax.coastlines(resolution='10m', linewidth=0.7, color='black', zorder=11)

    ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor='gray', alpha=0.5, zorder=12)

    gl = ax.gridlines(
        draw_labels=True,
        linestyle=':',
        linewidth=0.5,
        alpha=0.6,
        color='gray',
        zorder=13
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

    cbar = plt.colorbar(
        cf,
        ax=ax,
        pad=0.02,
        shrink=0.9,
        aspect=20,
        format='%d'
    )
    cbar.set_label("Elevation (m)", fontsize=10)
    cbar.ax.tick_params(labelsize=7)

    plt.savefig(f"../plots/bathymetry/{plot_name}", bbox_inches="tight", dpi=600)
    plt.close()

gebco_path = "../data/GEBCO_data/"

coords_coarse = Coordinates(36.0, 40.0, -12.0, -6.0)
coords_medium = Coordinates(37.0, 39.0, -9.5, -8.5)
coords_fine   = Coordinates(38.6, 38.8, -9.3, -9.0)

amr_ascii_convert(
    gebco_path + "gebco_coarse_data.nc",
    gebco_path + "gebco_coarse_data.asc",
    "bathymetry_1755_coarse.pdf",
    coords_coarse
)

amr_ascii_convert(
    gebco_path + "gebco_medium_data.nc",
    gebco_path + "gebco_medium_data.asc",
    "bathymetry_1755_medium.pdf",
    coords_medium
)

amr_ascii_convert(
    gebco_path + "gebco_fine_data.nc",
    gebco_path + "gebco_fine_data.asc",
    "bathymetry_1755_fine.pdf",
    coords_fine
)
