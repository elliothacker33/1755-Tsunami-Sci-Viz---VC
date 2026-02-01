# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Diogo Silva, Frederico Afonso, Tomás Pereira

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from pathlib import Path
from clawpack.pyclaw import Solution
from scipy.interpolate import RegularGridInterpolator
import matplotlib.ticker as ticker

# ====================================== #
#  Elsevier plotting configurations      #
# ====================================== #

width_inch = 190 / 25.4
height_inch = width_inch * 0.45

plt.rcParams.update({
     "figure.figsize": (width_inch, height_inch),
     "font.size": 10,
     "axes.labelsize": 11,
     "xtick.labelsize": 9,
     "ytick.labelsize": 9,
     "legend.fontsize": 8,
     "savefig.dpi": 600,
     "figure.autolayout": True,
     "lines.linewidth": 1.2
 })

gauge_names = {
    "00001": "Lisbon (Terreiro do Paço)",
    "00002": "Lisbon (Cascais Harbor)",
    "00003": "Tagus River (Cacilhas)",
    "00004": "Offshore",
    "00005": "Fault (Horsehoe - HSF)",
    "00006": "Lisbon (Forte do Bugio)"
}

# Style copied from paper in report
def marigrams_gauges(output_dir, plots_dir):

    output_path = Path(output_dir)
    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)

    gauge_files = sorted(output_path.glob("gauge*.txt"))

    amr_styles = {
        1: {'color': '#bdc3c7', 'lw': 0.7, 'alpha': 0.6, 'label': 'L1'},
        2: {'color': '#7f8c8d', 'lw': 0.8, 'alpha': 0.8, 'label': 'L2'},
        3: {'color': '#2c3e50', 'lw': 1.0, 'alpha': 0.9, 'label': 'L3'},
        4: {'color': '#000000', 'lw': 1.2, 'alpha': 1.0, 'label': 'L4'}
    }

    for gauge_file in gauge_files:

        gid = gauge_file.stem.replace('gauge', '')
        gname = gauge_names.get(gid, f"Gauge {gid}")

        try:
            data = np.loadtxt(gauge_file)
            if data.size == 0: continue
            levels, times, eta = data[:, 0].astype(int), data[:, 1] / 60.0, data[:, 5]
        except: continue

        plt.style.use(['science', 'no-latex'])
        fig, ax = plt.subplots(figsize=(190/25.4, 190/25.4 * 0.3))

        ax.axhline(0, color='black', lw=0.5)

        for level, style in amr_styles.items():
            mask = (levels == level)
            if np.any(mask):
                ax.plot(np.where(mask, times, np.nan),
                        np.where(mask, eta, np.nan),
                        color=style['color'],
                        linewidth=style['lw'],
                        alpha=style['alpha'],
                        zorder=level)

        ax.text(0.98, 0.85, f"{gid} - {gname}",
                transform=ax.transAxes,
                fontsize=11,
                fontweight='bold',
                fontstyle='italic',
                verticalalignment='top',
                horizontalalignment='right')

        ax.set_xlabel('Time (min)', fontsize=10)
        ax.set_ylabel('$\zeta$ (m)', fontsize=10)
        ax.set_xlim(0, times.max())
        ax.tick_params(direction='in', top=True, right=True)
        ax.grid(False)

        plt.savefig(plots_path / f"marigram_{gid}.pdf", bbox_inches='tight', dpi=600)
        plt.close()

def maximum_wave_height(output_dir, plots_dir):

    output_path = Path(output_dir)
    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)

    fort_files = sorted([f for f in output_path.glob("fort.q*") if any(c.isdigit() for c in f.name)])

    if not fort_files:
        print("No fort.q files found.")
        return

    lon_min, lon_max = -12.0, -6.0
    lat_min, lat_max = 36.0, 40.0

    resolution = 1000

    lon_fixed = np.linspace(lon_min, lon_max, resolution)
    lat_fixed = np.linspace(lat_min, lat_max, resolution)
    LON, LAT = np.meshgrid(lon_fixed, lat_fixed)

    max_eta_global = np.zeros_like(LON)
    tolerance = 1e-3

    for fort_file in fort_files:

        frame_num_str = ''.join(filter(str.isdigit, fort_file.name))
        if not frame_num_str:
            continue
        frame_num = int(frame_num_str)

        try:
            sol = Solution(frame_num, path=str(output_path), file_format='ascii')
        except Exception as e:
            print(f"Skipping frame {frame_num}: {e}")
            continue

        for state in sol.states:
            h = state.q[0, :, :]

            try:
                eta = state.q[3, :, :]
            except IndexError:
                print("Error in eta Index: {e}")
                continue

            eta_masked = np.where(h > tolerance, eta, -9999).T

            x_patch = state.grid.x.centers
            y_patch = state.grid.y.centers

            # Just fail-safe for other simulations (NOT NECCESSARY ON THIS SETTINGS)
            if (x_patch.max() < lon_min or x_patch.min() > lon_max or
                y_patch.max() < lat_min or y_patch.min() > lat_max):
                continue

            try:
                interp_func = RegularGridInterpolator(
                    (y_patch, x_patch),
                    eta_masked,
                    bounds_error=False,
                    fill_value=-9999
                )
            except ValueError as e:
                print(f"Interpolation error in frame {frame_num}: {e}")
                continue

            pts = np.array([LAT.ravel(), LON.ravel()]).T
            interpolated_values = interp_func(pts).reshape(LON.shape)

            valid_mask = interpolated_values > -9000
            max_eta_global[valid_mask] = np.maximum(
                max_eta_global[valid_mask],
                interpolated_values[valid_mask]
            )

    plt.style.use(['science', 'no-latex'])
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))

    from matplotlib.colors import PowerNorm
    cmap = plt.get_cmap('turbo')
    cmap.set_under('white', alpha=0)

    norm = PowerNorm(gamma=0.6, vmin=0.1, vmax=7.0)

    im = ax.pcolormesh(LON, LAT, max_eta_global,
                       cmap=cmap,
                       norm=norm,
                       shading='auto')

    cbar = fig.colorbar(im, ax=ax, extend='max', aspect=20, pad=0.02)
    cbar.set_label('$\zeta_{max}$ (m)')

    tick_values = [0, 1, 2, 3, 4, 5, 6, 7]
    cbar.set_ticks(tick_values)
    cbar.outline.set_linewidth(0.8)
    for val in tick_values:
        cbar.ax.axhline(val, color='black', linewidth=0.8)

    ax.set_xlabel('Longitude ($^\circ$)')
    ax.set_ylabel('Latitude ($^\circ$)')
    ax.set_aspect('equal')

    plt.savefig(plots_path / "maximum_wave_height_map.pdf", bbox_inches='tight')
    plt.close()

def velocity_heatmap(): pass

if __name__ == "__main__":
    print("Starting to plot results for 1755 Lisbon Tsunami")
    output_dir = "_output"
    plots_dir = "../plots"

    marigrams_gauges(output_dir, f"{plots_dir}/marigrams/")
    maximum_wave_height(output_dir,f"{plots_dir}/maximum_wave_height/")
