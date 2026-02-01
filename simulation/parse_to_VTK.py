# SPDX-License-Identifier: MIT
# Copyright (c) 2026
# Diogo Silva, Frederico Afonso, Tomás Pereira

# ========================================================================= #
# =============   Parse GeoClaw AMR output to VTK format   ================= #
# ========================================================================= #

import numpy as np
import os
import sys
from pathlib import Path
from pyevtk.hl import gridToVTK

def write_statistics(stats_dir, h, b, eta, vel_mag, u, v, dx, dy, step, sim_time):
    stat_path = os.path.join(stats_dir, f"stats_step_{step:04d}.txt")
    dry_tolerance = 1e-3
    sea_level = 0.0

    wet = h > dry_tolerance
    land = b > 0
    inundated = wet & land
    ocean_wet = wet & ~land
    max_eta = np.max(eta[wet]) if np.any(wet) else 0.0
    min_eta = np.min(eta[wet]) if np.any(wet) else 0.0
    max_crest = np.max(eta[ocean_wet] - sea_level) if np.any(ocean_wet) else 0.0

    mean_eta_pos = np.mean(eta[wet & (eta > sea_level)]) if np.any(wet & (eta > sea_level)) else 0.0

    max_inundation_depth = np.max(h[inundated]) if np.any(inundated) else 0.0
    max_runup_height = np.max(b[inundated]) if np.any(inundated) else 0.0

    area_per_cell = dx * dy
    total_area_inundated = np.sum(inundated) * area_per_cell

    max_vel = np.max(vel_mag[wet]) if np.any(wet) else 0.0

    momentum_flux = h * (vel_mag**2)
    max_flux = np.max(momentum_flux[wet]) if np.any(wet) else 0.0

    mean_u = np.mean(u[wet]) if np.any(wet) else 0.0
    mean_v = np.mean(v[wet]) if np.any(wet) else 0.0

    with open(stat_path, "w") as f:
        f.write("# ================================================= #\n")
        f.write(f"#           Statistics for 1755 - Step {step:04d}         #\n")
        f.write(f"#           Sim Time: {sim_time/60:10.2f} min             #\n")
        f.write("# ================================================= #\n\n")
        f.write(f"--- ELEVATION & WAVES ---\n")
        f.write(f"Absolute Max Eta (Wet):   {max_eta:10.4f} m\n")
        f.write(f"Max Wave Crest (Ocean):   {max_crest:10.4f} m\n")
        f.write(f"Max Wave Trough (Eta):    {min_eta:10.4f} m\n")
        f.write(f"Mean Positive Elevation:  {mean_eta_pos:10.4f} m\n\n")
        f.write(f"--- INUNDATION (Land) ---\n")
        f.write(f"Max Flow Depth (h):       {max_inundation_depth:10.4f} m\n")
        f.write(f"Max Run-up Elevation (B): {max_runup_height:10.4f} m\n")
        f.write(f"Total Area Inundated:     {total_area_inundated:10.4e} m^2\n\n")
        f.write(f"--- DYNAMICS ---\n")
        f.write(f"Max Velocity:             {max_vel:10.4f} m/s\n")
        f.write(f"Max Momentum Flux:        {max_flux:10.4f} m^3/s^2\n")
        f.write(f"Mean Direction (U, V):    {mean_u:.3f}, {mean_v:.3f}\n")

def write_pvd(vtk_dir, timesteps_dict):
    pvd_path = os.path.join(vtk_dir, "tsunami_1755.pvd")
    with open(pvd_path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1">\n')
        f.write('  <Collection>\n')
        for t in sorted(timesteps_dict.keys()):
            for vtr_file in timesteps_dict[t]:
                f.write(f'    <DataSet timestep="{t}" file="{vtr_file}"/>\n')
        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')

def geoclaw_to_vtk(stats_dir="_stats", output_dir="_output", vtk_dir="_vtk"):

    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(vtk_dir, exist_ok=True)
    output_path = Path(output_dir)

    fort_files = sorted(output_path.glob("fort.q*"))

    if not fort_files:
        return

    timesteps_collection = {}

    for step, fort_file in enumerate(fort_files):
        t_file = fort_file.with_name(fort_file.name.replace('.q', '.t'))
        sim_time = step * 300.0  # Default ~5 min

        if t_file.exists():
            try:
                with open(t_file, 'r') as f:
                    sim_time = float(f.readline().split()[0])
            except:
                pass

        print(f"[{step+1:3d}/{len(fort_files)}] t={sim_time/60:6.1f} min ... ", end='')

        grids = read_geoclaw_amr(fort_file)

        if not grids:
            print("Skipped")
            continue

        timesteps_collection[sim_time] = []

        step_h, step_b, step_eta, step_vel, step_u, step_v = [], [], [], [], [], []
        last_dx, last_dy = 0.0, 0.0

        for gidx, grid in enumerate(grids):
            x = grid["x"]
            y = grid["y"]
            z = np.array([0.0, 1.0])

            h = grid["h"]
            hu = grid["hu"]
            hv = grid["hv"]
            eta = grid["eta"]
            b = grid["b"]
            last_dx, last_dy = grid["dx"], grid["dy"]

            u = np.zeros_like(h)
            v = np.zeros_like(h)
            wet = h > 1.0e-3
            u[wet] = hu[wet] / h[wet]
            v[wet] = hv[wet] / h[wet]

            vel_mag = np.sqrt(u**2 + v**2)

            step_h.append(h.flatten())
            step_b.append(b.flatten())
            step_eta.append(eta.flatten())
            step_vel.append(vel_mag.flatten())
            step_u.append(u.flatten())
            step_v.append(v.flatten())

            h_3d = np.repeat(h[:, :, np.newaxis], 2, axis=2)
            eta_3d = np.repeat(eta[:, :, np.newaxis], 2, axis=2)
            b_3d = np.repeat(b[:, :, np.newaxis], 2, axis=2)
            u_3d = np.repeat(u[:, :, np.newaxis], 2, axis=2)
            v_3d = np.repeat(v[:, :, np.newaxis], 2, axis=2)
            vel_3d = np.repeat(vel_mag[:, :, np.newaxis], 2, axis=2)

            vtr_filename = f"tsunami_step{step:04d}_grid{grid['grid_number']:02d}_level{grid['level']}"
            vtr_path = os.path.join(vtk_dir, vtr_filename)

            gridToVTK(
                vtr_path,
                x, y, z,
                pointData={
                    "surface_elevation": eta_3d,
                    "water_depth": h_3d,
                    "bathymetry": b_3d,
                    "velocity_x": u_3d,
                    "velocity_y": v_3d,
                    "velocity_magnitude": vel_3d,
                }
            )

            timesteps_collection[sim_time].append(vtr_filename + ".vtr")

        write_statistics(stats_dir,
                         np.concatenate(step_h), np.concatenate(step_b),
                         np.concatenate(step_eta), np.concatenate(step_vel),
                         np.concatenate(step_u), np.concatenate(step_v),
                         last_dx, last_dy, step, sim_time)
        print("Done")

    write_pvd(vtk_dir, timesteps_collection)

def read_geoclaw_amr(filename):

    grids = []
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file")
        return grids

    idx = 0
    while idx < len(lines):
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            break

        try:
            grid_number = int(lines[idx].split()[0]); idx += 1
            level = int(lines[idx].split()[0]); idx += 1
            mx = int(lines[idx].split()[0]); idx += 1
            my = int(lines[idx].split()[0]); idx += 1
            xlow = float(lines[idx].split()[0]); idx += 1
            ylow = float(lines[idx].split()[0]); idx += 1
            dx = float(lines[idx].split()[0]); idx += 1
            dy = float(lines[idx].split()[0]); idx += 1

            while idx < len(lines) and not lines[idx].strip():
                idx += 1

            x = xlow + (np.arange(mx) + 0.5) * dx
            y = ylow + (np.arange(my) + 0.5) * dy

            h = np.zeros((mx, my))
            hu = np.zeros((mx, my))
            hv = np.zeros((mx, my))
            eta = np.zeros((mx, my))

            for j in range(my):
                for i in range(mx):
                    while idx < len(lines) and not lines[idx].strip():
                        idx += 1
                    if idx >= len(lines):
                        break
                    vals = lines[idx].split()
                    if len(vals) >= 4:
                        h[i, j] = float(vals[0])
                        hu[i, j] = float(vals[1])
                        hv[i, j] = float(vals[2])
                        eta[i, j] = float(vals[3])
                    idx += 1

            b = eta - h

            # Append patch
            grids.append({
                "grid_number": grid_number,
                "level": level,
                "x": x, "y": y,
                "h": h, "hu": hu, "hv": hv, "eta": eta, "b": b,
                "dx": dx, "dy": dy
            })

        except Exception as e:
            print(f"Error processing grid")
            break
    return grids

if __name__ == "__main__":
    stats_dir = sys.argv[1] if len(sys.argv) > 1 else "_stats"
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "_output"
    vtk_dir = sys.argv[3] if len(sys.argv) > 3 else "_vtk"

    geoclaw_to_vtk(stats_dir, output_dir, vtk_dir)
