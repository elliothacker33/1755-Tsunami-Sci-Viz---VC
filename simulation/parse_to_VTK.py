# SPDX-License-Identifier: MIT
# Copyright (c) 2026
# Diogo Silva, Frederico Silva, Tomás Pereira

# =============  Parse GeoClaw AMR output to VTK format  ================= #
# ============  Supports multiple AMR grids per fort.qXXXX  ============= #
# ====================================================================== #

import numpy as np
import os
from pathlib import Path
from pyevtk.hl import gridToVTK

# ----------------------------------------------------------------------
# Main conversion routine
# ----------------------------------------------------------------------
def geoclaw_to_vtk(output_dir, vtk_dir):

    os.makedirs(vtk_dir, exist_ok=True)
    output_path = Path(output_dir)

    # Procura todos os ficheiros fort.q* ordenados
    fort_files = sorted(output_path.glob("fort.q*"))
    if not fort_files:
        print("No fort.q files found.")
        return

    for step, fort_file in enumerate(fort_files):

        # Lê as grades AMR do ficheiro fort.qXXXX
        grids = read_geoclaw_amr(fort_file)

        for gidx, grid in enumerate(grids):

            x = grid["x"]  # coordenadas x dos pontos (tamanho mx+1)
            y = grid["y"]  # coordenadas y dos pontos (tamanho my+1)

            # Campos lidos para as células (mx, my)
            h = grid["h"]
            hu = grid["hu"]
            hv = grid["hv"]
            eta = grid["eta"]

            # ----------------------------------------------------------
            # Calcular velocidades u,v (evitar divisão por zero)
            # ----------------------------------------------------------
            u = np.zeros_like(h)
            v = np.zeros_like(h)

            wet = h > 1.0e-6
            u[wet] = hu[wet] / h[wet]
            v[wet] = hv[wet] / h[wet]

            vel_mag = np.sqrt(u**2 + v**2)

            # ----------------------------------------------------------
            # Dimensões da malha (pontos)
            # Para VTK espera-se (nx, ny, nz) pontos
            # Temos dados 2D, então z será dimensão unitária 1
            # ----------------------------------------------------------
            nx = x.size
            ny = y.size
            nz = 1  # 2D plano

            # Criar array z unidimensional (ex: 0)
            z = np.array([0.0])

            # ----------------------------------------------------------
            # Converter dados centrados em células para dados em pontos
            # Média simples para interpolar dados das células para os pontos
            # Para um campo 2D mx x my (células), cria um array (mx+1, my+1) para pontos
            # ----------------------------------------------------------
            def cell_to_point(data):
                # data shape (mx, my)
                mx, my = data.shape
                pdata = np.zeros((mx+1, my+1))
                # média dos valores das 4 células vizinhas para cada ponto interior
                pdata[1:-1,1:-1] = (data[:-1,:-1] + data[1:,:-1] + data[:-1,1:] + data[1:,1:]) / 4.0
                # Bordas (copiar valores da célula adjacente)
                pdata[0,1:-1] = data[0,:-1]
                pdata[-1,1:-1] = data[-1,:-1]
                pdata[1:-1,0] = data[:-1,0]
                pdata[1:-1,-1] = data[:-1,-1]
                # Cantos
                pdata[0,0] = data[0,0]
                pdata[0,-1] = data[0,-1]
                pdata[-1,0] = data[-1,0]
                pdata[-1,-1] = data[-1,-1]
                return pdata

            # Convertendo para dados em pontos
            h_points = cell_to_point(h)
            eta_points = cell_to_point(eta)
            u_points = cell_to_point(u)
            v_points = cell_to_point(v)
            vel_mag_points = cell_to_point(vel_mag)

            # ----------------------------------------------------------
            # Nome do ficheiro VTK de saída
            # ----------------------------------------------------------
            vtk_name = (
                f"tsunami_t{step:04d}"
                f"_grid{grid['grid_number']:03d}"
                f"_lvl{grid['level']}"
            )

            vtk_path = os.path.join(vtk_dir, vtk_name)

            # ----------------------------------------------------------
            # Escreve ficheiro VTK no formato RectilinearGrid
            # Coloca as coordenadas e os dados como pointData
            # ----------------------------------------------------------
            gridToVTK(
                vtk_path,
                x, y, z,
                pointData={
                    "water_depth": h_points[..., np.newaxis],
                    "eta": eta_points[..., np.newaxis],
                    "velocity_x": u_points[..., np.newaxis],
                    "velocity_y": v_points[..., np.newaxis],
                    "velocity_magnitude": vel_mag_points[..., np.newaxis],
                }
            )


# ----------------------------------------------------------------------
# GeoClaw AMR parser
# ----------------------------------------------------------------------
def read_geoclaw_amr(filename):

    grids = []

    with open(filename, "r") as f:
        lines = f.readlines()

    idx = 0
    nlines = len(lines)

    while idx < nlines:

        # ----------------------------------------------------------
        # Skip empty lines
        # ----------------------------------------------------------
        if not lines[idx].strip():
            idx += 1
            continue

        # ----------------------------------------------------------
        # Read grid header
        # ----------------------------------------------------------
        try:
            grid_number = int(lines[idx].split()[0]); idx += 1
            level       = int(lines[idx].split()[0]); idx += 1
            mx          = int(lines[idx].split()[0]); idx += 1
            my          = int(lines[idx].split()[0]); idx += 1
            xlow        = float(lines[idx].split()[0]); idx += 1
            ylow        = float(lines[idx].split()[0]); idx += 1
            dx          = float(lines[idx].split()[0]); idx += 1
            dy          = float(lines[idx].split()[0]); idx += 1
        except Exception:
            # If header parsing fails, abort cleanly
            break

        # ----------------------------------------------------------
        # Build grid coordinates for points (mx+1, my+1)
        # ----------------------------------------------------------
        x = np.linspace(xlow, xlow + mx * dx, mx + 1)
        y = np.linspace(ylow, ylow + my * dy, my + 1)

        # ----------------------------------------------------------
        # Allocate cell-centered arrays
        # Shape: (mx, my) to match data read from file
        # ----------------------------------------------------------
        h   = np.zeros((mx, my))
        hu  = np.zeros((mx, my))
        hv  = np.zeros((mx, my))
        eta = np.zeros((mx, my))

        # ----------------------------------------------------------
        # Read mx * my cell values
        # ----------------------------------------------------------
        for j in range(my):
            for i in range(mx):
                while not lines[idx].strip():
                    idx += 1

                vals = lines[idx].split()
                h[i, j]   = float(vals[0])
                hu[i, j]  = float(vals[1])
                hv[i, j]  = float(vals[2])
                eta[i, j] = float(vals[3])

                idx += 1

        grids.append({
            "grid_number": grid_number,
            "level": level,
            "mx": mx,
            "my": my,
            "x": x,
            "y": y,
            "h": h,
            "hu": hu,
            "hv": hv,
            "eta": eta,
        })

    return grids


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("Converting GeoClaw AMR output to VTK...\n")
    geoclaw_to_vtk("_output", "_vtk")
    print("Done.")

