# 1755 Lisbon Tsunami — Sci-Viz Reconstruction (GeoClaw + ParaView)
**with representative near-field flow analysis (OpenFOAM)**

This repository contains the code, configurations, and visualization assets used in our project:

**“A Sci-Viz Reconstruction and Analysis of the 1755 Lisbon Tsunami”**

The work is organized around two complementary components:

1. **Far-field tsunami modelling (GeoClaw)** — large-scale propagation from the Horseshoe Fault using AMR.
2. **Near-field hydrodynamics (OpenFOAM)** — representative, physically plausible estuarine flow scenarios (e.g., coherent structures and vortex shedding-like patterns) driven by tsunami-inspired inflow conditions.

> Note: the OpenFOAM part is **not** a direct reconstruction of 1755. It is a **representative** exploration of mechanisms that *could plausibly occur* under extreme forcing.

---

## License
MIT License — Copyright (c) 2026  
Diogo Silva, Frederico Silva, Tomás Pereira

---

## Repository structure

```
.
├── simulation/
│   ├── geoclaw/ # 1755 tsunami modelling (AMR, gauges, outputs)
│   ├── openfoam/
│   │   ├── case_S/ # 2D cylinder in channel (small preset)
│   │   ├── case_M/ # medium preset
│   │   ├── case_L/ # large preset
│   │   ├── case_S_2cyl/ # two cylinders variant
│   │   ├── case_3D/ # 3D slab variant
│   │   └── README.md # OpenFOAM case details (meshing/time/BCs)
│   ├── build.sh
│   ├── run_S.sh
│   ├── run_M.sh
│   └── run_L.sh
├── visualization/
│   ├── tsunami.paraview/ # ParaView assets for GeoClaw (optional organization)
│   ├── vortex_shedding.pvsm # ParaView state for OpenFOAM cases
│   └── ... # colormaps, screenshots, state variants
├── results/
│   ├── figures/ # exported frames/plots used in the report
│   └── animations/ # mp4/gif exports (optional)
├── report/
│   ├── main.tex
│   ├── mybibfile.bib
│   └── figures/
└── README.md
```

(Your exact folder names may differ — adjust this tree to match the repo.)

---

## 1) GeoClaw — 1755 Lisbon tsunami (far-field)

### Requirements
- A working GeoClaw installation (Clawpack)
- Python (as required by your GeoClaw workflow)
- ParaView (for Sci-Viz)

### What this module produces
- AMR outputs across refinement levels  
- **Gauges** for mareograms (time series of surface elevation)
- Derived fields for visualization (e.g., flow velocity magnitude)
- Maps such as **Maximum Wave Height (MWH)**

### Typical workflow
1. Configure the rupture/source parameters (Horseshoe Fault) and bathymetry inputs.
2. Run GeoClaw simulation.
3. Convert outputs to ParaView-friendly formats (`.pvd/.vtu` or equivalent).
4. Visualize:
   - AMR grid evolution
   - Free surface elevation
   - Flow velocity magnitude
   - MWH maps

> The report explains the modelling assumptions, AMR usage, gauges, and validation discussion.

---

## 2) OpenFOAM — representative near-field flow analysis

This module packages a ParaView-friendly OpenFOAM vortex-shedding case (2D cylinder in channel) with multiple presets. It is used as a controlled analogue for studying coherent structures and vortex dynamics, and it supports the “Further Discussions” part of the report.

### Requirements
- OpenFOAM installed and sourced  
  Example:
  ```bash
  source /opt/openfoam*/etc/bashrc
  ```
- ParaView for post-processing

### Build
```bash
./simulation/build.sh
```

### Run presets
```bash
./simulation/run_S.sh
./simulation/run_M.sh
./simulation/run_L.sh
```

Outputs are written into:

- `simulation/openfoam/case_S`
- `simulation/openfoam/case_M`
- `simulation/openfoam/case_L`

Each case writes time directories plus a `case.foam` file for ParaView.

Additional cases:

- `simulation/openfoam/case_S_2cyl` (two cylinders)
- `simulation/openfoam/case_3D` (3D slab)

### Configuration
Mesh and time settings live inside each case’s `system/` dictionaries.

See `simulation/openfoam/README.md` for per-case details (meshing, solver settings, BCs).

---

## 3) ParaView — visualization pipeline

### OpenFOAM (vortex shedding)
Open `case.foam` inside a chosen case directory.

Load the ParaView state:

- `visualization/vortex_shedding.pvsm`

Export frames/animations to:

- `results/figures/` and/or `results/animations/`

Typical pipeline elements:

- Vorticity (e.g., ωz)
- Stream tracers
- Velocity glyphs / vectors
- ROI zooms and composite layouts

### GeoClaw (tsunami)
The GeoClaw raw outputs may require conversion before ParaView can load them. The repository includes (or should include) scripts/tools used to convert the ASCII outputs into ParaView-readable formats.

If you have a dedicated converter script, link it here and give a one-command example.

---

## Reproducibility notes
- Far-field tsunami results depend strongly on AMR refinement and bathymetry resolution.
- Near-field OpenFOAM cases are intended as representative numerical experiments, not historical reconstructions.
- For the exact figures and animations used in the paper, see `results/` and `visualization/`.

---

## Citation
If you use this repository in academic work, please cite the associated report:

Diogo Silva, Frederico Afonso, Tomás Pereira.  
*A Sci-Viz Reconstruction and Analysis of the 1755 Lisbon Tsunami.*

---

## Acknowledgements
- GeoClaw / Clawpack community
- ParaView
- OpenFOAM
- GEBCO bathymetry dataset (2023)

---
