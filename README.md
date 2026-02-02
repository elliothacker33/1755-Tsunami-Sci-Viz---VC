MIT License

Copyright (c) 2026 Diogo Silva, Frederico Silva, Tomás Pereira

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# 1755-Tsunami-Sci-Viz---VC
This repository contains all the files and plots used to simulate the 1755 Lisbon Tsunami.

# vc-sciviz: von Karman vortex street (OpenFOAM)

This repo packages an OpenFOAM vortex-shedding case (2D cylinder in a channel)
with S/M/L presets and ParaView-friendly outputs.

## Requirements
- OpenFOAM installed and sourced (e.g. `source /opt/openfoam*/etc/bashrc`)
- ParaView for visualization

## Build
```
./simulation/build.sh
```

## Run presets
```
./simulation/run_S.sh
./simulation/run_M.sh
./simulation/run_L.sh
```

Outputs are written inside:
- `simulation/openfoam/case_S`
- `simulation/openfoam/case_M`
- `simulation/openfoam/case_L`

Each case writes time directories plus a `case.foam` file for ParaView.

Additional cases live in:
- `simulation/openfoam/case_S_2cyl` (two cylinders)
- `simulation/openfoam/case_3D` (3D slab)

## Configuration
- Mesh and time settings are defined in each case's `system/` dictionaries.
- See `simulation/openfoam/README.md` for details.

## ParaView
- Open `case.foam` inside the chosen case directory.
- Build your pipeline (vorticity, streamlines) and save it as
  `visualization/vortex_shedding.pvsm`.
- Export key frames or animations to `results/`.

## Notes
- Legacy Basilisk files remain in the repo for reference but are not used.
