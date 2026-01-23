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
