# OpenFOAM cases (cylinder vortex shedding)

This folder contains OpenFOAM cases:
- `case_S` (fast preview)
- `case_M` (report-quality)
- `case_L` (stress test)
- `case_S_2cyl` (two cylinders)
- `case_3D` (3D slab with front/back walls)

Each case includes a 2D channel with a cylinder (snappyHexMesh) and a
transient laminar run using `pimpleFoam`.

## Requirements
- OpenFOAM installed and sourced (e.g. `source /opt/openfoam*/etc/bashrc`)
- ParaView (or `paraFoam`) for visualization

## Run
From repo root:
```
./simulation/run_S.sh
./simulation/run_M.sh
./simulation/run_L.sh
```

Outputs are written inside each case directory. Open the case in ParaView
by loading the `case.foam` file (created after running).

Extra cases
```
./simulation/openfoam/case_S_2cyl/Allrun
./simulation/openfoam/case_3D/Allrun
```

## Clean
Inside a case folder:
```
./Allclean
```
