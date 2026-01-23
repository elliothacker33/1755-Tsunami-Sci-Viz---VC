# von Karman vortex street (OpenFOAM)

## Goals
- Simulate 2D flow around a cylinder and observe vortex shedding.
- Export data for ParaView and build a visualization pipeline.
- Assess visualization insight and performance as data volume grows.

## Simulation setup
- Solver: OpenFOAM `pimpleFoam` (laminar, transient incompressible).
- Mesh: `blockMesh` base + `snappyHexMesh` around a cylinder (D = 0.125).
- Domain: 8 x 1 channel, cylinder centered at the inlet (x = 0, y = 0).
- Inlet: uniform velocity U = 1; outlet: fixed pressure.
- Reynolds number: 160 (nu = 7.8125e-4).

## Instrumentation
- Outputs: `p`, `U`, and `vorticity` (function object in `controlDict`).
- Format: native OpenFOAM time directories (`case.foam` for ParaView).
- Output cadence: `writeInterval` in each case's `controlDict`.

## Visualization pipeline
- `vorticity` colormap (diverging) to highlight vortices.
- Contours of `vorticity` to track shedding.
- Streamlines built from `U`.

## Insights
- Onset of shedding (transient to periodic regime).
- Wake structure and alternating vortices.
- Dependence of vortex strength and wavelength on mesh/time resolution.

## Performance assessment
| Case | blockMesh (nx x ny) | snappy level | deltaT | endTime | Output size | Wall time |
|------|----------------------|-------------|--------|---------|-------------|-----------|
| S    | 160 x 40             | 2           | 0.01   | 6       |             |           |
| M    | 320 x 80             | 3           | 0.005  | 15      |             |           |
| L    | 640 x 160            | 4           | 0.0025 | 20      |             |           |

Notes:
- Track total file count and disk usage in each case directory.
- Comment on compute vs I/O costs as outputs increase.

## Results
- Image: `results/vortex.png`
- Animation: `results/vortex.mp4`

## Conclusions
- Summarize the main insights and how visualization supports them.
- Discuss scalability and practical output settings.
