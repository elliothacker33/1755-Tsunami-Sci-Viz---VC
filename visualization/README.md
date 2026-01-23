# ParaView pipeline (suggested)

1) Open `case.foam` from `simulation/openfoam/case_S` (or `case_M` / `case_L`).
2) Use `vorticity` for the main color map (diverging, centered at 0).
   - If missing, add a `Calculator` with expression `curl(U)` and name it `vorticity`.
3) Add a `Contour` on `vorticity` with a few levels to highlight vortices.
4) Add `Stream Tracer` using `U` and seed a line just upstream of the cylinder.
5) Save the state as `visualization/vortex_shedding.pvsm`.

Recommended view setup
- Camera: 2D, view along +Z
- Color maps: diverging for `vorticity`
- Use a `Clip` or `Slice` to focus on the wake region

Exports
- Image: `results/vortex.png`
- Animation: `results/vortex.mp4`

Scripted state
```
pvpython visualization/make_state.py simulation/openfoam/case_S/case.foam visualization/vortex_shedding.pvsm results/vortex.png
```
