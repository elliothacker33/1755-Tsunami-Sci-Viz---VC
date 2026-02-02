#!/usr/bin/env python3
import os
import sys

try:
    from paraview.simple import (
        _DisableFirstRenderCameraReset,
        Calculator,
        ColorBy,
        Contour,
        GetActiveViewOrCreate,
        GetColorTransferFunction,
        Hide,
        OpenFOAMReader,
        Render,
        ResetCamera,
        SaveScreenshot,
        SaveState,
        Show,
        CellDatatoPointData,
        Slice,
        StreamTracer,
        Tube,
    )
    try:
        from paraview.simple import AnnotateTimeFilter
    except Exception:
        AnnotateTimeFilter = None
except Exception as exc:
    raise SystemExit(f"ParaView Python not available: {exc}")


def main():
    if len(sys.argv) < 2:
        print("Usage: pvpython visualization/make_state.py /path/to/case.foam [output.pvsm] [screenshot.png]")
        return 1

    case_path = os.path.abspath(sys.argv[1])
    if not os.path.exists(case_path):
        raise SystemExit(f"Case file not found: {case_path}")

    default_state = os.path.join(os.path.dirname(__file__), "vortex_shedding.pvsm")
    state_path = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else default_state
    screenshot_path = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else None

    _DisableFirstRenderCameraReset()

    reader = OpenFOAMReader(FileName=case_path)
    reader.MeshRegions = ["internalMesh"]
    reader.CellArrays = ["U", "p", "vorticity"]
    reader.UpdatePipeline()

    slice_plane = Slice(Input=reader)
    slice_plane.SliceType = "Plane"
    slice_plane.SliceType.Origin = [3.5, 0.0, 0.0]
    slice_plane.SliceType.Normal = [0.0, 0.0, 1.0]

    cell_to_point = CellDatatoPointData(Input=slice_plane)
    cell_to_point.PassCellData = 1

    vort_calc = Calculator(Input=cell_to_point)
    vort_calc.ResultArrayName = "vortZ"

    def _has_array(src, name):
        info = src.GetPointDataInformation()
        if info and info.HasArray(name):
            return True
        info = src.GetCellDataInformation()
        return bool(info and info.HasArray(name))

    def _set_calc(expr):
        vort_calc.Function = expr
        try:
            vort_calc.UpdatePipeline()
        except Exception:
            return False
        return _has_array(vort_calc, "vortZ")

    if not _set_calc("vorticity_Z"):
        if not _set_calc("vorticity[2]"):
            _set_calc("mag(vorticity)")

    contour = Contour(Input=vort_calc)
    contour.ContourBy = ["POINTS", "vortZ"]
    contour.Isosurfaces = [-6, -4, -2, 2, 4, 6]

    stream = StreamTracer(Input=cell_to_point, SeedType="Line")
    stream.Vectors = ["POINTS", "U"]
    stream.SeedType.Point1 = [-0.4, -0.2, 0.0]
    stream.SeedType.Point2 = [-0.4, 0.2, 0.0]
    stream.SeedType.Resolution = 30
    stream.MaximumStreamlineLength = 10.0
    stream.IntegrationDirection = "FORWARD"

    tube = Tube(Input=stream)
    tube.Radius = 0.003
    tube.NumberofSides = 10

    time_annot = AnnotateTimeFilter(Input=vort_calc) if AnnotateTimeFilter else None

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [1600, 900]
    view.Background = [0.06, 0.07, 0.09]
    if hasattr(view, "BackgroundColorMode"):
        view.BackgroundColorMode = "Gradient"
    elif hasattr(view, "UseGradientBackground"):
        view.UseGradientBackground = 1
    if hasattr(view, "Background2"):
        view.Background2 = [0.015, 0.02, 0.03]

    vort_display = Show(vort_calc, view)
    vort_display.Representation = "Surface"
    ColorBy(vort_display, ("POINTS", "vortZ"))
    vort_lut = GetColorTransferFunction("vortZ")
    vort_lut.ApplyPreset("Cool to Warm", True)
    vort_lut.RescaleTransferFunction(-8, 8)
    vort_display.SetScalarBarVisibility(view, True)

    contour_display = Show(contour, view)
    contour_display.Representation = "Surface"
    contour_display.DiffuseColor = [0.95, 0.95, 0.95]
    contour_display.Opacity = 0.6

    tube_display = Show(tube, view)
    tube_display.Representation = "Surface"
    tube_display.DiffuseColor = [1.0, 0.9, 0.75]
    tube_display.Specular = 0.4
    tube_display.SpecularPower = 30.0

    if time_annot:
        time_display = Show(time_annot, view)
        time_display.FontSize = 16
        time_display.Color = [0.9, 0.9, 0.9]
        time_display.WindowLocation = "Upper Left Corner"

    Hide(reader, view)
    Hide(slice_plane, view)
    Hide(cell_to_point, view)

    ResetCamera(view)
    view.CameraParallelProjection = 1
    view.CameraPosition = [3.5, 0.0, 4.0]
    view.CameraFocalPoint = [3.5, 0.0, 0.0]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelScale = 0.65

    Render()
    SaveState(state_path)

    if screenshot_path:
        SaveScreenshot(screenshot_path, view, ImageResolution=view.ViewSize)

    print(f"Saved state: {state_path}")
    if screenshot_path:
        print(f"Saved screenshot: {screenshot_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
