#!/usr/bin/env pvpython
"""
ParaView pipeline for GeoClaw VTK outputs.

Usage:
  pvpython chile2010/paraview_pipeline.py
"""

from __future__ import print_function

import glob
import os
import sys

from paraview.simple import (  # type: ignore
    CellDatatoPointData,
    Calculator,
    ColorBy,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    GetTimeKeeper,
    OpenDataFile,
    Render,
    SaveScreenshot,
    SaveState,
    Show,
    Threshold,
    WarpByScalar,
)


# -----------------------------
# Config
# -----------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
VTK_DIR = os.environ.get("GEOCLAW_VTK_DIR") or os.environ.get(
    "GEoCLAW_VTK_DIR", os.path.join(SCRIPT_DIR, "vtk")
)
VTHB_PATTERN = os.path.join(VTK_DIR, "geoclaw*.vthb")

COLOR_BY = "eta"  # "eta" or "speed"
WARP_SCALE = 50.0  # exaggeration for eta (meters vs degrees)
WATER_DEPTH_MIN = 1.0  # meters, filters out dry cells

ETA_RANGE = (-2.0, 2.0)  # meters
SPEED_RANGE = (0.0, 5.0)  # m/s

SAVE_STATE = True
SAVE_SCREENSHOT = True
STATE_PATH = os.path.join(SCRIPT_DIR, "paraview_pipeline.pvsm")
SCREENSHOT_PATH = os.path.join(SCRIPT_DIR, "paraview_snapshot.png")
SCREENSHOT_SIZE = [1400, 900]


def main() -> int:
    files = sorted(glob.glob(VTHB_PATTERN))
    if not files:
        print("No .vthb files found in:", VTK_DIR)
        print("Expected pattern:", VTHB_PATTERN)
        return 1

    # Load time series
    reader = OpenDataFile(files)

    # Keep only wet cells (q_0 = depth)
    threshold = Threshold(Input=reader)
    threshold.Scalars = ["CELLS", "q_0"]
    try:
        threshold.ThresholdRange = [WATER_DEPTH_MIN, 1.0e9]
    except Exception:
        threshold.LowerThreshold = WATER_DEPTH_MIN
        threshold.UpperThreshold = 1.0e9
        threshold.ThresholdMethod = "Between"

    # Convert to point data for smooth coloring/warp
    cell_to_point = CellDatatoPointData(Input=threshold)

    # Compute speed = sqrt(hu^2 + hv^2) / h
    calc = Calculator(Input=cell_to_point)
    calc.ResultArrayName = "speed"
    calc.Function = "sqrt(q_1*q_1 + q_2*q_2) / q_0"

    # Warp by eta
    warp = WarpByScalar(Input=calc)
    warp.Scalars = ["POINTS", "q_3"]
    warp.ScaleFactor = WARP_SCALE

    # Render setup
    render_view = GetActiveViewOrCreate("RenderView")
    display = Show(warp, render_view)
    display.Representation = "Surface"

    if COLOR_BY == "speed":
        ColorBy(display, ("POINTS", "speed"))
        lut = GetColorTransferFunction("speed")
        lut.RescaleTransferFunction(SPEED_RANGE[0], SPEED_RANGE[1])
    else:
        ColorBy(display, ("POINTS", "q_3"))
        lut = GetColorTransferFunction("q_3")
        lut.RescaleTransferFunction(ETA_RANGE[0], ETA_RANGE[1])

    # Keep opacity sane
    otf = GetOpacityTransferFunction("q_3")
    otf.RescaleTransferFunction(ETA_RANGE[0], ETA_RANGE[1])

    render_view.ResetCamera()
    Render()

    # Save state and a snapshot of the last timestep
    if SAVE_STATE:
        SaveState(STATE_PATH)

    if SAVE_SCREENSHOT:
        time_keeper = GetTimeKeeper()
        if time_keeper.TimestepValues:
            render_view.ViewTime = time_keeper.TimestepValues[-1]
        Render()
        SaveScreenshot(SCREENSHOT_PATH, render_view, ImageResolution=SCREENSHOT_SIZE)

    print("Done.")
    print("State:", STATE_PATH)
    print("Screenshot:", SCREENSHOT_PATH)
    return 0


if __name__ == "__main__":
    sys.exit(main())
