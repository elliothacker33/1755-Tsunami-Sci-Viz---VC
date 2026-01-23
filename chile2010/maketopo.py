# make_dtopo_1755.py
import os
import numpy as np
from clawpack.geoclaw import dtopotools

# -----------------------------
# 1) Output file
# -----------------------------
outdir = os.path.join(os.path.dirname(__file__), "dtopo")
os.makedirs(outdir, exist_ok=True)
dtopo_fname = os.path.join(outdir, "dtopo_1755.tt3")

# -----------------------------
# 2) Define the fault (Okada)
# -----------------------------
# NOTE: These values are a *starting point* (you MUST adjust to your chosen 1755 scenario).
# Units:
#  - lon/lat in degrees
#  - depth in meters (positive down)
#  - length/width in meters
#  - slip in meters
# Angles in degrees: strike, dip, rake

fault = dtopotools.Fault()

subfault = dtopotools.SubFault()
subfault.longitude = -10.0      # center lon (example)
subfault.latitude  =  36.5      # center lat (example)
subfault.depth     =  15000.0   # top depth (m), example
subfault.strike    =  60.0      # example
subfault.dip       =  20.0      # example
subfault.rake      =  90.0      # thrust
subfault.length    =  180e3     # 180 km
subfault.width     =   80e3     # 80 km
subfault.slip      =   10.0     # 10 m
subfault.mu        =  3.0e10    # shear modulus (Pa), ok default
subfault.coordinate_specification = "top center"

fault.subfaults = [subfault]

# -----------------------------
# 3) Define the dtopo grid + times
# -----------------------------
dtopo = dtopotools.DTopography()

# Your topo domain is [-16,-6] x [34,42]. Use the same (or a slightly larger) box.
x1, x2 = -16.0, -6.0
y1, y2 =  34.0, 42.0

# Resolution for deformation grid (coarser than topo is fine)
# ~0.01 deg ~ 1.1 km in latitude
dx = dy = 0.01

mx = int(round((x2 - x1) / dx)) + 1
my = int(round((y2 - y1) / dy)) + 1

dtopo.x = np.linspace(x1, x2, mx)
dtopo.y = np.linspace(y1, y2, my)

# Times: for an earthquake source you usually want *instantaneous* seafloor deformation.
# A common trick: provide t=[0, 1] with same deformation (a step).
dtopo.times = np.array([0.0, 1.0])

# -----------------------------
# 4) Compute deformation from fault and write TT3
# -----------------------------
fault.create_dtopography(dtopo.x, dtopo.y, dtopo.times, verbose=True)
dtopo = fault.dtopo
dtopo.write(dtopo_fname, dtopo_type=3)

print("Wrote:", dtopo_fname)
print("Grid:", mx, "x", my, "times:", dtopo.times)
