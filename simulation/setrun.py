# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Diogo Silva, Frederico Silva, Tomás Pereira

# ============ Build GeoClaw Sea Simulation ============================== #
# ==== Authors: Diogo Silva, Frederico Silva, Tomás Pereira ============== #
# ======================================================================== #

# === This file is adapted from the original geoclaw 2010 Chile Tsunami == #

from __future__ import absolute_import
from __future__ import print_function

import sys
from clawpack.geoclaw import dtopotools
import numpy as np
import matplotlib.pyplot as plt
import scienceplots


hsf_fault = {
     "longitude": -9.91,
     "latitude": 35.74,
     "rake": 90,
     "length": 165000.0,
     "width": 70000.0,
     "dip": 35,
     "slip": 10.7,
     "strike": 42.1,
     "depth": 4000.0,
     "subfaults": [
         {"slip": 8.5, "offset": -55000.0},
         {"slip": 10.7, "offset": 0.0},
         {"slip": 9.2, "offset": 55000.0}
     ]
 }


mpf_fault = {
     "longitude": -9.89,
     "latitude": 36.57,
     "rake": 90,
     "length": 129000.0,
     "width": 70000.0,
     "dip": 35,
     "slip": 8.0,
     "strike": 20.0,
     "depth": 4000.0,
     "subfaults": [
         {"slip": 6.5, "offset": -43000.0},
         {"slip": 8.0, "offset": 0.0},
         {"slip": 7.0, "offset": 43000.0}
     ]
 }

# List of default faults that can be selected
faults = [(hsf_fault, "hsf"), (mpf_fault, "mpf")]



# =========================== #
#       Set run simulator     #
# =========================== #
def setrun(claw_pkg="geoclaw", fault_params=None, fault_name=None):

    from clawpack.clawutil import data
    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    # NUmber of dimnesions (2D + 1 Dimension(Bathymetry))

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters
    #------------------------------------------------------------------
    probdata = rundata.new_UserData(name='probdata', fname='setprob.data')

    #------------------------------------------------------------------
    # Standard Clawpack parameters
    #------------------------------------------------------------------
    clawdata = rundata.clawdata

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -12.0   # west longitude
    clawdata.upper[0] = -6.0    # east longitude

    clawdata.lower[1] = 36.0    # south latitude
    clawdata.upper[1] = 40.0    # north latitude

    # AMR 1 Level resolution
    clawdata.num_cells[0] = 45
    clawdata.num_cells[1] = 30

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Number of waves
    clawdata.num_waves = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # ----------------------
    # Simulation time:
    # ----------------------
    clawdata.num_output_times = 24
    clawdata.t0 = 0.0
    clawdata.tfinal = 7200.0
    clawdata.output_style = 1
    clawdata.output_t0 = True

    clawdata.output_format = 'ascii'
    clawdata.output_q_components = 'all'
    clawdata.output_aux_components = 'none'
    clawdata.output_aux_onlyonce = True

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Numerical method
    clawdata.order = 2
    clawdata.dimensional_split = 'unsplit'
    clawdata.transverse_waves = 2
    clawdata.limiter = ['mc', 'mc', 'mc']
    clawdata.use_fwaves = True

    clawdata.source_split = 'godunov'
    clawdata.cfl_desired = 0.8
    clawdata.cfl_max = 1.0

    clawdata.steps_max = 20000

    clawdata.dt_variable = True
    clawdata.dt_initial = 2.0
    clawdata.dt_max = 1.e99

    rundata = setgeo(rundata, fault_params, fault_name)

    return rundata

# ------------------ #
# Amr and general -- #
# ------------------ #
def setgeo(rundata, fault_params, fault_name):
    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # ------------------------ #
    # General Earth Settings - #
    # ------------------------ #

    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1000.0

    rundata.amrdata.amr_levels_max = 4

    # ============================================================================ #
    #                     AMR REFINEMENT ZONES - LISBON 1755 TSUNAMI              #
    # ============================================================================ #

    rundata.amrdata.refinement_ratios_x = [4, 5, 4]
    rundata.amrdata.refinement_ratios_y = [4, 5, 4]
    rundata.amrdata.refinement_ratios_t = [4, 5, 4]

    rundata.amrdata.flag_richardson = False
    rundata.amrdata.flag2refine = True
    rundata.amrdata.regrid_interval = 3
    rundata.amrdata.regrid_buffer_width = 2
    rundata.amrdata.clustering_cutoff = 0.7
    rundata.amrdata.verbosity_regrid = 0

    rundata.amrdata.aux_type = ['center', 'capacity', 'center']

    rundata.regiondata.regions = []

    # ============================================================================ #
    # LEVEL 4 - CRITICAL COASTAL ZONES
    # ============================================================================ #

    # Region 4.1: Lisbon Metropolitan Area
    # Includes Baixa, Alfama, Belém
    rundata.regiondata.regions.append([
        4, 4,
        1200.0, 1e9,
        -9.25, -9.05,
        38.68, 38.76
    ])

    # Region 4.2: Cascais Bay
    rundata.regiondata.regions.append([
        4, 4,
        1200.0, 1e9,
        -9.45, -9.38,
        38.68, 38.73
    ])

    # Region 4.3: Lagos (Algarve)
    rundata.regiondata.regions.append([
        4, 4,
        1800.0, 1e9,
        -8.70, -8.65,
        37.08, 37.12
     ])

    # ============================================================================ #
    # LEVEL 3 - NEAR-SHORE ZONES
    # ============================================================================ #

    # Region 3.1: Portuguese Coast - Lisbon District
    rundata.regiondata.regions.append([
        3, 3,
        900.0, 1e9,
        -9.50, -8.90,
        38.50, 39.00
    ])

    # Region 3.2: Setúbal Peninsula
    rundata.regiondata.regions.append([
        3, 3,
        1200.0, 1e9,
        -9.20, -8.80,
        38.40, 38.60
    ])

    # Region 3.3: Algarve Coast
    rundata.regiondata.regions.append([
        3, 3,
        1500.0, 1e9,
        -9.00, -7.50,
        36.95, 37.20
    ])

    # Region 3.4: Cadiz Bay (Spain)
    rundata.regiondata.regions.append([
        3, 3,
        1800.0, 1e9,
        -6.50, -6.10,
        36.40, 36.60
    ])

    # Region 3.5: Sines Harbor
    rundata.regiondata.regions.append([
        3, 3,
        1500.0, 1e9,
        -8.90, -8.85,
        37.93, 37.97
    ])

    # Region 3.6: Peniche
    rundata.regiondata.regions.append([
        3, 3,
        1200.0, 1e9,
        -9.40, -9.35,
        39.33, 39.37
    ])

    # ============================================================================ #
    # LEVEL 2 - SOURCE AND PROPAGATION ZONES
    # ============================================================================ #

    # Region 2.1: Horseshoe Fault Zone
    rundata.regiondata.regions.append([
        2, 2,
        0.0, 450.0,
        -11.50, -8.50,
        34.50, 37.50
    ])

    # Region 2.2: Continental Shelf (Portugal)
    rundata.regiondata.regions.append([
        2, 2,
        300.0, 1e9,
        -10.50, -8.50,
        36.50, 39.50
    ])

    # Region 2.3: Gulf of Cadiz
    rundata.regiondata.regions.append([
        2, 2,
        900.0, 1e9,
        -9.50, -6.00,
        36.00, 37.50
    ])

    # ============================================================================ #
    # GAUGE LOCATIONS (Validation Points)
    # ============================================================================ #

    rundata.gaugedata.gauges = []

    # Gauge 1: Lisbon - Terreiro do Paço
    rundata.gaugedata.gauges.append([1, -9.133, 38.708, 0.0, 1e9])

    # Gauge 2: Cascais Harbor
    rundata.gaugedata.gauges.append([2, -9.420, 38.697, 0.0, 1e9])

    # Gauge 3: Tagus River - Cacilhas
    rundata.gaugedata.gauges.append([3, -9.155, 38.685, 0.0, 1e9])

    # Gauge 4: Setúbal
    rundata.gaugedata.gauges.append([4, -8.892, 38.524, 0.0, 1e9])

    # Gauge 5: Sines
    rundata.gaugedata.gauges.append([5, -8.868, 37.950, 0.0, 1e9])

    # Gauge 6: Lagos (Algarve)
    rundata.gaugedata.gauges.append([6, -8.674, 37.101, 0.0, 1e9])

    # Gauge 7: Faro (Algarve)
    rundata.gaugedata.gauges.append([7, -7.935, 37.015, 0.0, 1e9])

    # Gauge 8: Cadiz (Spain)
    rundata.gaugedata.gauges.append([8, -6.288, 36.533, 0.0, 1e9])

    # Gauge 9: Offshore
    rundata.gaugedata.gauges.append([9, -10.000, 37.500, 0.0, 1e9])

    # Gauge 10: Source Region (Near fault)
    rundata.gaugedata.gauges.append([10, -10.500, 36.000, 0.0, 1e9])

    # Gauge 11: Forte do Bugio
    rundata.gaugedata.gauges.append([11, -9.31, 38.65, 0.0, 1e9])

    topo_path = '../data/'

    rundata.topo_data.topofiles = []
    rundata.topo_data.topofiles.append([3, 1, 1, 0.0, 1e9, topo_path + 'GEBCO_data/gebco_coarse_data.asc'])
    rundata.topo_data.topofiles.append([3, 1, 2, 0.0, 1e9, topo_path + 'GEBCO_data/gebco_medium_data.asc'])
    rundata.topo_data.topofiles.append([3, 1, 3, 0.0, 1e9, topo_path + 'GEBCO_data/gebco_fine_data.asc'])

    rundata.dtopo_data.dtopofiles = []
    rundata.dtopo_data.dtopofiles.append(
        [3, topo_path + f'topology/topo_lisbon1755_{fault_name}.tt3']
    )

    rundata.dtopo_data.dt_max_dtopo = 1.0

    return rundata


# ============================= #
# Menu - Run simulation         #
# ============================= #
# 1. Horshe-Shoe Fault          #
# reports of magnitude          #
# 2. MPF - Marques de Pombal    #
# ============================= #

if __name__ == '__main__':

    # Menu Header
    print("============================================\n")
    print("                                            \n")
    print("               Run Simulation               \n")
    print("                                            \n")
    print("             1. Horse-Shoe Fault            \n")
    print("             2. Marques de Pombal Fault     \n")
    print("                                            \n")
    print("============================================\n")

    # Menu question
    try:
        user_input = input("Insert topology you want: ")
        topo = int(user_input)
    except ValueError:
        print("Insert a valid number, considering the options available\n")

    if 1 <= topo <= len(faults):
        fault_params, fault_name = faults[topo - 1]
        rundata = setrun(fault_params=fault_params, fault_name=fault_name)
        rundata.write()
