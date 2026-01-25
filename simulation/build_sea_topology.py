# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Diogo Silva, Frederico Silva, Tomás Pereira

# ============ Build GeoClaw Sea Toplogy (et Baptista)  ================== #
# ==== Authors: Diogo Silva, Frederico Silva, Tomás Pereira ============== #
# ======================================================================== #

from clawpack.geoclaw import dtopotools
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

# =========================== #
#   Make topology of seabed    #
# =========================== #
# The topology defined used 16 seabed faults #
# We use this paper as reference Far field tsunami simulations of the 1755 Lisbon earthquake: Implications for tsunami
# hazard to the U.S. East Coast and the Caribbean
# ============================================== #

# ====================================== #
#  Elsevier plotting configurations      #
# ====================================== #

width_inch = 190 / 25.4
height_inch = width_inch * 0.5

plt.rcParams.update({
    "figure.figsize": (width_inch, height_inch),
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 10,
    "savefig.dpi": 600,
    "figure.autolayout": True
})

# Define faults with multiple subfaults
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

# ======================= #
# Build topology function #
# ======================= #
# It's important to choose an arc-minutes resolution for the grid #

def build_topology(fault_params, fault_name):

    # Fault configuration with multiple subfaults
    fault1 = dtopotools.Fault()

    subfaults_list = []
    segment_length = fault_params['length'] / len(fault_params['subfaults'])

    for i, subfault_data in enumerate(fault_params['subfaults']):
        sub = dtopotools.SubFault()
        sub.strike = fault_params['strike']
        sub.length = segment_length
        sub.width = fault_params['width']
        sub.depth = fault_params['depth']
        sub.slip = subfault_data['slip']
        sub.rake = fault_params['rake']
        sub.dip = fault_params['dip']

        # Calculate position along strike
        strike_rad = np.radians(sub.strike)
        offset = subfault_data['offset']
        cos_lat = np.cos(np.radians(fault_params['latitude']))

        dx = offset * np.sin(strike_rad) / (111000.0 * cos_lat)
        dy = offset * np.cos(strike_rad) / 111000.0

        sub.longitude = fault_params['longitude'] + dx
        sub.latitude = fault_params['latitude'] + dy
        sub.coordinate_specification = 'top center'

        subfaults_list.append(sub)

    fault1.subfaults = subfaults_list
    fault1.rupture_type = "static"

    x = np.linspace(-12, -6, 300)
    y = np.linspace(34, 40, 200)
    time_rupture = [0., 1.]
    fault1.create_dtopography(x, y, time_rupture, verbose=True)

    # Topology file
    dtopo = fault1.dtopo
    topo_name = "topo_lisbon1755_" + fault_name + ".tt3"
    dtopo.write(f"../data/topology/{topo_name}", dtopo_type=3)
    print("Created a file with topologic deformation\n")

    # Plotting
    plotting_topological(fault1, dtopo, x, y, fault_name, fault_params)

# ========================================= #
# All of this pictures are used on the report
# ========================================= #
def plotting_topological(fault, dtopo, x, y, fault_name, fault_params):

    plt.style.use(['science', 'no-latex'])

    # 1. Plot - Deformation Contour lines (2D)
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))

    dz = dtopo.dZ[-1,:,:]

    # Levels
    v_max = np.max(dz)
    v_min = np.min(dz)
    v_limit = max(abs(v_max), abs(v_min))
    v_limit = np.ceil(v_limit) if v_limit > 0 else 1
    levels = np.linspace(-v_limit, v_limit, 21)

    im = ax.contourf(x, y, dz, levels=levels, cmap='RdBu_r', extend='both')

    ax.contour(x, y, dz, levels=np.linspace(0.1, v_limit, 10), colors='black', linestyles='dashed',  linewidths=0.4, alpha=0.8)

    ax.contour(x, y, dz, levels=np.linspace(-v_limit, -0.05, 5), colors='blue', linestyles='dashed', linewidths=0.4, alpha=0.8)

    ax.set_facecolor('#fdfdfd')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Deformation (m)')
    cbar.set_ticks([-v_limit, -v_limit/2, 0, v_limit/2, v_limit])

    ax.set_xlabel(r'Longitude ($^\circ$)')
    ax.set_ylabel(r'Latitude ($^\circ$)')

    ax.grid(True, linestyle=':', alpha=0.4)
    ax.set_aspect('equal')

    plot_name = "fault_deformation_" + fault_name + ".pdf"
    plt.savefig(f"../plots/fault_deformation/{plot_name}", bbox_inches='tight')
    plt.close()

    # 2. Plot the contours on a 3D surface
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(width_inch, height_inch))
    ax = fig.add_subplot(111, projection='3d')

    levels_pos = np.linspace(0.1, v_limit, 10)
    levels_neg = np.linspace(-v_limit, -0.05, 5)

    for level in levels_pos:
        ax.contour(X, Y, dz, levels=[level], colors='red', linewidths=1.5, alpha=0.8)

    for level in levels_neg:
        ax.contour(X, Y, dz, levels=[level], colors='blue', linewidths=1.5, alpha=0.8)

    ax.contour(X, Y, dz, levels=[0], colors='black', linewidths=2.5)
    ax.plot_surface(X, Y, np.zeros_like(dz), alpha=0.1, color='gray')

    ax.set_xlabel('Longitude (°)', fontsize = 9)
    ax.set_ylabel('Latitude (°)', fontsize = 9)
    ax.set_zlabel('Deformation (m)', fontsize =9)
    ax.zaxis.label.set_rotation(90)
    ax.set_box_aspect(None, zoom=0.7)
    ax.grid(True, linestyle=':', alpha=0.4)

    import matplotlib.patches as mpatches

    red_patch = mpatches.Patch(color='red', alpha=0.6,   label='Uplift')
    blue_patch = mpatches.Patch(color='blue', alpha=0.6, label='Subsidence')
    black_patch = mpatches.Patch(color='black', alpha=0.6, label='Zero deformation')

    ax.legend(handles=[red_patch, blue_patch, black_patch], loc='upper right', fontsize = 8)
    ax.view_init(elev=25, azim=135)

    plot_name = "fault_deformation_3D_" + fault_name + ".pdf"
    plt.savefig(f"../plots/fault_deformation/{plot_name}", bbox_inches='tight')
    plt.close()

    # 3. Validation statistics
    validation_stats(dtopo, dz, x, y, fault_params, fault_name)


# ========================================= #
# Validation statistics                     #
# ========================================= #
def validation_stats(dtopo, dz, x, y, fault_params, fault_name):

    # 1. Histogram
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))

    n, bins, patches_hist = ax.hist(dz.flatten(), bins=50,
                                      color='blue',
                                      edgecolor='black',
                                      linewidth=0.8,
                                      alpha=0.2)


    ax.set_xlabel('Deformation (m)')
    ax.set_ylabel('Frequency')
    ax.set_title('Deformation Distribution')
    ax.grid(True, linestyle=':', alpha=0.3, axis='y')

    plt.savefig(f"../plots/validation/histogram_{fault_name}.pdf", bbox_inches='tight')
    plt.close()

    # 2. E-W Profile
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))

    mid_lat = dz.shape[0] // 2
    ax.plot(x, dz[mid_lat, :], 'b-', linewidth=1.5, label='E-W Profile')
    ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.fill_between(x, 0, dz[mid_lat, :], where=(dz[mid_lat, :] > 0),
                     color='red', alpha=0.2, label='Uplift')
    ax.fill_between(x, 0, dz[mid_lat, :], where=(dz[mid_lat, :] < 0),
                     color='blue', alpha=0.2, label='Subsidence')

    ax.set_xlabel(r'Longitude ($^\circ$)')
    ax.set_ylabel('Deformation (m)')
    ax.set_title('East-West Profile')
    ax.grid(True, linestyle=':', alpha=0.3)
    ax.legend(loc='best', framealpha=0.9)

    plt.savefig(f"../plots/validation/profile_ew_{fault_name}.pdf", bbox_inches='tight')
    plt.close()

    # 3. N-S Profile
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))

    mid_lon = dz.shape[1] // 2
    ax.plot(y, dz[:, mid_lon], 'r-', linewidth=1.5, label='N-S Profile')
    ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.fill_between(y, 0, dz[:, mid_lon], where=(dz[:, mid_lon] > 0),
                     color='red', alpha=0.2, label='Uplift')
    ax.fill_between(y, 0, dz[:, mid_lon], where=(dz[:, mid_lon] < 0),
                     color='blue', alpha=0.2, label='Subsidence')

    ax.set_xlabel(r'Latitude ($^\circ$)')
    ax.set_ylabel('Deformation (m)')
    ax.set_title('North-South Profile')
    ax.grid(True, linestyle=':', alpha=0.3)
    ax.legend(loc='best', framealpha=0.9)

    plt.savefig(f"../plots/validation/profile_ns_{fault_name}.pdf", bbox_inches='tight')
    plt.close()



# ============================= #
# ============================= #
# Menu - Choose Fault Topology  #
# ============================= #
# 1. Horshe-Shoe Fault - Normally goes more in accord with the
# reports of magnitude          #
# 2. MPF - Marques de Pombal    #
# ============================= #

while True:

    # Menu Header
    print("============================================\n")
    print("                                            \n")
    print("               Choose Topology              \n")
    print("                                            \n")
    print("             1. Horse-Shoe Fault            \n")
    print("             2. Marques de Pombal Fault     \n")
    print("             3. Exit Menu                   \n")
    print("                                            \n")
    print("============================================\n")

    # Menu question
    try:
        user_input = input("Insert topology you want: ")
        topo = int(user_input)
    except ValueError:
        print("Insert a valid number, considering the options available\n")
        continue

    if topo == 3:
        break

    if 1 <= topo <= len(faults):
        build_topology(*faults[topo - 1])
