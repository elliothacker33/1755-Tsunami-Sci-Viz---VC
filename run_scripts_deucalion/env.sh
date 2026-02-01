#!/bin/bash

# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Diogo Silva, Frederico Afonso, Tomás Pereira
#
# ======================================================================================================== #
# ============ This file is configured to setup a python virtual environment to run all tests ============ #
# ======================================================================================================== #

#SBATCH --job-name=setup_env
#SBATCH --partition=dev-arm
#SBATCH --account=f202500010hpcvlabuminhoa
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:02:00
#SBATCH --output=slurm_bin/job_%j.out
#SBATCH --error=slurm_bin/job_%j.err

# ===================================== #
#            Load modules               #
# ===================================== #

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_PATH="$PROJECT_ROOT/venv_deucalion_os"

 if [ ! -d "$VENV_PATH" ]; then
     echo "Error: Python virtual environment was not found."
     echo "Creating a virtual environment in: $VENV_PATH"
     python3 -m venv "$VENV_PATH"
     source "$VENV_PATH/bin/activate"
     pip install --upgrade pip
     pip install "numpy<2"
     pip install matplotlib scienceplots
     pip install setuptools wheel meson-python
     pip install --no-build-isolation clawpack -v
     pip install netcdf4 cmocean cartopy xarray
 else
     echo "Virtual environment already exists."
     source "$VENV_PATH/bin/activate"
fi


