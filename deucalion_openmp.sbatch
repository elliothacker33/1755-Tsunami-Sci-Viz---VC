#!/bin/bash
#SBATCH -A f202500010hpcvlabuminhoa
#SBATCH -p normal-arm
#SBATCH -t 00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=tests/slurm_logs/test_out.o%j
#SBATCH --error=tests/slurm_logs/test_err.e%j
#SBATCH --exclusive
#SBATCH --acctg-freq=energy=10

set -euo pipefail

REPO_DIR="$HOME/1755-Tsunami-Sci-Viz---VC/chile2010"
CLAW_DIR="$HOME/clawpack"
VENV_DIR="$REPO_DIR/.venv"

module purge
# module load gcc python  # adjust for Deucalion modules

if [ -d "$VENV_DIR" ]; then
  source "$VENV_DIR/bin/activate"
fi

export PYTHONPATH="$CLAW_DIR:$PYTHONPATH"
export CLAW="$CLAW_DIR"
export CLAW_PYTHON="$(which python)"
export FC=gfortran

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export MPLBACKEND=Agg

export FFLAGS="-O2 -fopenmp"
export LFLAGS="-fopenmp"

cd "$REPO_DIR"

python maketopo.py
make clobber
make .output
