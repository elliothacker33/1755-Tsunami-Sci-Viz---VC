#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CASE_DIR="$ROOT_DIR/simulation/openfoam/case_L"

if [[ ! -x "$CASE_DIR/Allrun" ]]; then
  echo "OpenFOAM case not found: $CASE_DIR" >&2
  exit 1
fi

echo "Running OpenFOAM L case. Logs will stream here and be saved in $CASE_DIR/log.*"

"$CASE_DIR/Allrun"
