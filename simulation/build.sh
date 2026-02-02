#!/usr/bin/env bash
set -euo pipefail

if ! command -v blockMesh >/dev/null 2>&1; then
  echo "OpenFOAM not in PATH. Source your OpenFOAM bashrc first." >&2
  exit 1
fi

echo "OpenFOAM environment detected."
