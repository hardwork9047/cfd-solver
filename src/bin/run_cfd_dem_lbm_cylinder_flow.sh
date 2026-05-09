#!/usr/bin/env bash
set -euo pipefail

# Fluid-only cylinder flow through the same solver path used by run_lbm_dem.py.
# Extra command-line options are forwarded to the Python module.
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/cfd-solver-matplotlib}"
PYTHONPATH="${PYTHONPATH:-}:src" python3 -m cfd_dem_lbm.cylinder_flow \
  --fluid-method lbm-trt-guo \
  --flow-condition fixed-pressure \
  --result-tag cylinder_flow_shared_solver \
  "$@"
