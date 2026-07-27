#!/bin/bash
export SHELL="/bin/bash"
set +euo pipefail  
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw v10_14_00d00 -q e26:prof
set -euo pipefail
lar "$@"
