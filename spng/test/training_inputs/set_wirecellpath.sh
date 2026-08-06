#!/bin/bash
# Automatically determine paths based on script location

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Navigate up to the toolkit root (3 levels up from spng/test/training_inputs)
TOOLKIT_ROOT="$( cd "${SCRIPT_DIR}/../../.." && pwd )"

# Build WIRECELL_PATH using relative paths
export WIRECELL_PATH="${SCRIPT_DIR}:${TOOLKIT_ROOT}/wire-cell-data:${TOOLKIT_ROOT}/spng/cfg:${TOOLKIT_ROOT}/spng/test:${TOOLKIT_ROOT}/build/root:$WIRECELL_PATH"