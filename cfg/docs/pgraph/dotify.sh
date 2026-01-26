#/bin/bash
mydir="$(realpath "$(dirname "$BASH_SOURCE")")"
export WIRECELL_PATH="$(realpath "$mydir/../..")"

out=$1
base=${out%.*}
wcpy pgraph dotify --no-params --no-services --graph-options rankdir=TB  --graph-options bgcolor=transparent $base.jsonnet $out

