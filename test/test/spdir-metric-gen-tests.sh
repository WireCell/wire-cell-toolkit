#!/bin/bash

# This is a tad kludgey but we want a main BATS file to run the "spdir-metric"
# test separately for each detector so waf will run each in parallel.  But, we
# do not want to edit each one just to place the detector name where needed.
# So, we let sed do the editing.

tdir="$(dirname "$(realpath "$BASH_SOURCE")")"

for det in pdsp uboone
do
    sed -e s/@DETECTOR@/$det/g < $tdir/spdir-metric.tmpl > $tdir/test-spdir-metric-${det}.bats
done
