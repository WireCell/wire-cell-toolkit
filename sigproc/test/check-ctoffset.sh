#!/bin/bash
# Drift-distance closure check of the OmnibusSigProc collection time offset
# (ctoffset) for FD-VD.  See check-ctoffset.jsonnet for the method.
#
# usage: check-ctoffset.sh [OUTDIR]     (default: ./check-ctoffset-out)
#
# Needs wire-cell in PATH and WIRECELL_PATH holding the WCT cfg/ and
# wire-cell-data (wires, field and electronics responses, noise spectra).
#
# Runs both drift volumes x {ctoffset from the formula, the PD-VD data value
# 4 us} x {no noise, noise} in parallel (~2 min), then checks:
#   - formula: |bias| <= 0.2 cm and spread over distance <= 0.1 cm
#   - 4 us:    bias of about +1.4 cm (shows the check detects a wrong offset)
# Exits non-zero if any check fails.  Output: OUTDIR/{*.tar.bz2,*.log,results.txt}

set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outdir="${1:-check-ctoffset-out}"
mkdir -p "$outdir"

pids=()
for drift in bottom top; do
    for ct in formula 4; do
        for noise in false true; do
            name="$drift-ct$ct-noise$noise"
            wire-cell -l "$outdir/$name.log" -L info \
                      -A drift=$drift -A ctoffset=$ct -A noise=$noise \
                      -A outname="$outdir/$name.tar.bz2" \
                      -c "$here/check-ctoffset.jsonnet" > "$outdir/$name.out" 2>&1 &
            pids+=($!)
        done
    done
done
status=0
for pid in "${pids[@]}"; do
    wait $pid || status=1
done
if [ $status -ne 0 ]; then
    echo "check-ctoffset: a wire-cell job failed, see $outdir/*.out" >&2
    exit 1
fi

analyze() {
    local ret=0
    for drift in bottom top; do
        python3 "$here/check-ctoffset.py" --max-bias 0.2 --max-spread 0.1 $drift \
                "$outdir/$drift-ctformula-noisefalse.tar.bz2" \
                "$outdir/$drift-ctformula-noisetrue.tar.bz2" || ret=1
        # The wrong (PD-VD data) offset must show the expected ~+1.4 cm bias:
        # require it to fail a 1.0 cm bias cut.
        echo "-- ctoffset=4us ($drift): the bias cut below is EXPECTED to FAIL"
        if python3 "$here/check-ctoffset.py" --max-bias 1.0 $drift \
                   "$outdir/$drift-ct4-noisefalse.tar.bz2" \
                   "$outdir/$drift-ct4-noisetrue.tar.bz2"; then
            echo "check-ctoffset: ctoffset=4us unexpectedly closes for $drift drift"
            ret=1
        fi
    done
    if [ $ret -eq 0 ]; then echo "check-ctoffset: PASS"; else echo "check-ctoffset: FAIL"; fi
    return $ret
}
analyze > "$outdir/results.txt" 2>&1
status=$?
cat "$outdir/results.txt"
exit $status
