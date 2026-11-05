#!/usr/bin/env bats

# pr_ml doc 54: the truth labeler v2 (truth/python/labels_v2) on synthetic events drawn on the SBND wires (selftest_v2.py).
# Needs the SBND wires file (WIRECELL_PATH or /home/xqian/toolkit-dev/wire-cell-data); skipped when absent.

bats_load_library wct-bats.sh

@test "truth labels_v2 selftests (a)-(f) pass" {
    local wires=""
    for d in ${WIRECELL_PATH//:/ } /home/xqian/toolkit-dev/wire-cell-data; do
        [ -f "$d/sbnd-wires-geometry-v0206.json.bz2" ] && wires="$d/sbnd-wires-geometry-v0206.json.bz2" && break
    done
    [ -n "$wires" ] || skip "no SBND wires file"
    run python3 "$BATS_TEST_DIRNAME/../python/labels_v2/selftest_v2.py" "$wires"
    echo "$output"
    [ "$status" -eq 0 ]
    [[ "$output" == *"selftest_v2: PASS"* ]]
}
