#!/usr/bin/env bats

# pr_ml doc 33 perf round 1: steiner_truth.py's vectorized paths equal the legacy forms they replace (truth/test/steiner_truth_equiv.py).

bats_load_library wct-bats.sh

@test "truth steiner_truth vectorized paths equal the legacy forms" {
    run python3 "$BATS_TEST_DIRNAME/steiner_truth_equiv.py"
    echo "$output"
    [ "$status" -eq 0 ]
    [[ "$output" != *FAIL* ]]
}
