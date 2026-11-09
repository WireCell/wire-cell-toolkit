#!/usr/bin/env bats

# Checks of the SNB pointing proxy and the MARLEY summary on synthetic
# inputs (python, numpy only).

# bats file_tags=trig

bats_load_library wct-bats.sh

@test "trig snb pointing proxy" {
    run python3 "$(srcdir trig)/test/check_snbpoint.py"
    echo "$output"
    [[ "$status" -eq 0 ]]
    [[ "$output" =~ "checks pass" ]]
}
