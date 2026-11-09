#!/usr/bin/env bats

# Run the SNB TA study job for one frame of 39Ar background (no MARLEY
# needed) and check its truth tables with the metrics script.

# bats file_tags=trig

bats_load_library wct-bats.sh

setup_file () {
    cd_tmp
    local cfg="$(srcdir trig)/cfg/snb-ta-study.jsonnet"
    export WIRECELL_PATH="$(srcdir trig)/cfg:$(srcdir)/cfg${WIRECELL_PATH:+:$WIRECELL_PATH}"
    run_idempotently -s "$cfg" -t snb.npz -t wire-cell.log -- \
        command wire-cell -l wire-cell.log -L debug -c "$cfg" \
        -A outfile=snb.npz -A nframes=1 -A radiologicals=Ar39
}

@test "trig snb ta study makes truth tables" {
    cd_tmp file
    [[ -s snb.npz ]]
    run grep -c "TATruth:> frame 0" wire-cell.log
    [[ "$output" == "1" ]]
}

@test "trig snb ta metrics runs on them" {
    cd_tmp file
    run python3 "$(srcdir trig)/pytools/snb_ta_metrics.py" snb.npz --json metrics.json
    echo "$output"
    [[ "$status" -eq 0 ]]
    [[ "$output" =~ "Ar39" ]]
    # 39Ar alone makes many TAs per APA, in the tens of kHz
    run python3 -c "import json; r=json.load(open('metrics.json'))['results'][0]; print(r['rate_background_hz'] > 1e4)"
    [[ "$output" == "True" ]]
}
