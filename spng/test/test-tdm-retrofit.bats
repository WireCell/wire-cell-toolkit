#!/usr/bin/env bats

# This runs tests covering the "retrofitting" of the SPNG TDM.

bats_load_library wct-bats.sh


run_fans () {
    
    local anode="$1" ; shift
    local paradigm="$1" ; shift

    local depos="$(resolve_file test/data/muon-depos.npz)"
    local cfg="$(relative_path test-frame-fan-outin.jsonnet)"
    local name="fans-anode${anode}-${paradigm}"
    local log="$name.log"
    local out="$name.npz"


    run_idempotently -s $cfg -s $depos -t $log -t $out -- \
                     wire-cell -l $log -L debug \
                       -A anodeid=$anode -A input=$depos -A output=$out -A paradigm=$paradigm $cfg

}

# APA 0 does not get any activity
@test "spng fans empty orig" {
    cd_tmp
    run_fans 0 orig
}
@test "spng fans empty tdm" {
    cd_tmp
    run_fans 0 tdm
}
# APA 3 does.
@test "spng fans depos orig" {
    cd_tmp
    run_fans 3 orig
}
@test "spng fans depos tdm" {
    cd_tmp
    run_fans 3 tdm
}
