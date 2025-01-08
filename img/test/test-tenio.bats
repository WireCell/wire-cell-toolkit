#!/usr/bin/env bats

# This file is generated.  Any edits may be lost.
# See img/test/test-tenio.org for source to modify.

# standard WCT Bats support library
bats_load_library "wct-bats.sh"

data_ext=npz
# zip:22, gzip:29, tar:10240, tar.gz:122, no-PC-zip:371
empty_data=500

function make_dag () {
    local src=$1; shift
    local tgt=$1; shift
    #declare -a args=( -A "detector=pdsp" -S "anode_iota=[0]" )
    declare -a args=( -A "detector=pdsp" )
    if [ "$src" != "depo" ] ; then
        args+=( -A "infiles=apa-%(anode)s-${src}.${data_ext}" )
    fi
    args+=( -A "outfiles=apa-%(anode)s-${tgt}.${data_ext}")
    local cfg_file="$(relative_path tenio-${src}-${tgt}.jsonnet)"
    run_idempotently -s "$cfg_file" -t "dag-${tgt}.json" -- \
        wcsonnet -o "dag-${tgt}.json" "${args[@]}" "$cfg_file"
    [[ -s dag-${tgt}.json ]]
    run_idempotently -s "dag-${tgt}.json" -t "dag-${tgt}.png" -- \
        dotify_graph "dag-${tgt}.json" "dag-${tgt}.png"
}

@test "compile configuration for depo to adc" {
    cd_tmp file
    make_dag depo adc 
}

@test "compile configuration for adc to acp" {
    cd_tmp file
    make_dag adc acp 
}

@test "compile configuration for adc to sig" {
    cd_tmp file
    make_dag adc sig 
}

@test "compile configuration for sig to img" {
    cd_tmp file
    make_dag sig img 
}

@test "compile configuration for sig to scp" {
    cd_tmp file
    make_dag sig scp 
}

@test "compile configuration for img to ptc" {
    cd_tmp file
    make_dag img ptc 
}

@test "compile configuration for img to brf" {
    cd_tmp file
    make_dag img brf 
}

function run_dag () {
    local src=$1; shift
    local tgt=$1; shift
    run_idempotently -s apa-0-${src}.${data_ext} -t apa-0-${tgt}.${data_ext} -- \
        wire-cell -l dag-${tgt}.log -L debug dag-${tgt}.json
    local warnings=$(grep '\bW\b' dag-${tgt}.log)
    echo "$warnings" 1>&3
    local errors=$(grep '\bE\b' dag-${tgt}.log)
    [[ -z "$errors" ]]
}

@test "run wire-cell stage depo to adc" {
    cd_tmp file
    run_dag depo adc 
}

@test "run wire-cell stage adc to acp" {
    cd_tmp file
    run_dag adc acp 
}

@test "run wire-cell stage adc to sig" {
    cd_tmp file
    run_dag adc sig 
}

@test "run wire-cell stage sig to img" {
    cd_tmp file
    run_dag sig img 
}

@test "run wire-cell stage sig to scp" {
    cd_tmp file
    run_dag sig scp 
}

@test "run wire-cell stage img to ptc" {
    cd_tmp file
    run_dag img ptc 
}

@test "run wire-cell stage img to brf" {
    cd_tmp file
    run_dag img brf 
}

@test "no unexpected empty files" {
    cd_tmp file

    for tier in adc acp sig ; do
        for apa in 0 1 2 3 4 5 ; do
            file_larger_than apa-${apa}-${tier}.${data_ext} ${empty_data}
        done
    done
    for tier in img ptc ; do
        for apa in 0 2 3 5 ; do # some apas are "empty"
            file_larger_than apa-${apa}-${tier}.${data_ext} ${empty_data}
        done
    done
}
