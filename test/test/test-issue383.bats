#!/usr/bin/env bats

# Sim signal has doubled in master after 0.29.0. 

# for debugging, run like:
#
# $ WCT_BATS_TMPDIR=$(pwd)/junk \
#   WCT_BATS_LOG_SINK=terminal \
#   WCT_BATS_LOG_LEVEL=debug \
#     bats test/test/test-issue383.bats

bats_load_library wct-bats.sh


@test "get depos file" {
    cd_tmp file
    local depourl="https://www.phy.bnl.gov/~bviren/tmp/wcttest/data_repo/ad-hoc/depos.tar.bz2"
    local depofile=$(download_file "$depourl")
    debug "Using depo input file: $depofile"

    local cfg_file="$(relative_path issue383.jsonnet)"
    debug "Using config input file: $cfg_file"

    local refver="0.29.3"
    local refurl="https://www.phy.bnl.gov/~bviren/tmp/wcttest/data_repo/ad-hoc/issue383-${refver}.tar.bz2"
    local reffile=$(download_file "$refurl")
    debug "Reference version: $refver using file: $reffile"

    local this="issue383"
    local thisfile="$(realpath "${this}.tar.gz")"

    run_idempotently -s "$cfg_file" -s "$depofile" -t "$thisfile" -t "log" -- \
                     wire-cell -l log -L debug -c "$cfg_file" -V input=$depofile -V output="$thisfile"

    export MPLBACKEND=Agg

    run_idempotently -s "$thisfile" -t "${this}.pdf" -- \
                     wirecell-plot frame -n wave -o "${this}.pdf" "$thisfile"
    run_idempotently -s "$reffile" -t "${refver}.pdf" -- \
                     wirecell-plot frame -n wave -o "${refver}.pdf" "$reffile"

    local combo="${this}-${refver}"
    run_idempotently -s $thisfile -t "ch700-${combo}.pdf" -- \
                     wirecell-plot comp1d -n wave -o "ch700-${combo}.pdf" --chmin 700 --chmax 701 "$reffile" "$thisfile"
    run_idempotently -s $thisfile -t "ch1230-${combo}.pdf" -- \
                     wirecell-plot comp1d -n wave -o "ch1230-${combo}.pdf" --chmin 1230 --chmax 1231 "$reffile" "$thisfile"
}
