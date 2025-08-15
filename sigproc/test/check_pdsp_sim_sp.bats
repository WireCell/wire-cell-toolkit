#!/usr/bin/env bats

# composited check for sim and sigproc using PDSP

# bats file_tags=sim,sigproc,PDSP,history

bats_load_library wct-bats.sh

# The intention is to run this test in multiple releases and compare across releases.
# bats test_tags=history,plots,implicit
@test "composited check for sim and sigproc using PDSP" {

    cd_tmp

    local name="check_pdsp_sim_sp"
    local outfile="${name}.tar.bz2" # format with support going back the longest
    local cfgfile="${BATS_TEST_FILENAME%.bats}.jsonnet"
    local depofile=( $(input_file depos/cosmic-500-1.npz) )
    yell "depofile: ${depofile}"

    wire-cell -l "$logfile" -L debug \
        -V input=$depofile \
        -V output="$outfile" \
        -c "$cfgfile"
    yell "outfile: ${outfile}"

    saveout -c history "$outfile"

    for what in wave
    do
        local pout="${name}-comp1d-${what}.png"
        check wcpy plot comp1d \
                -o $pout \
                -t 'orig' -n $what \
                --chmin 700 --chmax 701 -s \
                "${outfile}"
        [[ -s "$pout" ]]
        saveout -c plots "$pout"
    done
}

