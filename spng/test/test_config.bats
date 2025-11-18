#!/usr/bin/env bats

bats_load_library wct-bats.sh

@test "spng tdm configs detsim" {
    cd_tmp file

    local cfgdir=$(srcdir spng)/cfg
    # make sure we find it via relative path
    local cfile=spng/test-detsim.jsonnet
    local ifile=$cfgdir/$cfile
    local ofile=test-detsim.json
    run_idempotently -s $ifile -t $ofile -- \
                     wcsonnet -o $ofile -P $cfgdir \
                     -A input_filename=depos.npz \
                     -A output_filename_pattern=frames-anode%d.npz \
                     $cfile
    file_larger_than $ofile 0

}
