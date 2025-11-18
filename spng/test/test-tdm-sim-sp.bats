# This runs sim + sp in a few stages using tdm SP


bats_load_library wct-bats.sh

@test "spng sim depos" {
    cd_tmp file

    local cfgdir=$(srcdir spng)/cfg

    local depofile=( $(resolve_file test/data/muon-depos.npz) )
    [[ -f $depofile ]]

    local cfile=spng/test-detsim.jsonnet
    local ifile=$cfgdir/$cfile
    [[ -f $ifile ]]
    local jfile=test-detsim.json

    # Explicitly compile to json.  Normally we let wire-cell do this internally
    run_idempotently -s $ifile -t $jfile -- \
                     wcsonnet -o $jfile -P $cfgdir \
                     -A input_filename=$depofile \
                     -A output_filename_pattern=frames-anode%d.npz \
                     $cfile
    file_larger_than $jfile 0
    [ -z "$(grep Pnode: $jfile)" ]

    run_idempotently -s $depofile -s $jfile -t sim.log \
                     -t frames-anode0.npz \
                     -t frames-anode1.npz \
                     -t frames-anode2.npz \
                     -t frames-anode3.npz \
                     -- \
                     wire-cell -L debug -l sim.log -c $jfile

}
