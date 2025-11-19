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

        
    run_idempotently -s $ifile \
                     -t $jfile \
                     -t muon-anode0.npz \
                     -t muon-anode1.npz \
                     -t muon-anode2.npz \
                     -t muon-anode3.npz \
                     -- \
                     wire-cell \
                     -L debug -l muon.log \
                     -P $cfgdir \
                     -A engine=TbbFlow \
                     -A input=$depofile \
                     -A output=muon-anode%d.npz \
                     spng/test-detsim.jsonnet

    # Per-anode files all about the same size, mostly noise.
    file_larger_than muon-anode0.npz 16000000
    file_larger_than muon-anode1.npz 16000000
    file_larger_than muon-anode2.npz 16000000
    file_larger_than muon-anode3.npz 16000000

    # the muon-depos land in tpc3.
    [ -z "$(grep 'DepoTransform:tpc0.*output: call=0, ndepos_used=8087' muon.log)" ]
    [ -z "$(grep 'DepoTransform:tpc1.*output: call=0, ndepos_used=8087' muon.log)" ]
    [ -z "$(grep 'DepoTransform:tpc2.*output: call=0, ndepos_used=8087' muon.log)" ]
    [ -n "$(grep 'DepoTransform:tpc3.*output: call=0, ndepos_used=8087' muon.log)" ]

}
