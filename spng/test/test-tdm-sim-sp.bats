# This runs sim + sp in a few stages using tdm SP


bats_load_library wct-bats.sh

@test "spng sim depos" {
    cd_tmp file

    local cfgdir=$(srcdir spng)/cfg

    local og_depofile=( $(resolve_file test/data/muon-depos.npz) )
    [[ -f $og_depofile ]]

    local depofile=muon-depos-moved.npz
        
    run_idempotently -s $og_depofile -t $depofile -- \
                     wcpy gen shift-depos -c '-2*m,3*m,1*m' $og_depofile $depofile

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

    [ -n "$(grep 'DepoTransform:tpc0.*output: call=0, ndepos_used=7512,' muon.log)" ]
    [ -n "$(grep 'DepoTransform:tpc1.*output: call=0, ndepos_used=0,' muon.log)" ]
    [ -n "$(grep 'DepoTransform:tpc2.*output: call=0, ndepos_used=6793,' muon.log)" ]
    [ -n "$(grep 'DepoTransform:tpc3.*output: call=0, ndepos_used=0,' muon.log)" ]

}

@test "spng sim spng" {
    cd_tmp file
    local cfgdir=$(srcdir spng)/cfg

    local cfg="$cfgdir/spng/test-spng-tpc.jsonnet"

    for tpcid in 0 1 2 3
    do
        local log="test-spng-tpc${tpcid}.log"
        local inf="muon-anode${tpcid}.npz"
        local outf="muon-signals-anode${tpcid}.pkl"

        run_idempotently -s "$cfg" \
                         -t "$inf" \
                         -t "$log" \
                         -t "$outf" \
                         -- \
                         wire-cell -L debug  -l $log \
                         "$cfg" \
                         -A "tpcid=$tpcid"  \
                         -A "input=$inf"  \
                         -A "output=$outf" \
                         -A device=cpu -A engine=Pgrapher

    done
}
