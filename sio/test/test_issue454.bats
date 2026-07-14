#!/usr/bin/env bats

# Test for https://github.com/WireCell/wire-cell-toolkit/issues/454

bats_load_library wct-bats.sh

@test "round trip depos npz to tar to tar" {
    cd_tmp file

    local orignpz="$(resolve_file muon-depos.npz)"
    [[ -f "$orignpz" ]]

    local cfg="$(relative_path check-dfs2.jsonnet)"
    local log1=npz-to-tar.log
    local log2=tar-to-npz.log

    run_idempotently -s "$orignpz" -s "$cfg" -t "muon-depos.tar" -t "$log1" -- \
                     wire-cell -l "$log1" -L debug \
                     -c "$cfg" -A tarin="$orignpz" -A tarout="muon-depos.tar"
    
    [[ -n "$(grep '<DepoFileSource:> EOS at call=1' $log1)" ]]
    [[ -n "$(grep '<DepoFileSink:> EOS at call=1' $log1)" ]]

    run_idempotently -s "$muon-depos.tar" -s "$cfg" -t "muon-depos2.npz" -t "$log2" -- \
                     wire-cell -l "$log2" -L debug \
                     -c "$cfg" -A tarin="muon-depos.tar" -A tarout="muon-depos2.npz"
    
    # EOS at call=0 on failure
    [[ -n "$(grep '<DepoFileSource:> EOS at call=1' $log2)" ]]
    [[ -n "$(grep '<DepoFileSink:> EOS at call=1' $log2)" ]]

    
}
