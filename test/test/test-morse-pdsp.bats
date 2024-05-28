#!/usr/bin/env bats

bats_load_library "wct-bats.sh"

setup_file () {
    cd_tmp
    echo 'depo-morse --wires pdsp --output depos.npz' > depo-morse-args-tmp
    mv_if_diff depo-morse-args-tmp depo-morse-args 
}

@test "morse gen depos" {
    cd_tmp file
    read -a args < depo-morse-args

    run_idempotently -s depo-morse-args -t depos.npz -- \
                     wirecell-gen "${args[@]}"
    file_larger_than depos.npz 23
}

@test "morse run wct" {
    cd_tmp file
    local cfg_file="$(relative_path spdir-metric.jsonnet)"

    run_idempotently -s "$cfg_file" -s depos.npz -t wct.log -t drift.npz -t splat.npz -t digit.npz -t signal.npz -- \
                     wire-cell \
                     -l wct.log -L debug \
                     -A input="depos.npz" \
                     --tla-code output="{drift:\"drift.npz\",splat:\"splat.npz\",sim:\"digit.npz\",sp:\"signal.npz\"}" \
                     -A detector=pdsp \
                     -A variant=morse_nominal \
                     -A tasks="drift,splat,sim,sp" \
                     "$cfg_file"
}

@test "morse make plots" {
    cd_tmp file

    run_idempotently -s signal.npz -t peaks.json -- \
                     wirecell-gen morse-summary -o peaks.json signal.npz

    run_idempotently -s peaks.json -t splat-params.json -- \
                     wirecell-gen morse-splat -o splat-params.json peaks.json

    run_idempotently -s signal.npz -t peaks.pdf -- \
                     wirecell-gen morse-plots -o peaks.pdf signal.npz

}
