#!/usr/bin/env bats

bats_load_library "wct-bats.sh"

function gen_fname () {
    local tier="$1" ; shift
    local num="$1" ; shift
    local ext="${1:-npz}"
    printf "%s-%s-%02.0f.%s" "pdsp" "$tier" "$num" "$ext"
}

function generate_depos () {
    local num="1"
    local fname="$(gen_fname depos $num)"
    local pname="$(gen_fname params $num json)"

    # APA0, Face 0 is approx at x=[0,-3.6m], y=[0,6m], z=[0,2.6m]
    # local corner='-1*m,3*m,1.0*m'
    local corner='0,0,0'
    # local diagonal='-1*m,0*m,0.6*m'
    local diagonal='-3.6*m,6*m,2.6*m'


    # wcpy gen depo-lines --track-info "$pname" -T 0.0 -c "$corner" -d "$diagonal" -o "$fname" --seed $num
    wcpy gen depo-point -n 50000 -p '-1*m,3*m,1.3*m' -o "$fname"
    touch "$pname"

    # avoid empty zip of 22 bytes
    file_larger_than "$fname" 23
}

@test "sss generate depos" {
    cd_tmp file
    generate_depos 
}


@test "sss run job" {
    cd_tmp file

    local cfg_file="$(relative_path spdir-metric.jsonnet)"
    local num=1
    local depos="$(gen_fname depos "$num")"
    # save out the intermediate drifted depos for input to metric 
    local drifts="$(gen_fname drifts "$num")"
    local splats="$(gen_fname splats "$num")"
    local signals="$(gen_fname signals "$num")"
    local log="$(gen_fname wct "$num" log)"

    # sort of pointless until depo-lines is also run idempotently
    run_idempotently -s "$depos" -t "$drifts" -t "$splats" -t "$signals" -- \
                     wire-cell -c "$cfg_file" \
                     -l "$log" -L debug \
                     -A variant=simple \
                     -A input="$depos" \
                     --tla-code output="{drift:\"$drifts\",splat:\"$splats\",sp:\"$signals\"}" \
                     -A detector=pdsp \
                     -A tasks="drift,splat,sim,sp"

}

@test "sss plots" {
    cd_tmp file
    local num=1

    for thing in splats signals
    do

        wirecell-plot frame -t '*' -n wave \
                      -o "$(gen_fname "$thing" "$num" pdf)" \
                      --single "$(gen_fname "$thing" "$num")"
    done

}
