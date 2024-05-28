#!/usr/bin/env bats
# -*- shell-script -*-

# Caution: if this file has .bats extension it is fully generated.
# Any edits are subject to loss.
# Instead, edit spdir-metric.tmpl and run spdir-metric-gen-tests.sh.

bats_load_library "wct-bats.sh"

DETECTOR=pdsp

# The detector variant.  
VARIANT="spdir"

function gen_fname () {
    local tier="$1" ; shift
    local num="$1" ; shift
    local ext="${1:-npz}"
    printf "%s-%s-%02.0f.%s" "$DETECTOR" "$tier" "$num" "$ext"
}


# The MicroBooNE conventional directions are defined in X-axis of figure 30 of
# the the MB SP paper with angle definitions in figure 29.
#
# Two big caveats in the MB convention:
#
# 1) It relies on U/V wire direction being mirror-symmetric about W.
#
# 2) It samples U/V wire direction space differently than W.  Event when the
# "primed" (U/V) angle theta_xz' matches the W angle theta_xz, the theta_y' and
# theta_y differ.  And so it is not trivial to make U/V vs W comparisons.
mb_theta_xz_deg=( 0  1  3  5 10 20 30 45 60 75 80 82 84 89 )
mb_theta_y_deg=( 90 90 90 90 90 90 90 90 90 90 90 90 90 90 )

# In future, perhaps we pick others.  For now, accept MB's
theta_xz_deg=( ${mb_theta_xz_deg[@]} )
theta_y_deg=( ${mb_theta_y_deg[@]} )
ndirs="${#theta_xz_deg[@]}"
directions=$(seq 0 $(( $ndirs - 1 )))

function generate_depos ()
{
    local num="$1" ; shift
    local theta_xz=${theta_xz_deg[$num]}'*deg'
    local theta_y=${theta_y_deg[$num]}'*deg'

    echo "generate_depos: $num $theta_xz $theta_y"

    local fname="$(gen_fname depos $num)"
    local pname="$(gen_fname params $num json)"

    # note: MB only select based on collection plane=2
    wcpy gen detlinegen \
         -d $DETECTOR --apa 0 --plane 2 --angle-coords=wire-plane \
         --theta_xz="$theta_xz" --theta_y="$theta_y" \
         -o "$fname" -m "$pname"

    # avoid empty zip of 22 bytes
    file_larger_than "$fname" 23
}


@test "spdir metric pdsp generate depos" {
    cd_tmp file

    for num in $directions
    do
        generate_depos $num
    done
}


@test "spdir metric pdsp sim+nf+sp one job graph viz" {
    cd_tmp file
    local cfg_file="$(relative_path spdir-metric.jsonnet)"
    local num="${directions[0]}"
    local depos="$(gen_fname depos "$num")"
    # save out the intermediate drifted depos for input to metric 
    local drifts="$(gen_fname drifts "$num")"
    local splats="$(gen_fname splats "$num")"
    local signals="$(gen_fname signals "$num")"
    local pdf="$(gen_fname graph "$num" pdf)"
    local log="$(gen_fname wct "$num" log)"
    wirecell-pgraph dotify \
                    -A input="$depos" \
                    -A output="{drift:\"$drifts\",splat:\"$splats\",sp:\"$signals\"}" \
                    -A detector=$DETECTOR \
                    -A variant="$VARIANT" \
                    "$cfg_file" "$pdf"
}

@test "spdir metric pdsp sim+nf+sp run job" {
    cd_tmp file
    local cfg_file="$(relative_path spdir-metric.jsonnet)"
    for num in ${directions[*]}
    do
        local depos="$(gen_fname depos "$num")"
        # save out the intermediate drifted depos for input to metric 
        local drifts="$(gen_fname drifts "$num")"
        local splats="$(gen_fname splats "$num")"
        local signals="$(gen_fname signals "$num")"
        local log="$(gen_fname wct "$num" log)"
        wire-cell -c "$cfg_file" \
                  -l "$log" -L debug \
                  -A input="$depos" \
                  --tla-code output="{drift:\"$drifts\",splat:\"$splats\",sp:\"$signals\"}" \
                  -A detector=$DETECTOR \
                  -A variant="$VARIANT"
                    
    done
}


# FIXME
#
# We do not yet have a command to generate a metric from each frame file.
# When we do, fix this function.  For now, it does some bogosity.
#
# This function should consume params.json
function calcualte_metrics ()
{
    local jname="$(gen_fname metrics $num json)"
    cat <<EOF > $jname
{
  "num": $num
  "theta": $num,
  "phi": $num,
  "nchan": 0,
  "sigtot": 0,
  "eff": [],
  "bias": [],
  "res": []
}
EOF
    # This metrics.json should produce rich but summary data.  Eg, the
    # eff/bias/res should be arrays the hit channels.
}
# FIXME

@test "spdir metric pdsp metrics" {
    cd_tmp file
    for num in ${directions[*]}
    do
        calcualte_metrics "$num"
    done
    
}


# FIXME: we also need a canonical plotter that consumes all metric.json.  
@test "spdir metric pdsp plots" {
    cd_tmp file
    for num in ${directions[*]}
    do
        ls -l "$(gen_fname metrics $num json)" >> "$(gen_fname plots all txt)"
    done
}
