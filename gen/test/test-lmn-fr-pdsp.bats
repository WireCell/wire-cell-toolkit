#!/usr/bin/env bats

# Run things for LMN FR resampling test for detector:
DETECTOR=pdsp

# The field sampling period in ns which we want.
FIELDS_TDE_TICK=100

# ADC sampling period in ns that we want.
ADC_TDE_TICK=500

# The field sampling period in ns that we will resample.
FIELDS_BDE_TICK=64

# Original ADC sampling period in ns that we will resample.
ADC_BDE_TICK=512



bats_load_library wct-bats.sh


setup_file () {
    cd_tmp

    # prime idempotent running
    cat <<EOF > depo-line-args-tmp
depo-line \
--electron-density -5000 \
--first -3578.36*mm,76.1*mm,3.3497*mm \
--last -1.5875*mm,6060*mm,2303.37*mm \
--sigma 1*mm,1*mm \
--output depos.npz
EOF
    mv_if_diff depo-line-args-tmp depo-line-args

    cat <<EOF > preplots-args-tmp
lmn-fr-plots  \
--period 100*ns \
--zoom-start 50*us \
--plane 0 \
--detector-name pdsp \
-o lmn-fr-pdsp-0.pdf \
-O lmn-fr-pdsp-0.org
EOF
    mv_if_diff preplots-args-tmp preplots-args
}


@test "lmn fr preplots" {
    cd_tmp file

    read -a args < preplots-args

    run_idempotently -s preplots-args -t lmn-fr-pdsp-0.pdf -t lmn-fr-pdsp-0.org -- \
                     wirecell-resp "${args[@]}"

}


fields_original="$(wirecell-util resolve -k fields "$DETECTOR" | head -1)"
fields_name="$(basename "$fields_original" .json.bz2)"

# These file names MUST match those set by the "resample_*" PDSP variants.
fields_tde="${fields_name}-${FIELDS_TDE_TICK}ns.json.bz2"
fields_bde="${fields_name}-${FIELDS_BDE_TICK}ns.json.bz2"


@test "lmn fr resample fr" {
    cd_tmp file

    # Produce field file names that match the resample_* PDSP variants.

    # default PDSP fr.period is not exactly 100ns so we rewrite it with 100ns
    # forced and also put the result locally.
    run_idempotently -s "$files_original" -t "$fields_tde" -- \
                     wirecell-resp condition \
                     -P "$FIELDS_TDE_TICK"'*ns' \
                     -o "$fields_tde" \
                     "$fields_original"

    # Resample 100ns TDE fields to construct a 64ns BDE fields
    run_idempotently -s "$fields_tde" -t "$fields_bde" -- \
                     wirecell-resp resample  \
                     -t "$FIELDS_BDE_TICK"'*ns' \
                     -o "$fields_bde" \
                     "$fields_tde"
}

@test "lmn fr gen depos" {
    cd_tmp file

    read -a args < depo-line-args

    run_idempotently -s depo-line-args -t depos.npz -- \
                     wirecell-gen "${args[@]}"
}

tier_name () {
    local tier="$1" ; shift
    local tick="$1" ; shift
    local ext="${1:-npz}"
    echo "${DETECTOR}-${tier}-${tick}ns.${ext}"
}

@test "lmn fr sim" {
    cd_tmp file

    # We simulate with at 100ns->500ns and 64ns->512ns.
    for tick in $ADC_TDE_TICK $ADC_BDE_TICK
    do
        log="$(tier_name wct-sim $tick log)"
        out="$(tier_name digits $tick)"
        run_idempotently -s depos.npz -t "$log" -t "$out" -- \
                         wire-cell -l "$log" -L debug \
                         -A input="depos.npz" \
                         --tla-code output="{sim:\"$out\"}" \
                         -A detector="${DETECTOR}" \
                         -A variant="resample_${tick}ns" \
                         -A tasks="drift,sim" \
                         layers/omnijob.jsonnet
    done
}

@test "lmn fr resample digits" {
    cd_tmp file

    # We symlink the TDE digits to a new name to match later pattern
    ln -sf \
       $(tier_name digits $ADC_TDE_TICK) \
       $(tier_name digits-$ADC_TDE_TICK $ADC_TDE_TICK) 
       

    # Resample BDE->TDE.  Note, eventually this should be done inside NF in
    # order to save one FFT round trip per channel.
    digin=$(tier_name digits $ADC_BDE_TICK)
    digout=$(tier_name digits-$ADC_BDE_TICK $ADC_TDE_TICK)
    run_idempotently -s $digin -t $digout -- \
                     wirecell-util resample \
                     -t "$ADC_TDE_TICK"'*ns' \
                     -o $digout \
                     $digin
}

# fixme: should move this into a hear document
# select_channels_near='335:350,1500:1548,2080:2100'
select_channels_near='348,1547,2080'
# select_channels_far='0:20,1200:1220,2540:2560'
select_channels_far='0,1200,2559'

digit_trange_near='200*us,300*us'
digit_trange_far='2500*us,2600*us'

signal_trange_near='180*us,210*us'
signal_trange_far='2450*us,2490*us'

@test "lmn fr adc plots" {
    cd_tmp file

    dig1=$(tier_name digits $ADC_TDE_TICK)
    dig2=$(tier_name digits $ADC_BDE_TICK)
    dig3=$(tier_name digits-$ADC_BDE_TICK $ADC_TDE_TICK)

    out=channels-adc-near.pdf
    run_idempotently -s "$dig1" -s "$dig2" -s "$dig3" -t "$out" -- \
                     wirecell-plot channels \
                     -o "$out" \
                     --trange "$digit_trange_near" \
                     --channel "$select_channels_near" \
                     "$dig1" "$dig2" "$dig3"
                  
    out=channels-adc-far.pdf
    run_idempotently -s "$dig1" -s "$dig2" -s "$dig3" -t "$out" -- \
                     wirecell-plot channels \
                     -o "$out" \
                     --trange "$digit_trange_far" \
                     --channel "$select_channels_far" \
                     "$dig1" "$dig2" "$dig3"
                  
}


@test "lmn fr sigproc" {
    cd_tmp file

    local name=pdsp

    # Do normal signal processing using normal tick on normally sampled and
    # resampled digits
    for tick in $ADC_TDE_TICK $ADC_BDE_TICK 
    do
        log="$(tier_name wct-sp-$tick $ADC_TDE_TICK log)"
        digin=$(tier_name digits-$tick $ADC_TDE_TICK)
        sigout=$(tier_name signals-$tick $ADC_TDE_TICK)
        run_idempotently -s "$digin" -t "$log" -t "$sigout" -- \
                         wire-cell -l "$log" -L debug \
                         -A input="$digin" \
                         --tla-code output='{sp:"'$sigout'"}' \
                         -A detector="${DETECTOR}" \
                         -A variant="resample_${ADC_TDE_TICK}ns" \
                         -A tasks="sp" \
                         layers/omnijob.jsonnet
    done

}

@test "lmn fr sig plots" {
    cd_tmp file

    sig1=$(tier_name signals-$ADC_TDE_TICK $ADC_TDE_TICK)
    sig2=$(tier_name signals-$ADC_BDE_TICK $ADC_TDE_TICK)

    out=channels-sig-near.pdf
    run_idempotently -s "$sig1" -s "$sig2" -t "$out" -- \
                     wirecell-plot channels \
                     -o "$out" \
                     --trange "$signal_trange_near" \
                     --channel "$select_channels_near" \
                     "$sig1" "$sig2"

    out=channels-sig-far.pdf
    run_idempotently -s "$sig1" -s "$sig2" -t "$out" -- \
                     wirecell-plot channels \
                     -o "$out" \
                     --trange "$signal_trange_far" \
                     --channel "$select_channels_far" \
                     "$sig1" "$sig2"
                  
}
