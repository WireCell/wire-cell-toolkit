#!/usr/bin/env bats

# fixme: this is first developed outside of check-test branch.  It
# should be refactored to use that Bats support, once merged.
# It does assume a new bats such as the one provided by check-test.

#   coherent/incoherent
#   vin = voltage level
#   adc = ADC counts
#   dac = ADC scaled to voltage
#   spectra output of NoiseModeler
# test-noise-roundtrip-inco-vin.npz
# test-noise-roundtrip-inco-dac.npz
# test-noise-roundtrip-inco-adc.npz
# test-noise-roundtrip-inco-spectra.json.bz2
# test-noise-roundtrip-cohe-vin.npz
# test-noise-roundtrip-cohe-dac.npz
# test-noise-roundtrip-cohe-adc.npz
# test-noise-roundtrip-cohe-spectra.json.bz2

# The bug 202 used "normal" distribution at mean=0.04 instead of
# mean=0.0.  
#
# The way noise generaion is supposed to work: For each frequency bin
# we draw 2 normal (mean=0, sigma=1) random numbers and multiply each
# by the input real-valued "mode" (not technically "mean") spectral
# amplitude at that frequency.  We then interpret the two numbers as
# the real and imaginary parts of the sampled spectrum at that
# frequency.  This is equivalent to a Rayleigh-sampled radius and a
# uniform-sampled phase angle.
# 
# The bug shifted the two sampled "normals" from being distributed
# around mean=0 to being distributed around mean=0.04.  This shifts
# the center of the complex distribution from 0 to 0.04*(1+i).  Thus
# the "mode" of the sampled (real part) distribution should be
# increased by the bug.  Right?  What is the mean distance from origin
# to a circle that is offset from the origin by less than one radius?


declare -A bugmean=( [nobug]=0.0 [smlbug]=0.04 [bigbug]=0.08 )

setup_file () {
    # cd_tmp

    local cfgfile="${BATS_TEST_FILENAME%.bats}.jsonnet"
    local base="$(basename ${BATS_TEST_FILENAME} .bats)"

    for bugsiz in nobug smlbug bigbug
    do
        for trunc in round floor
        do
            wrkdir="${trunc}-${bugsiz}"
            if [ -d "$wrkdir" ] ; then
                echo "already have $wrkdir" 1>&3
                continue
            else
                mkdir -p $wrkdir
                cd $wrkdir
                local bug=${bugmean[$bugsiz]}
                echo "BUG=$bug" 1>&3
                run bash -c "wcsonnet -A bug202=$bug -A round=$trunc $cfgfile > ${base}.json"
                echo "$output"
                [[ "$status" -eq 0 ]]
                if [ "$trunc" = "round" ] ; then
                    [[ "$(grep -c '"round" : true' ${base}.json)" -eq 3 ]]
                else
                    [[ "$(grep -c '"round" : false' ${base}.json)" -eq 3 ]]
                fi
                cd ..
            fi
        done
    done

    for bugsiz in nobug smlbug bigbug
    do
        for trunc in round floor
        do
            wrkdir="${trunc}-${bugsiz}"
            cd $wrkdir
            if [ -f "${base}.log" ] ; then
                echo "Already ran wire-cell for $wrkdir" 1>&3
            else
                run wire-cell -l "${base}.log" -L debug "${base}.json"
                echo "$output"
                [[ "$status" -eq 0 ]]
                if [ "$trunc" = "round" ] ; then
                    [[ -n "$(grep round=1 ${base}.log)" ]]
                else
                    [[ -n "$(grep round=0 ${base}.log)" ]]
                fi
            fi
            cd ..
        done
    done
}

@test "plot adc" {
    # cd_tmp file
    local mydir="$(dirname ${BATS_TEST_FILENAME})"
    local base="$(basename ${BATS_TEST_FILENAME} .bats)"

    for bugsiz in nobug smlbug bigbug
    do
        for trunc in round floor
        do
            wrkdir="${trunc}-${bugsiz}"
            cd $wrkdir
            for name in inco cohe empno
            do
                nbase="${base}-${name}-adc"
                if [ -f "${nbase}.npz" ] ; then
                    python $mydir/check_noise_roundtrip.py plot \
                           -o "${nbase}-plot.png" "${nbase}.npz"
                fi
            done
            cd ..
        done
    done
    # archive these
}


@test "plot spectra" {
    # cd_tmp file
    local mydir="$(dirname ${BATS_TEST_FILENAME})"
    local jsonfile="$(basename ${BATS_TEST_FILENAME} .bats).json"
    local base="$(basename ${BATS_TEST_FILENAME} .bats)"

    for bugsiz in nobug smlbug bigbug
    do
        for trunc in round floor
        do
            wrkdir="${trunc}-${bugsiz}"
            cd $wrkdir
            for name in inco empno
            do
                local nbase="test-noise-roundtrip-${name}-spectra"
                if [ ! -f "${nbase}.json.bz2" ] ; then
                    continue
                fi
                wirecell-sigproc plot-noise-spectra -z \
                                 "${nbase}.json.bz2" \
                                 "${nbase}-output.pdf"
                if [ "$name" = "inco" ] ; then
                    python $mydir/check_noise_roundtrip.py \
                           configured-spectra \
                           -n $name \
                           "$jsonfile" \
                           ${nbase}-input.pdf
                fi
                if [ "$name" = "empno" ] ; then
                    # warning: hard-coding empno spectra out of laziness
                    # and assuming it's the one used in the Jsonnet!
                    wirecell-sigproc plot-noise-spectra -z \
                                     protodune-noise-spectra-v1.json.bz2 \
                                     "${nbase}-input.pdf"
                fi

            done
            cd ..
        done
    done
}

@test "compare round vs floor no bug" {
    # cd_tmp file

    local base="$(basename ${BATS_TEST_FILENAME} .bats)"

    for bl in ac none
    do
        for thing in wave spec
        do
            for grp in $(seq 10)
            do
                local num=$(( $grp - 1 ))
                local chmin=$(( $num * 256 ))
                local chmax=$(( $grp * 256 ))

                wirecell-plot comp1d \
                    --transform $bl \
                    --single \
                    --tier '*' -n $thing --chmin $chmin --chmax $chmax \
                    --output comp1d-floor+round-nobug-${thing}-${bl}-grp${grp}.png \
                    floor-nobug/test-noise-roundtrip-inco-adc.npz \
                    round-nobug/test-noise-roundtrip-inco-adc.npz
            done
        done
    done
}

@test "compare bug and no bug with round" {

    # cd_tmp file
    local base="$(basename ${BATS_TEST_FILENAME} .bats)"

    for bl in ac none
    do
        for thing in wave spec
        do
            for grp in $(seq 10)
            do
                local num=$(( $grp - 1 ))
                local chmin=$(( $num * 256 ))
                local chmax=$(( $grp * 256 ))

                wirecell-plot comp1d \
                    --transform $bl \
                    --single \
                    --tier '*' -n $thing --chmin $chmin --chmax $chmax \
                    --output comp1d-round-bugs-${thing}-${bl}-grp${grp}.png \
                    round-bigbug/test-noise-roundtrip-inco-adc.npz \
                    round-smlbug/test-noise-roundtrip-inco-adc.npz \
                    round-nobug/test-noise-roundtrip-inco-adc.npz
            done
        done
    done
}

@test "compare bug and no bug with floor" {

    # cd_tmp file
    local base="$(basename ${BATS_TEST_FILENAME} .bats)"

    for bl in ac none
    do
        for thing in wave spec
        do
            for grp in $(seq 10)
            do
                local num=$(( $grp - 1 ))
                local chmin=$(( $num * 256 ))
                local chmax=$(( $grp * 256 ))

                wirecell-plot comp1d \
                    --transform $bl \
                    --single \
                    --tier '*' -n $thing --chmin $chmin --chmax $chmax \
                    --output comp1d-floor-bugs-${thing}-${bl}-grp${grp}.png \
                    floor-bigbug/test-noise-roundtrip-inco-adc.npz \
                    floor-smlbug/test-noise-roundtrip-inco-adc.npz \
                    floor-nobug/test-noise-roundtrip-inco-adc.npz 
            done
        done
    done
}

# @test "no weird endpoints of mean waveform" {
#     local mydir="$(dirname ${BATS_TEST_FILENAME})"
#     declare -a line

#     ## good lines are:
#     # U t 4096 0.02282684326171875 0.06892286163071752 4096 1291 185 13 2 0 0 0 0 0 0
#     # U c 800 0.02282684326171875 0.4955734777788893 800 376 0 0 0 0 0 0 0 0 0
#     # V t 4096 0.01374267578125 0.03856538100586203 4096 1334 192 13 1 0 0 0 0 0 0
#     # V c 800 0.01374267578125 0.49734063258678657 800 386 0 0 0 0 0 0 0 0 0
#     # W t 4096 0.009838104248046875 0.016252509189760622 4096 1328 177 14 0 0 0 0 0 0 0
#     # W c 960 0.009838104248046875 0.49116500530978957 960 463 0 0 0 0 0 0 0 0 0
#     ## bad lines are like:
#     # U t 4096 0.23914672851562502 1.838315680187899 4096 30 22 15 15 14 13 13 13 11 10
#     # U c 800 0.239146728515625 0.4333089727023359 800 201 0 0 0 0 0 0 0 0 0
#     # V t 4096 0.24620086669921876 1.0176595165745093 4096 30 22 15 15 14 13 13 13 11 10
#     # V c 800 0.24620086669921876 0.43083546962103825 800 198 0 0 0 0 0 0 0 0 0
#     # W t 4096 0.23283335367838542 0.20440386815418213 4096 39 26 15 14 14 14 13 11 10 10
#     # W c 960 0.23283335367838542 0.4332668683930272 960 249 0 0 0 0 0 0 0 0 0    

#     while read line ; do
#         parts=($line)

#         echo "$line"
#         if [ "${parts[1]}" = "c" ] ; then
#             continue
#         fi
#         # require nothing past 5 sigma
#         [[ "${parts[10]}" = "0" ]]
        
#         # mean should be small.  dumb subst for floating point check.
#         [[ -n "$(echo ${parts[3]} | grep '^0\.0')" ]]

#         # rms should be small.  dumb subst for floating point check.
#         [[ -n "$(echo ${parts[4]} | grep '^0\.0')" ]]
        
#     done < <(python $mydir/check_noise_roundtrip.py stats test-noise-roundtrip-inco-adc.npz)
    
    
# }

