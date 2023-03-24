#!/usr/bin/env bats

# Initial test of generating files for input to possible GNNs for
# charge solving and deghosting tasks.

# Note: this test may break for a while as we consolidate the testing
# system and test data repo.

##############
# Fixme: these are temporary replacements for what are provided in
# wct-bats.sh, once we merge check-test we should remove them.
function top () {
    dirname $(dirname $(dirname $(realpath $BATS_TEST_FILENAME)))
}    
function tname () {
    basename "$BATS_TEST_FILENAME" .bats
}
function tdir () {
    dirname "$BATS_TEST_FILENAME"
}
function resolve_path () {
    local want="$1"; shift
    if [[ "$want" =~ ^/.* ]] ; then
        echo $want
    fi
    for pathlst in "$(tdir)" $@ ; do
        for maybe in $(echo ${pathlst} | tr ":" "\n")
        do
            if [ -f "${maybe}/$want" ] ; then
                echo "Found: ${maybe}/$want" 1>&3
                echo "${maybe}/$want"
                return
            else
                echo "Not found: ${maybe}/$want" 1>&3                
            fi
        done
    done
}
function saveout_path () {
    src="$1" ; shift
    tgt="$1"
    name="$(tname)"
    if [ -z "$tgt" ] ; then
        tgt="$(basename $src)"
    fi
    echo "$(top)/build/output/${name}/${tgt}"
}
function saveout () {
    src="$1" ; shift
    tgt="$(saveout_path $src $1)"
    mkdir -p "$(dirname $tgt)"
    cp "$src" "$tgt"
}
function cd_tmp () {

    pwd
    if [ -n "$WCTEST_TMPDIR" ] ; then
        mkdir -p "$WCTEST_TMPDIR"
        cd "$WCTEST_TMPDIR"
        return
    fi

    local which="${1:-test}"
    case $which in
        suite) cd "$BATS_SUITE_TMPDIR";;
        file) cd "$BATS_FILE_TMPDIR";;
        run) cd "$BATS_RUN_TMPDIR";;
        *) cd "$BATS_TEST_TMPDIR";; # "test"
    esac
}
##############

bimg="blobs-img.npz"
btru="blobs-tru.npz"

function setup_file () {

    echo "me= $BATS_TEST_FILENAME"
    echo "top=$(top)"
    # fixme: replace this with WCTEST_{IN,OUT}PUT to allow for variant
    # test.
    local depos="$(resolve_path test/data/muon-depos.npz $(top))"
    [[ -n "$depos" ]]
    [[ -s "$depos" ]]

    local log="wire-cell.log"
    local cfg="$(resolve_path test-uboone-img.jsonnet)"
    [[ -s "$cfg" ]] # temporary test of internal resolve_path
    
    cd_tmp file

    local cfgjson=cfg.json
    local cfgpdf=graph.pdf

    local tlas="-A depos=$depos -A outimg=$bimg -A outtru=$btru"
    ### anything we should test for these intermediate outputs?
    # tlas = "$tlas -A drifted=drifted.npz -A outadc=adc.npz -A outsig=sig.npz"
    if [ -f $cfgjson ] ; then
        echo "Reusing $cfgjson"
    else
        run wcsonnet -o $cfgjson $tlas $cfg
        echo "$output"
        [[ "$status" -eq 0 ]]
    fi
    [[ -s $cfgjson ]]

    if [ -f $cfgpdf ] ; then
        echo "Reusing $cfgpdf"
    else
        run wirecell-pgraph dotify $cfgjson $cfgpdf
        echo "$output"
        [[ "$status" -eq 0 ]]
    fi
    [[ -s $cfgpdf ]]
    saveout $cfgpdf

    if [ -f $log ] ; then
        echo "Reusing $log"
    else
        cmd="wire-cell -l $log -L debug -c $cfgjson"
        echo $cmd
        run $cmd
        echo "$output"
        [[ "$status" -eq 0 ]]
    fi
    saveout $log
    # saveout blobs-img.npz
    # saveout blobs-tru.npz
    [[ -s $bimg ]]
    [[ -s $btru ]]
}                              

@test "same blobs in zip" {
    cd_tmp file

    run bash -c "unzip -v $bimg | grep cluster_ | awk '{print $8}'"
    echo "$output"
    [[ "$status" -eq 0 ]]
    local limg="$output"

    run bash -c "$unzip -v $btru | grep cluster_ | awk '{print $8}'"
    echo "$output"
    [[ "$status" -eq 0 ]]
    local ltru="$output"

    [[ "$(limg)" = "$(ltru)" ]]
}

@test "same blobs in npz" {
    cd_tmp file

    run wirecell-util ls $bimg
    echo "$output"
    [[ "$status" -eq 0 ]]
    local limg="$output"

    run wirecell-util ls $btru
    echo "$output"
    [[ "$status" -eq 0 ]]
    local ltru="$output"

    [[ "$(limg)" = "$(ltru)" ]]
}
