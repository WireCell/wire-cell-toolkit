#!/usr/bin/env bats

tstdir="$(realpath $BATS_TEST_DIRNAME)"
pkgdir="$(dirname $tstdir)"
pkg="$(basename $pkgdir)"
topdir="$(dirname $pkgdir)"
blddir="$topdir/build"
bindir="$blddir/$pkg"

setup_file () {
    export outdir=$(mktemp -d /tmp/wct-util-test-tiling-uboone.XXXXXX)

    run $bindir/check_act2viz -o $outdir/blobs.svg -w microboone-celltree-wires-v2.1.json.bz2 -n 0.01 -d $outdir/blobs.txt $tstdir/activities3.txt 
    echo "$output"
    [ "$status" -eq 0 ]
}

@test "reproduce act2vis blob finding" {

    run diff $tstdir/activities3-act2viz-uboone.txt $outdir/blobs.txt
    echo "$output"
    [ "$status" -eq 0 ]
    [ -z "$output" ] 
}

@test "reproduce full blob finding" {

    run diff $tstdir/activities3-full.txt $outdir/blobs.txt
    echo "$output"
    [ "$status" -eq 0 ]
    [ -z "$output" ] 
}

@test "no missing bounds" {
    run grep -E 'pind:\[([0-9]+),\1\]' $outdir/blobs.txt
    [ -z "$output" ]
}

teardown_file () {
    if [ -n "$outdir" ] ; then
       rm -rf "$outdir"
    fi
}
