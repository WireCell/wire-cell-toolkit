#!/bin/bash

# Test new hio using just HDF5 C lib.  See:
#  https://github.com/WireCell/wire-cell-toolkit/issues/192
#
# This script requires "adc.npz" and "sig.npz" files in WCT "frame file format"
# (not tdm) to be in PWD.

mydir="$(dirname $(realpath $BASH_SOURCE))"

set -e

for gzqual in 1 2 3 4 # 0 5 9
do
    for chunk in 32 64 128 256 512
    do
        post="${gzqual}-${chunk}"

        time wire-cell -L debug -l adc-${post}.log --tla-code tags='["orig0"]' \
                  -A digitize=yes -A infile=adc.npz -A outfile=adc-${post}.hdf \
                  -A chunk=$chunk -A gzip=$gzqual \
                  $mydir/test-frame-npz-to-hio.jsonnet
        wire-cell -L debug -l sig-${post}.log --tla-code tags='["gauss0"]' \
                  -A digitize=no -A infile=sig.npz -A outfile=sig-${post}.hdf \
                  -A chunk=$chunk -A gzip=$gzqual \
                  $mydir/test-frame-npz-to-hio.jsonnet

    done
done

for gzqual in 0 9
do
    for chunk in 32 512
    do

        for ext in pdf png
        do
            python $mydir/frame-npz-to-hio.py adc-${post}.hdf adc-${post}.${ext} '100/frame_orig0'
            python $mydir/frame-npz-to-hio.py sig-${post}.hdf sig-${post}.${ext} '100/frame_gauss0'
        done
    done
done
