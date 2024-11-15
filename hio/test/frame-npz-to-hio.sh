#!/bin/bash

# Test new hio using just HDF5 C lib.  See:
#  https://github.com/WireCell/wire-cell-toolkit/issues/192
#
# This script requires "adc.npz" and "sig.npz" files in WCT "frame file format"
# (not tdm) to be in PWD.  Also "pip install h5py"

mydir="$(dirname $(realpath $BASH_SOURCE))"

set -e
set -x

# Test generating HDF5
for gzqual in 0 9 # 0-9 possible, 0 is no compression, 9 is most.
do
    for chunk in 0 32 64 128 256 512  # Note, chunk=0 implies gz=0.
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

# Test consuming HDF5
for gzqual in 0 9
do
    for chunk in 0 32 512
    do
        post="${gzqual}-${chunk}"

        for ext in pdf png
        do
            python $mydir/frame-npz-to-hio.py adc-${post}.hdf adc-${post}.${ext} '100/frame_orig0'
            python $mydir/frame-npz-to-hio.py sig-${post}.hdf sig-${post}.${ext} '100/frame_gauss0'
        done
    done
done
