#!/usr/bin/env python
'''
Make a plot with file produced by:

  wire-cell -L debug -l stderr --tla-code tags='["orig0"]' \
    -A infile=frames.npz -A outfile=frames.hdf hio/test/test-frame-npz-to-hio.jsonnet
'''
import sys
import h5py
import numpy
import matplotlib.pyplot as plt

def main(infile, outfile, datapath='100/frame_orig0', *args):
    fp = h5py.File(infile,'r')

    a = fp[datapath][:]  # load and convert to numpy array
    print(f'{a.shape=} {a.dtype=}')
    m = numpy.median(a, axis=1)
    b = (a.T-m).T
    plt.imshow(b, aspect='auto',interpolation="none")
    plt.colorbar()
    plt.xlabel("ticks")
    plt.ylabel("channel index")
    plt.savefig(outfile, dpi=300)

if '__main__' == __name__:
    main(*sys.argv[1:])
