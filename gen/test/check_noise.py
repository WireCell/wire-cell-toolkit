#!/usr/bin/env python

# this holds a few ugly hacks.  nicer stuff in wire-cell-python.

import io
import sys
import numpy
import pathlib
import tarfile
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    tfpath = sys.argv[1]
except IndexError:
    p = pathlib.Path(__file__)
    tfpath = p.parent.parent.parent / "build" / "gen" / "test_noise.tar"
assert(tfpath.exists())


tf = tarfile.open(tfpath)
def get(name):
    dat = tf.extractfile(name + ".npy").read()
    a = io.BytesIO()
    a.write(dat)
    a.seek(0)
    return numpy.load(a)



ss_freq = get("ss_freq")
ss_spec = get("ss_spec")

freqs = get("freqs")
spec = get("true_spectrum")
aspec = get("average_spectrum")

pdfpath = tfpath.parent / (tfpath.stem + ".pdf")
with PdfPages(pdfpath) as pdf:

    fig,ax = plt.subplots(1,1)
    ax.plot(freqs, spec)
    ax.plot(ss_freq, ss_spec)
    pdf.savefig(fig)
    plt.close()

    fig,ax = plt.subplots(1,1)
    for ind in range(5):
        w = get(f'example_wave{ind}')
        plt.plot(w)
    pdf.savefig(fig)
    plt.close()

    fig,ax = plt.subplots(1,1)
    ax.plot(freqs, aspec, label="mean")
    ax.plot(freqs, spec, label="true")
    fig.legend()
    pdf.savefig(fig)
    plt.close()

    
print(str(pdfpath))
