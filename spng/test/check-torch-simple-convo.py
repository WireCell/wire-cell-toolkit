#!/usr/bin/env python
'''
Make some plots from test-torch-simple-convo.npz
'''

import sys
import numpy
import matplotlib.pyplot as plt
from wirecell.util.plottools import pages


def title(fig, name, arr):
    tit = f'{name} {arr.shape} {arr.dtype}'
    print(tit)
    fig.suptitle(tit)


def plot_image_colorbar(fig, ax, arr, title=""):
    if title:
        ax.set_title(title)
    im = ax.imshow(arr)
    fig.colorbar(im, ax=ax)

def plot_interval(name, arr):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    title(fig, name, arr)
    plot_image_colorbar(fig,ax,arr)


def plot_complex(name, arr):
    fig, axes = plt.subplots(nrows=2, ncols=2)
    title(fig, name, arr)

    plot_image_colorbar(fig, axes[0,0], numpy.absolute(arr), "magnitude")
    plot_image_colorbar(fig, axes[1,0], numpy.angle(arr), "phase")
    plot_image_colorbar(fig, axes[0,1], numpy.real(arr), "real")
    plot_image_colorbar(fig, axes[1,1], numpy.imag(arr), "imaginary")



def plot_all(pdf, **arrs):
    print(arrs.keys())
    with pages(pdf) as out:
        for name, arr in arrs.items():
            plt.clf();
            if arr.dtype == numpy.dtype("float64"):
                plot_interval(name, arr)
            elif arr.dtype == numpy.dtype("complex64"):
                plot_complex(name, arr)
            else:
                print(f'Unknown dtype: {arr.dtype}')
            out.savefig()


def main(npz, pdf):

    with numpy.load(npz) as arrs:
        plot_all(pdf, **arrs)

if '__main__' == __name__:
    main(*sys.argv[1:])
