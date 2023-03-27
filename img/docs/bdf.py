#!/usr/bin/env python
'''
Make plots for the blob-depo-fill.org presentation
'''
from math import sqrt, pi, cos, sin, tan
import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

def subplots(nrows=1, ncols=1, **kwds):
    return plt.subplots(nrows, ncols, tight_layout=True, **kwds)

def gaussian1d(x, mu, sig):
    x = numpy.array(x)
    return 1.0/(sqrt(2*pi)*sig) * numpy.exp(-0.5*((x-mu)/sig)**2)

def plot_sliced_depos(slices, mus, sigs, tps=100):
    '''
    Slices is a linspace of slice times, (mu,sig) describes depos
    '''
    ls = numpy.linspace(slices[0], slices[-1], tps*len(slices))

    fig,ax = subplots()

    for mu,sig in zip(mus, sigs):
        g = gaussian1d(ls, mu, sig)
        ax.plot(ls, g)

    for t in slices:
        ax.axvline(t)

    ax.set_title("Cartoon depos and time slices")
    ax.set_xlabel("Drift dimension [arb]")
    ax.set_ylabel("Depo charge [arb]")

def plot_slices():
    tmin = 0
    tmax = 3
    slices = numpy.linspace(tmin, tmax,7)
    ndepos=6
    mus = numpy.random.uniform(tmin, tmax, ndepos)
    sigs = numpy.random.uniform(0.1, 1, ndepos)
    plot_sliced_depos(slices, mus, sigs)
    plt.savefig("bdf-slices.pdf")

def draw_wires(ax, ang, pitch, off, N=10):
    pang = ang+pi/2.0
    wvec = numpy.array([cos(ang), sin(ang)])
    pvec = numpy.array([cos(pang), sin(pang)])
    for ind in range(N):
        p1 = (pitch*ind + off)*pvec
        p2 = p1 + wvec
        ax.axline(p1, p2)
        p1 = (pitch*(ind-N) + off)*pvec
        p2 = p1 + wvec
        ax.axline(p1, p2)


def plot_blobs():
    angs = (numpy.array([0., 60., 120.]) + 90) * pi / 180.0
    offs = numpy.array([0.5, 0., 0.])
    pitch = 1.0
    fig,ax = subplots()
    lim = 5
    ax.set_xlim((-lim, lim))
    ax.set_ylim((-lim, lim))
    for ang,off in zip(angs,offs):
        draw_wires(ax, ang, pitch, off)

    ndepos = 10
    dx = numpy.random.uniform(-lim, lim, ndepos)
    dy = numpy.random.uniform(-lim, lim, ndepos)
    dr = numpy.random.uniform(0, 0.1*lim, ndepos)

    depos = list()
    dq = list()
    for sigma in [3,2,1]:
        s = [Circle((x,y), r, alpha=0.3) for x,y,r in zip(dx,dy,dr*sigma)]
        depos += s
        dq += [sigma]*len(dr)
        
    pc = PatchCollection(depos, cmap='viridis')
    pc.set_array(dq)
    ax.add_collection(pc)
    ax.set_aspect('equal')
    ax.set_title("Cartoon blobs and wires")
    ax.set_xlabel("Z [arb]")
    ax.set_ylabel("Y [arb]")
    plt.savefig("bdf-blobs.pdf")
    

        
