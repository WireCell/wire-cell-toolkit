#!/usr/bin/env python
import numpy
import numpy.ma as ma

import scipy
import matplotlib.pyplot as plt
import click
from pathlib import Path
from wirecell import units
from wirecell.util.cli import log, context
from wirecell.util.plottools import pages
from wirecell.util.peaks import (
    baseline_noise, select_activity,
    find1d as find_peaks,
    gauss as gauss_func
)

@context("ssss-pdsp")
def cli(ctx):
    '''
    The ssss-pdsp test
    '''
    pass


@cli.command("log-summary")
@click.argument("log")
def depo_summary(log):
    '''
    Print summary of wct log
    '''
    want = "<AnodePlane:0> face:0 with 3 planes and sensvol: "
    for line in open(log).readlines():
        if want in line:
            n = line.find(": [(") + len(want)
            off = len(want) + line.find(want) + 1
            parts = [p.strip() for p in line[off:-2].split(" --> ")]
            p1 = numpy.array(list(map(float, parts[0][1:-1].split(" "))))/units.mm
            p2 = numpy.array(list(map(float, parts[1][1:-1].split(" "))))/units.mm
            print(f'''
#+caption: Sensitive volume diagonal and ideal line track endpoints.
#+name: tab:diagonal-endpoints
| x (mm) | y (mm) | z (mm) |
|--------|--------|--------|
| {p1[0]}| {p1[1]}| {p1[2]}|
| {p2[0]}| {p2[1]}| {p2[2]}|
            ''')
            return

def load_depos(fname):
    fp = numpy.load(fname)
    d = fp["depo_data_0"]
    i = fp["depo_info_0"]
    ind = i[:,2] == 0
    return d[ind,:]

@cli.command("depo-summary")
@click.argument("depos")
@click.argument("drift")
def depo_summary(depos, drift):
    '''
    Print summary of depos
    '''
    d1 = load_depos(depos)
    d2 = load_depos(drift)
    das = [d1,d2]

    tmins = numpy.array([numpy.min(d[:,0]) for d in das])/units.us
    # print(f'{tmins=}')

    tmaxs = numpy.array([numpy.max(d[:,0]) for d in das])/units.us
    # print(f'{tmaxs=}')

    xmins = numpy.array([numpy.min(d[:,2]) for d in das])/units.mm
    xmaxs = numpy.array([numpy.max(d[:,2]) for d in das])/units.mm

    # [t/x][min/max][depos/drift]

    # copy-paste into org and hit tab
    print(f'''
#+caption: Depo t/x ranges:
#+name: tab:depo-tx-ranges
|         |  min | max | units |
|---------|------|-----|-------|
| depos t | {tmins[0]:.2f}|{tmaxs[0]:.3f} | us | 
| drift t | {tmins[1]:.2f}|{tmaxs[1]:.3f} | us | 
| depos x | {xmins[0]:.2f}|{xmaxs[0]:.2f} | mm | 
| drift x | {xmins[1]:.2f}|{xmaxs[1]:.2f} | mm | 
    ''')

def load_frame(fname, tunit=units.us):
    '''
    Load a frame with time values in explicit units.
    '''
    fp = numpy.load(fname)
    f = fp["frame_*_0"]
    t = fp["tickinfo_*_0"]
    c = fp["channels_*_0"]

    c2 = numpy.array(c)
    numpy.sort(c2)
    assert numpy.all(c == c2)

    cmin = numpy.min(c)
    cmax = numpy.max(c)
    nch = cmax-cmin+1
    ff = numpy.zeros((nch, f.shape[1]), dtype=f.dtype)
    for irow, ch in enumerate(c):
        ff[cmin-ch] = f[irow]

    t0 = t[0]/tunit
    tick = t[1]/tunit
    tf = t0 + f.shape[1]*tick

    # edge values
    extent=(t0, tf+tick, cmin, cmax+1)
    origin = "lower"            # lower flips putting [0,0] at bottom
    ff = numpy.flip(ff, axis=0)
    return (ff, extent, origin)

@cli.command("frame-summary")
@click.argument("files", nargs=-1)
def frame_summary(files):
    data = dict()
    for fname in files:
        p = Path(fname)
        data[p.stem] = load_frame(fname)


    rowpat = '| {name} | {start} | {duration:.0f} | {npos} | {vmin:.1f} | {vmax:.1f} |'
    rows = [
        '#+caption: Frame times in microsecond.',
        '#+name: tab:frame times',
        '| frame | start | duration | n>0 | min | max |',
    ]
    for name, (f,e,o) in sorted(data.items()):
        p = dict(
            name = name,
            start = e[0],
            duration = e[1]-e[0],
            npos = numpy.sum(f > 0),
            vmin = numpy.min(f),
            vmax = numpy.max(f),
        )
        row = rowpat.format(**p)
        rows.append(row)

    print('\n'.join(rows))
    
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

def dump_frame(name,f):
    print(name,f.shape,numpy.sum(f))

def plot_ticks(ax, t0, tf, drift, f1, f2, channel_ranges, speed = 1.565*units.mm/units.us):
    nticks = f1.shape[1]
    times = numpy.linspace(t0,tf,nticks+1, endpoint=True)
    deltat = times[1]-times[0]
    print(f'time binning: [{t0:.1f},{tf:1f}] us x {nticks}')
    dump_frame('splat',f1)
    dump_frame('signal',f2)
    
    print(f'warning: using speed of {speed/(units.mm/units.us)} mm/us')
    ts = (drift[:,0] + 10*units.cm / speed)/units.us
    qs = -1*drift[:,1]

    sigs = drift[:,5] / speed / units.us
    nsig=5
    nspread = numpy.array((nsig/deltat)*sigs, dtype=int)

    gt = numpy.zeros(nticks)
    for t,q,n in zip(ts,qs,nspread):
        tbin = int((t-t0)/deltat)
        gt[tbin-n:tbin+n+1] += q/(2*n+1)

    # h = numpy.histogram(ts, bins=times, weights=qs)
    # ax.plot(h[1][:-1], h[0], label='depo charge')

    for p,chans in zip("UVW",channel_ranges):
        val1 = f1[chans,:].sum(axis=0)
        val2 = f2[chans,:].sum(axis=0)
        ax.plot(times[:-1], val1, label=p+' splat')
        ax.plot(times[:-1], val2, label=p+' signal')
    val1 = f1.sum(axis=0)
    val2 = f2.sum(axis=0)
    ax.plot(times[:-1], val2/3.0, label='total signal / 3')
    ax.plot(times[:-1], gt, label='depo charge')
    ax.plot(times[:-1], val1/3.0, label='total splat / 3')



    # ax.set_xlim((1400,1600))
    # ax.set_xlim((0,200))
    ax.legend(loc='right')



def plot_frame(gs, f, e, o="lower", channel_ranges=None, which="splat", tit=""):
    t0,tf,c0,cf=e

    gs = GridSpecFromSubplotSpec(2,2, subplot_spec=gs,
                  height_ratios = [5,1], width_ratios = [6,1])                                 


    fax = plt.subplot(gs[0,0])
    tax = plt.subplot(gs[1,0], sharex=fax)
    cax = plt.subplot(gs[0,1], sharey=fax)

    cax.set_xlabel(which)
    fax.set_ylabel("channel")
    if which=="signal":
        tax.set_xlabel("time [us]")

    if tit:
        plt.title(tit)
    plt.setp(fax.get_xticklabels(), visible=False)
    plt.setp(cax.get_yticklabels(), visible=False)
    if which=="splat":
        plt.setp(tax.get_xticklabels(), visible=False)

    im = fax.imshow(f, extent=e, origin=o, aspect='auto', vmax=500, cmap='hot_r')

    tval = f.sum(axis=0)
    t = numpy.linspace(e[0],e[1],f.shape[1]+1,endpoint=True)
    tax.plot(t[:-1], tval)
    if channel_ranges:
        for p,chans in zip("UVW",channel_ranges):
            val = f[chans,:].sum(axis=0)
            c1 = chans.start
            c2 = chans.stop
            print(c1,c2,numpy.sum(val))
            tax.plot(t[:-1], val, label=p)
            fax.plot([t0,tf], [c1,c1])
            fax.text(t0 + 0.1*(tf-t0), c1 + 0.5*(c2-c1), p)
        fax.plot([t0,tf], [c2-1,c2-1])
        tax.legend()
    
    cval = f.sum(axis=1)
    c = numpy.linspace(e[2],e[3],f.shape[0]+1,endpoint=True)
    cax.plot(cval, c[:-1])

    return im

def plot_zoom(f1, f2, ts, cs, tit, o="lower"):

    e = (ts[0]*0.5, ts[1]*0.5, cs[0], cs[1])
    ts = slice(*ts)
    cs = slice(*cs)

    f1 = f1[cs, ts]
    f2 = f2[cs, ts]

    fig, axes = plt.subplots(3,1, sharex=True, sharey=True)
    plt.suptitle(tit)

    vmax=2000
    args=dict(extent=e, origin=o, aspect='auto', vmin=-vmax, vmax=vmax, cmap='Spectral')

    im0 = axes[0].imshow(f1, **args)
    im1 = axes[1].imshow(f2, **args)
    im2 = axes[2].imshow(f1-f2, **args)
    axes[2].set_xlabel('time [us]')
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    fig.colorbar(im1, cax=cbar_ax)

def smear_tick(f1, sigma=3.0, nsigma=6):
    print(f'smear: tot sig {numpy.sum(f1)}')
    gx = numpy.arange(-sigma*nsigma, sigma*nsigma, 1)
    if gx.size%2:
        gx = gx[:-1]
    gauss = numpy.exp(-0.5*(gx/sigma)**2)
    gauss /= numpy.sum(gauss)
    smear = numpy.zeros_like(f1[0])
    hngx = gx.size//2
    smear[:hngx] = gauss[hngx:]
    smear[-hngx:] = gauss[:hngx]
    nticks = f1[0].size
    totals = numpy.sum(f1, axis=1)
    for irow in range(f1.shape[0]):
        wave = f1[irow]
        tot = numpy.sum(wave)
        wave = scipy.signal.fftconvolve(wave, smear)[:nticks]
        f1[irow] = wave * tot / numpy.sum(wave)
    return f1
    

def plot_channels(f1, f2, chans, t0, t1, tit=""):
    x = numpy.arange(t0, t1)

    if tit:
        plt.title(tit)
    for c in chans:
        w1 = f1[c][t0:t1]
        w2 = f2[c][t0:t1]
        tot1 = numpy.sum(w1)
        tot2 = numpy.sum(w2)
        plt.plot(x, w1, label=f'ch{c} q={tot1:.0f} splat')
        plt.plot(x, w2, label=f'ch{c} q={tot2:.0f} signal')
    plt.legend()
    plt.xlabel("ticks")


@cli.command("plots")
@click.option("--channel-ranges", default="0,800,1600,2560",
              help="comma-separated list of channel idents defining ranges")
@click.option("--smear", default=0.0, help="Apply post-hoc smearing across tick dimension")
@click.option("--scale", default=0.0, help="Apply post-hoc multiplicative scale to waveforms")
@click.option("-o",'--output', default='plots.pdf')
@click.argument("depos")
@click.argument("drift")
@click.argument("splat")
@click.argument("signal")
def plots(channel_ranges, smear, scale, depos, drift, splat, signal, output, **kwds):
    '''
    Make plots comparing {depos,drift,signal,splat}.npz 
    '''
    # dumb sanity check user gave us files in the right order
    assert "depos" in depos
    assert "drift" in drift
    assert "splat" in splat
    assert "signal" in signal
    

    if channel_ranges:
        channel_ranges = list(map(int,channel_ranges.split(",")))
        channel_ranges = list(zip(channel_ranges[:-1], channel_ranges[1:]))

    fig = plt.figure()

    # d1 = load_depos(depos)
    d2 = load_depos(drift)
    f1,e1,o1 = load_frame(splat)
    print(f'extent={e1}')
    if smear:
        f1 = smear_tick(f1, smear)
    if scale:
        f1 *= scale

    f2,e2,o2 = load_frame(signal)

    with pages(output) as out:

        # plt.suptitle(f'{smear=:.1f} {scale=:.1f}')
        pgs = GridSpec(1,2, figure=fig, width_ratios = [7,0.2])
        gs = GridSpecFromSubplotSpec(2, 1, pgs[0,0])
        im1 = plot_frame(gs[0], f1, e1, o1, channel_ranges, which="splat")
        im2 = plot_frame(gs[1], f2, e2, o2, channel_ranges, which="signal")
        fig.colorbar(im2, cax=plt.subplot(pgs[0,1]))
        plt.tight_layout()
        out.savefig()

        plt.clf()

        fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
        plot_ticks(ax, e1[0], e1[1], d2, f1, f2, channel_ranges)
        out.savefig()

        plt.clf()
        plot_zoom(f1, f2, [0,400], [1400,1600],
                  f"splat - signal difference, begin of track, V-plane")
        out.savefig()

        plt.clf()
        plot_zoom(f1, f2, [4400, 4800], [1100, 1300],
                  f"splat - signal difference, end of track, V-plane")
        out.savefig()

        plt.clf()
        plot_zoom(f1, f2, [100,200], [1525,1560],
                  tit=f"splat - signal difference, begin of track, V-plane")
        out.savefig()

        plt.clf()
        plot_zoom(f1, f2, [4600, 4800], [1190, 1230],
                  tit=f"splat - signal difference, end of track, V-plane")
        out.savefig()

        plt.clf()
        plot_channels(f1, f2, [1540, 1530, 1520, 1510], 100, 350,
                      tit=f'V-plane start')
        out.savefig()

        plt.clf()
        plot_channels(f1, f2, [1200, 1210, 1220, 1230], 4550, 4750,
                      tit=f'V-plane end')
        out.savefig()

def relbias(a,b):
    rb = numpy.zeros_like(a)
    ok = b>0
    rb[ok] = a[ok]/b[ok] - 1
    return rb

@cli.command("plot-ssss")
@click.option("--channel-ranges", default="0,800,1600,2560",
              help="comma-separated list of channel idents defining ranges")
@click.option("--nsigma", default=3.0,
              help="Relative threshold on signal in units of number of sigma of noise width")
@click.option("-o",'--output', default='plots.pdf')
@click.argument("depos")
@click.argument("drift")
@click.argument("splat")
@click.argument("signal")
def plot_ssss(channel_ranges, depos, drift, splat, signal,
              nsigma, output, **kwds):
    '''
    Plot splat and sim+SP signals.
    '''
    # dumb sanity check user gave us files in the right order
    assert "depos" in depos
    assert "drift" in drift
    assert "splat" in splat
    assert "signal" in signal

    if channel_ranges:
        channel_ranges = list(map(int,channel_ranges.split(",")))
        channel_ranges = [slice(*cr) for cr in zip(channel_ranges[:-1], channel_ranges[1:])]

    fig = plt.figure()

    # d1 = load_depos(depos)
    d2 = load_depos(drift)

    f1,e1,o1 = load_frame(splat)
    f2,e2,o2 = load_frame(signal)

    tick_us = (e2[1] - e2[0])/(f1.shape[1]+1)
    print(f'{tick_us=}')


    with pages(output) as out:

        pgs = GridSpec(1,2, figure=fig, width_ratios = [7,0.2])
        gs = GridSpecFromSubplotSpec(2, 1, pgs[0,0])
        im1 = plot_frame(gs[0], f1, e1, o1, channel_ranges, which="splat")
        im2 = plot_frame(gs[1], f2, e2, o2, channel_ranges, which="signal")
        fig.colorbar(im2, cax=plt.subplot(pgs[0,1]))
        plt.tight_layout()
        out.savefig()

        byplane = list()

        # Per channel range plots.
        for pln, ch in enumerate(channel_ranges):
            letter = "UVW"[pln]

            spl = select_activity(f1, ch, nsigma)
            sig = select_activity(f2, ch, nsigma)
            print(f'selection: {ch} {spl.selection.shape} {sig.selection.shape}')

            # Find the bbox that bounds the biggest splat object.
            biggest = spl.plats.sort_by("sums")[-1]
            bbox = spl.plats.bboxes[biggest]
            print(f'{bbox=}')

            byplane.append((spl, sig, bbox))

            spl_act = spl.thresholded[bbox]
            sig_act = sig.thresholded[bbox]

            print('shapes',spl_act.shape, sig_act.shape)


            # bias of first w.r.t. second
            bias1 = relbias(sig_act, spl_act)
            bias2 = relbias(spl_act, sig_act)

            plt.clf()
            fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
            plt.suptitle(f'{letter}-plane')
            args=dict(aspect='auto')
            im1 = axes[0,0].imshow(sig_act, **args)
            fig.colorbar(im1, ax=axes[0,0])
            im2 = axes[0,1].imshow(spl_act, **args)
            fig.colorbar(im2, ax=axes[0,1])

            args = dict(args, cmap='jet', vmin=-50, vmax=50)

            im3 = axes[1,0].imshow(100*bias1, **args)
            fig.colorbar(im3, ax=axes[1,0])

            im4 = axes[1,1].imshow(100*bias2, **args)
            fig.colorbar(im4, ax=axes[1,1])

            axes[0,0].set_title(f'signal {nsigma=}')
            axes[0,1].set_title(f'splat {nsigma=}')

            axes[1,0].set_title(f'splat/signal - 1 [%]')
            axes[1,1].set_title(f'signal/splat - 1 [%]')

            chan_tit = 'chans (rel)'
            tick_tit = 'ticks (rel)'
            axes[0,0].set_ylabel(chan_tit)
            axes[1,0].set_ylabel(chan_tit)
            axes[1,0].set_xlabel(tick_tit)
            axes[1,1].set_xlabel(tick_tit)

            fig.subplots_adjust(right=0.85)

            plt.tight_layout()
            out.savefig()

        plt.clf()
        plt.suptitle('By channel')
        fig, axes = plt.subplots(nrows=2, ncols=3, sharey="row")
        for pln, (spl, sig, bbox) in enumerate(byplane):

            spl_qch = numpy.sum(spl.activity[bbox], axis=1)
            spl_tot = numpy.sum(spl_qch)
            sig_qch = numpy.sum(sig.activity[bbox], axis=1)
            sig_tot = numpy.sum(sig_qch)                

            print('ch tots:', spl_qch.shape, sig_qch.shape)

            ### MD SP paper definitions:
            # Charge bias - the mean fractional difference between the total
            # deconvolved charge and the true charge. In terms of a track (line
            # charge), the mean value of the distribution of each wire’s
            # integrated deconvolved charges from a range of wires will be used
            # to calculate the charge bias with respect to the true charge
            # within one wire pitch.
            #
            # Charge resolution - the standard deviation of the total
            # deconvolved charge relative to the mean deconvolved charge. In
            # analogy to the definition of charge bias, in terms of a track
            # (line charge), the distribution of each wire’s integrated
            # deconvolved charges from a range of wires will apply.
            #
            # Charge inefficiency - specifically associated with tracks (line
            # charges), the fraction of the wires which have ZERO deconvolved
            # charge (no ROI found). Note that these wires will not be involved
            # in the calculation of charge bias or charge resolution.
            ####
            
            # However, while a channel that has no signal is "inefficient" it is
            # possible for a channel to be "over efficient" by having signal
            # where there is no splat.  This has as much to do with (smeared)
            # "splat" being an imperfect "true signal" as it has anything to do
            # with "real" sim+SP signals.  So, we interpret MB's "ROI" above to
            # mean channels that have both signal and splat and then inefficient
            # channels have only splat and over efficient channels have only
            # signal.

            # either-or, exclude channels where both are zero
            eor   = numpy.logical_or (spl_qch  > 0, sig_qch  > 0)
            # both are nonzero
            both  = numpy.logical_and(spl_qch  > 0, sig_qch  > 0)
            nosig = numpy.logical_and(spl_qch  > 0, sig_qch == 0)
            nospl = numpy.logical_and(spl_qch == 0, sig_qch  > 0)

            neor = numpy.sum(eor)
            # inefficiency percent
            nnosig = numpy.sum(nosig)
            ineff = 100*nnosig/neor
            # over efficiency percent
            nnospl = numpy.sum(nospl)
            oveff = 100*nnospl/neor
            print(f'{neor=} {nnosig=} {nnospl=}')

            nbins = 50
            reldiff = 100.0*(spl_qch[both] - sig_qch[both])/spl_qch[both],
            bln = baseline_noise(reldiff, nbins, nbins//2)
            counts, edges = bln.hist

            # we plot:
            # 1. The total charge on each channel for signal and splat
            # 2. The histogram of their difference

            letter = "UVW"[pln]

            ax1,ax2 = axes[:,pln]

            ax1.plot(sig_qch, label='signal')
            ax1.plot(spl_qch, label='splat')
            ax1.set_xlabel('chans (rel)')
            ax1.set_ylabel('electrons')
            ax1.set_title(f'{letter} ineff={ineff:.1f}%')
            ax1.legend()

            ax2.step(edges[:-1], counts, label='data')
            model = gauss_func(edges[:-1], bln.A, bln.mu, bln.sigma)
            ax2.plot(edges[:-1], model, label='fit')
            ax2.set_title(f'mu={bln.mu:.2f}%\nsig={bln.sigma:.2f}%')
            ax2.set_xlabel('difference [%]')
            ax2.set_ylabel('counts')
            ax2.legend()

        plt.suptitle('(splat - signal) / splat')
        plt.tight_layout()
        out.savefig()

def main():
    cli()

if '__main__' == __name__:
    main()
