#!/usr/bin/env python
'''
Functions for use in adc-noise-sig.smake
'''
import os
import numpy
import zipfile
from collections import defaultdict, namedtuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from numpy.ma import masked_array, masked_where, masked_invalid

# jsonnet -J cfg -e 'local a=import "test/pdsp-adc.jsonnet"; a.voltageperadc'
voltageperadc=3.4179687499999998e-10
# Why is toffset not exactly tbin?
toffset=873

def load_array(fname, index=0, with_name=False):
    if type(fname) != str:
        return fname

    fp = numpy.load(fname)
    name = list(fp.keys())[index]
    arr = fp[name]
    if with_name:
        return arr, name
    return arr

def transform_array(orig, copy, trans=numpy.round, dtype='f2'):
    arr, name = load_array(orig, with_name=True)
    arr = numpy.array(trans(arr), dtype=dtype)
    numpy.savez_compressed(copy, **{name:arr})

adcmm=10
voltmm=adcmm*voltageperadc

def coalesce(fnames, out):
    '''
    Bring arrays in many files to one file to help reduce snakemake overhead
    '''
    with zipfile.ZipFile(out, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:
        for fname in fnames:
            base = os.path.splitext(fname)[0]
            aname = "_".join(base.split("/")[-5:])
            print (aname)
            aname += ".npy"
            # not much faster
            # afile = f'/dev/shm/{aname}'
            afile = aname

            fp = numpy.load(fname);
            arr = fp[list(fp.keys())[0]]
            numpy.save(afile, arr)
            zf.write(afile, aname)
            os.remove(afile)
        

    # arrs = dict();
    # for fname in fnames:
    #     fp = numpy.load(fname);

    #     # .../adc-noise-sig/UVCGAN/fake/b/258/vlt.npz
    #     base = os.path.splitext(fname)[0]
    #     aname = "_".join(base.split("/")[-5:])
    #     arrs[aname] = fp[list(fp.keys())[0]]
    # numpy.savez_compressed(out, **arrs)

def coal_aname(wc, **kwds):
    '''
    Form coalesced array name based on object and any overriding keywords.
    '''
    d = wc._asdict()
    d.update(kwds)
    return "{generator}_{source}_{domain}_{event}_{tier}".format(**d)

def coal_dict(aname):
    parts = aname.split("_")
    return dict(generator=parts[0],
                source=parts[1],
                domain=parts[2],
                event=parts[3],
                tier=parts[4]);
Coal = namedtuple("Coal", "generator source domain event tier")
def coal_tuple(aname):
    if "/" in aname:            # path
        parts = os.path.splitext(aname)[0][-5:]
    else:                       # array name
        parts = aname.split("_")

    return Coal(*parts)

    

def plot_params(generator, source, tier, domain):
    'Return dictionary of plotting parameters based wildcards'

    ## gleen pedestal via:
    # wcsonnet cfg/test/pdsp-adc.jsonnet | jq '.adcpeds[0]'

    # Someone decided to be creative with array names...
    array_name = dict(
        UVCGAN="arr_0",
        ACLGAN="data",
        UGATIT="data",
        CycleGAN="x",
    )[generator];
    if source == "real":
        array_name = "arr_0"

    pp = dict(
        source=source, tier=tier, domain=domain,
        pedestal=0,
        #cmap='seismic',
        cmap='tab20',
        tickslice=slice(0, 256),
        aname = dict(raw=array_name, rnd=array_name, vlt='frame_raw_0',
                     src='frame_*_0', ref='frame_*_0', noi='frame_*_0', dig='frame_*_0',
                     sig='frame_gauss0_0')[tier],
        tier_title = dict(raw='Raw',
                          rnd='Rounded',
                          vlt='Scaled',
                          src='Loaded',
                          ref='Reframed',
                          noi='Noise addded',
                          dig='Digitized, no pedestal',
                          sig='Signal')[tier],
        tier_units = dict(raw='ADC',
                          rnd='ADC',
                          vlt='voltage',
                          src='voltage',
                          ref='voltage',
                          noi='voltage',
                          dig='ADC',
                          sig='electrons')[tier],
        source_title = dict(real="Sim",
                            fake="GAN",
                            diff="%-diff",
                            xddp="%-diff",
                            pimp="Paired metric"
                            )[source],
        domain_title = dict(a="q1d", b="2d",
                            ab="q1d / 2d", ba="2d / q1d")[domain],
    )

    if tier in ('raw', 'rnd', 'dig'):
        pp['vmin'] = -adcmm
        pp['vmax'] = adcmm
    if tier in ('sig', ):
        pp['vmin'] = 0
        pp['vmax'] = 2000
        pp['masked_value'] = 0
        pp['cmap'] = 'jet'
    if tier in ('vlt','ref', 'noi', 'src'):
        pp['vmin'] = -voltmm
        pp['vmax'] = voltmm

    if tier in ('ref','noi', 'dig', 'sig'):
        pp['tickslice'] = slice(toffset, toffset+256)
    if tier in ('dig'):
        pp['pedestal']=2350

    ts = pp['tickslice']
    pp['tickrange'] = f'{ts.start}:{ts.stop}'
    pp['zoomslice'] = (slice(0,256), ts)
    pp['zoomrange'] = f'0:256,{ts.start}:{ts.stop}'

    if source == "diff":
        pp['pedestal']=0
        pp['vmin'] = -10.0      # percent
        pp['vmax'] = 10.0
        pp['cmap']='jet' #'tab20'
    if source == "xddp":
        pp['vmin'] = -50.0      # percent
        pp['vmax'] = 50.0
        pp['cmap']='jet'
    if source == "pimp":        # paired image metric plot
        pp['vmin'] = -50.0      # percent
        pp['vmax'] = 50.0
        pp['cmap']='jet'
    return pp
    

def render_params(w):
    pp = plot_params(w.generator, w.source, w.tier, w.domain)
    args = "-c {cmap} -a '{aname}' --vmin {vmin} --vmax {vmax} -z '{zoomrange}' -b {pedestal}"
    return args.format(**pp)

masked_value = numpy.NaN

def make_diffarr(wc, arrs, pars):
    aarr, barr = arrs
    apar, bpar = pars
    rmf = (barr - aarr)
    den = 0.5*(numpy.abs(barr-bpar["pedestal"])+numpy.abs(aarr-apar["pedestal"]))
    nz = den != 0
    da = numpy.zeros_like(rmf)
    da[nz] = (100.0*rmf[nz]) / den[nz]

    if wc.tier in ("dig",):
        da[den==0] = 0          # turn divide by zero into zero 
    else:                       # o.w. supress very small denominator
        nz = den > 0.0001*numpy.max(den)
        da[~nz] = masked_value
    return da

def plot_diff(wc, files, params):
    apar,bpar,dpar = params
    aarr = load_array(files[0])
    barr = load_array(files[1])
    darr = make_diffarr(wc, [aarr,barr], [apar,bpar])
    
    fig, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)

    common = dict(aspect='auto')
    stit = '{tier_title} ({domain_title}) [{tier_units}] (sample {event})'.format(event=wc.event, **dpar)
    fig.suptitle(stit)
    ims = list()
    for ind, (pp,aa) in enumerate(zip(params,[aarr,barr,darr])):
        aa = aa[pp['zoomslice']] - pp['pedestal']
        if 'masked_value' in pp:
            aa[aa==pp['masked_value']] = numpy.NaN

        ax = axes[ind]
        ax.get_xaxis().set_visible(False);
        ax.get_yaxis().set_visible(False);
        tit = pp['source_title']
        ax.set_title(tit)
        im = ax.imshow(masked_invalid(aa),
                       vmin=pp['vmin'], vmax=pp['vmax'],
                       cmap=pp['cmap'], aspect='auto')
        ims.append(im)
    # fig.subplots_adjust(left=0.2, right=0.8)
    # cax_val = fig.add_axes([0.05, 0.15, 0.05, 0.7])
    # cax_dif = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # fig.colorbar(ims[0], cax=cax_val)
    # fig.colorbar(ims[2], cax=cax_dif)
    fig.colorbar(ims[0], ax=axes[0:2], location='left')
    fig.colorbar(ims[2], ax=axes[2:])

    # plt.savefig(files[2], bbox_inches='tight')
    plt.savefig(files[2])

def plot_frdp(wc, images, plot):
    arrs = numpy.load(images)
    fake = arrs[coal_aname(wc, source="fake")]
    real = arrs[coal_aname(wc, source="real")]

    fp = plot_params(wc.generator, "fake", wc.tier, wc.domain)
    rp = plot_params(wc.generator, "real", wc.tier, wc.domain)
    dp = plot_params(wc.generator, "diff", wc.tier, wc.domain)
    plot_diff(wc, [fake,real,plot], [fp, rp, dp])

def plot_xddp(wc, afile, bfile, pfile):
    arrs = numpy.load(images)
    a = arrs[coal_aname(wc, domain="a")]
    b = arrs[coal_aname(wc, domain="b")]

    ap = plot_params(wc.generator, wc.source, wc.tier, "a")
    bp = plot_params(wc.generator, wc.source, wc.tier, "b")
    dp = plot_params(wc.generator, "xddp", wc.tier, "ab")
    plot_diff(wc, [a, b, pfile], [ap, bp, dp])


def pimp(wc, coal, pfile):
    fp = plot_params(wc.generator, 'fake', wc.tier, wc.domain)
    rp = plot_params(wc.generator, 'real', wc.tier, wc.domain)
    # dp = plot_params(wc.generator, "pimp", wc.tier, wc.domain)
    evs = list()                # event numbers
    fss = list()
    rss = list()
    l1s = list()

    l2as = list()

    # arrs["real"]["<eve>"]
    arrs = defaultdict(dict)

    fp = numpy.load(coal)
    for aname in fp.keys():
        cobj = coal_tuple(aname)
        if cobj.tier != wc.tier:
            continue
        if cobj.domain != wc.domain:
            continue
        if cobj.generator != wc.generator:
            continue
        if cobj.source not in ("real", "fake"):
            continue

        arrs[cobj.source][cobj.event] = (cobj, fp[aname])
    
    reve = list(arrs["real"].keys())
    reve.sort()
    feve = list(arrs["fake"].keys())
    feve.sort()
    if reve != feve:
        raise ValueError("mismatch in event numbers")

    for eve in reve:
        fo,fa = arrs["fake"][eve]
        ro,ra = arrs["real"][eve]

        evs.append(eve)

        fss.append(numpy.sum(fa))
        rss.append(numpy.sum(ra))
        l1s.append(numpy.sum(numpy.abs(fa-ra)))
        l2 = numpy.sum((fa-ra)**2)
        l2as.append((l2, fa, ra))

    l2s = [t[0] for t in l2as]

    # get images for largest l2
    l2as.sort()
    ml2,fml2,rml2 = l2as[0]

    
    numpy.savez_compressed(pfile,
                           l1=numpy.array(l1s),
                           l2=numpy.array(l2s),
                           fml2=fml2,
                           rml2=rml2,
                           ev=reve,
                           fsum=numpy.array(fss),
                           rsum=numpy.array(rss))
    
    
def plot_pimp(wc, ifile, ofile):
    arrs = numpy.load(ifile)

    fig = plt.figure()
    gs = GridSpec(2, 6, figure=fig)

    axes = [
        fig.add_subplot(gs[0, 0:2]),
        fig.add_subplot(gs[0, 2:4]),
        fig.add_subplot(gs[0, 4:6])
    ]

    #fig, axes2 = plt.subplots(nrows=2, ncols=3)

    stit = f'{wc.tier} ({wc.generator}, domain "{wc.domain}")'
    fig.suptitle(stit)

    axes[0].set_title('Total')
    axes[0].hist(-arrs['fsum'], label='GAN')
    axes[0].hist(-arrs['rsum'], label='Sim')
    axes[0].legend()

    axes[1].set_title('L1(GAN,Sim)')
    axes[1].hist(arrs['l1'])

    axes[2].set_title('L2(GAN,Sim)')
    axes[2].hist(arrs['l2'])

    axes = [
        fig.add_subplot(gs[1, 0:3]),
        fig.add_subplot(gs[1, 3:6])
    ]

    fml2 = arrs['fml2']
    rml2 = arrs['rml2']

    if fml2.shape == (256,256): # patch
        im = axes[0].imshow(fml2)
        plt.colorbar(im);
        im = axes[1].imshow(rml2)
        plt.colorbar(im);
    else:                       # full
        im = axes[0].imshow(fml2[:256,1000:1256])
        plt.colorbar(im);
        im = axes[1].imshow(rml2[:256,1000:1256])
        plt.colorbar(im);
    plt.savefig(ofile)

def main_coal_file(coalfile, *evefiles):
    coalesce(evefiles, coalfile)


def main_pimp(*evefiles):
    parts = os.path.splitext(evefiles[0])[0].split("/")[-5:]
    parts[3] = "all";
    wc = Coal(*parts)
    print (wc)

    name = f'{wc.generator}-{wc.domain}-{wc.tier}'
    coalfile = name + '-coal.npz'
    pimpfile = name + '-pimp.npz'

    if os.path.exists(coalfile):
        print(f"file exists, not remaking: {coalfile}")
    else:
        coalesce(evefiles, coalfile)
        print(f'coalesced to {coalfile}')

    if os.path.exists(pimpfile):
        print(f"file exists, not remaking: {pimpfile}")
    else:
        pimp(wc, coalfile, pimpfile)
        print(f'pimped to {pimpfile}')

    plotfile = name + "-pimp.pdf"
    plot_pimp(wc, pimpfile, plotfile)
    print(f'plotted to {plotfile}')    

if '__main__' == __name__:
    import sys

    meth = globals()["main_" + sys.argv[1]]
    got = meth(*sys.argv[2:])
    print(got)

    # coalesce(sys.argv[2:], sys.argv[1])
    # /home/bv/wrk/wct/point-cloud/toolkit/build/sigproc/test/adc-noise-sig/ACLGAN/pimp/a/sig.npz
    
