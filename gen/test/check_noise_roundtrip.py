#!/usr/bin/env python
import numpy
import functools
import click
cmddef = dict(context_settings = dict(help_option_names=['-h', '--help']))

@click.group(**cmddef)
@click.pass_context
def cli(ctx):
    '''
    wirecell-plot command line interface
    '''
    ctx.ensure_object(dict)


def calc_stats(x):
    n = x.size
    mu = numpy.mean(x)
    arel = numpy.abs(x-mu)
    rms = numpy.sqrt( (arel**2).mean() )
    outliers = [sum(arel >= sigma*rms) for sigma in range(0,11)]
    return [n,mu,rms]+outliers
        

@cli.command("stats")
@click.option("-a", "--array", default="frame_*_0", help="array name")
@click.option("-c", "--channels", default="800,800,960", help="comma list of channel counts per plane in u,v,w order")
@click.argument("npzfile")
def stats(array, channels, npzfile):
    '''
    Return (print) stats on the time distribution of a frame.

    '''
    fp = numpy.load(npzfile)
    array = fp[array]
    
    channels = [int(c) for c in channels.split(',')]
    chan0=0
    for chan, letter in zip(channels,"UVW"):
        plane = array[chan0:chan0+chan,:]
        plane = (plane.T - numpy.median(plane, axis=1)).T

        tsum = plane.sum(axis=0)/plane.shape[0]
        csum = plane.sum(axis=1)/plane.shape[1]
        # print(plane.shape, tsum.size, csum.size, (chan0, chan0+chan), numpy.sum(plane))

        print(' '.join([letter, 't'] + list(map(str,calc_stats(tsum)))))
        print(' '.join([letter, 'c'] + list(map(str,calc_stats(csum)))))

        chan0 += chan

# Decorator for a CLI that is common to a couple commands.
def frame_to_image(func):
    from matplotlib import colormaps

    @click.option("-a", "--array", default="frame_*_0", help="array name")
    @click.option("-c", "--channels", default="800,800,960", help="comma list of channel counts per plane in u,v,w order")
    @click.option("-f", "--format", default=None, help="Output file format, def=auto")
    @click.option("-o", "--output", default=None, help="Output file, def=stdout")
    @click.option("-C", "--cmap", default="gist_ncar", help="Color map name def=gist_ncar")
    @click.argument("npzfile")
    @functools.wraps(func)
    def wrapper(*args, **kwds):

        fmt = kwds["format"]
        out = kwds["output"]
        if fmt is None:
            if out is None or len(out.split(".")) == 1:
                kwds["format"] = png
            else:
                kwds["format"] = out.split(".")[-1]
        if out is None:
            kwds["output"] = "/dev/stdout"

        channels = list(map(int, kwds["channels"].split(",")))
        if len(channels) != 3:
            raise click.BadParameter("must give 3 channel group counts")
        chrange = list()
        for i,c in enumerate(channels):
            if not i:
                chrange.append((0, c))
            else:
                l = chrange[i-1][1]
                chrange.append((l, l+c))
        kwds["channels"] = chrange
        
        fname = kwds.pop("npzfile")
        kwds['fname'] = fname
        fp = numpy.load(fname)
        
        aname = kwds.pop("array")
        if aname in fp:
            kwds["array"] = fp[aname]
        else:
            have = '", "'.join(fp.keys())
            raise click.BadParameter(f'array "{aname}" not in "{fname}".  try: "{have}"')
        kwds['aname'] = aname

        kwds["cmap"] = colormaps[kwds["cmap"]]

        return func(*args, **kwds)
    return wrapper



@cli.command("plot")
@frame_to_image
def plot(array, channels, cmap, format, output, aname, fname):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize
    
    # layout = "constrained"
    layout = "tight"            # this gives warning but closer to what I want.
    fig = plt.figure(layout=layout, figsize=(10,8))
    fig.suptitle(f'array "{aname}" from {fname}\nand time/channel projected means')
    # base + 
    # [0:mean, 1:image, 2:colorbar]
    # [      , 4:mean,            ]
    # x3

    nplns=3
    nrows=2
    ncols=2
    gridspec = GridSpec(nplns*nrows+1, ncols,
                        figure=fig,
                        height_ratios=[1,3,1,3,1,3,1], width_ratios=[1,30],
                        left=0.05, right=0.95, hspace=0.0001, wspace=0.0001)
    def gs(pln, row, col):
        return gridspec[(pln*nrows+1)*ncols + row*ncols + col]

    steerx = None
    steerc = None
    steert = None

    normalizer=Normalize(numpy.min(array), numpy.max(array))
    cb = cm.ScalarMappable(norm=normalizer)
    
    aximgs = list()
    for pln, letter in enumerate("UVW"):

        pgs = lambda r,c: gs(pln, r, c)

        base = pln*6

        if steerx is None:
            aximg = fig.add_subplot(pgs(0, 1))
            steerx = aximg
        else:
            aximg = fig.add_subplot(pgs(0, 1), sharex=steerx)

        aximg.set_axis_off()
        aximgs.append(aximg)

        if steerc is None:
            axmu0 = fig.add_subplot(pgs(0, 0), sharey=aximg)
            steerc = axmu0
        else:
            axmu0 = fig.add_subplot(pgs(0, 0), sharey=aximg, sharex=steerc)

        if steert is None:
            axmu1 = fig.add_subplot(pgs(1, 1), sharex=steerx)
            steert = axmu1
        else:
            axmu1 = fig.add_subplot(pgs(1, 1), sharex=steerx, sharey=steert)
        axmus = [axmu1, axmu0]

        if pln == 0:
            plt.setp( axmu1.get_xticklabels(), visible=False)
            axmu0.xaxis.tick_top()
            axmu0.tick_params(axis="x", labelrotation=90)
        if pln == 1:
            plt.setp( axmu1.get_xticklabels(), visible=False)
            axmu0.set_ylabel('channels')
            plt.setp( axmu0.get_xticklabels(), visible=False)
        if pln == 2:
            axmu1.set_xlabel("ticks")
            plt.setp( axmu0.get_xticklabels(), visible=False)

        axmu0.ticklabel_format(useOffset=False)
        axmu1.ticklabel_format(useOffset=False)

        crange = channels[pln]
        achans = numpy.array(range(*crange))
        aticks = numpy.array(range(array.shape[1]))
        xses = [aticks, achans]
        plane = array[achans,:]
        im = aximg.imshow(plane, cmap=cmap, norm=normalizer,
                          extent=(aticks[0], aticks[-1], crange[1], crange[0]),
                          interpolation="none", aspect="auto")

        for axis in [0,1]:
            mu = plane.sum(axis=axis)/plane.shape[axis]
            axmu = axmus[axis]
            xs = xses[axis]

            if axis: 
                axmu.plot(mu, xs)
            else:
                axmu.plot(xs, mu)
                      
            
    axcb = fig.add_subplot(gridspec[1])
    fig.colorbar(cb, cax=axcb, ax=aximgs, cmap=cmap, location='top')

    fig.savefig(output, format=format, bbox_inches='tight')

@cli.command("image")
@frame_to_image
def image(array, channels, cmap, format, output, aname, fname):
    '''
    Dump array to image, ignoring channels.
    '''
    import matplotlib.image

    matplotlib.image.imsave(output, array, format=format, cmap=cmap)



@cli.command("configured-spectra")
@click.option("-t","--type", default="GroupNoiseModel", help="component type name")
@click.option("-n","--name", default="inco", help="component instance name")
@click.argument("cfgfile")
@click.argument("output")
def configured_spectra(type, name, cfgfile, output):
    from wirecell.util import jsio
    from wirecell.sigproc.noise.plots import plot_many

    got = None
    for one in jsio.load(cfgfile):
        if one['type'] == type and one['name'] == name:
            got=one['data']['spectra']
            break
    if got is None:
        raise click.BadParameter(f'failed to find node {type}:{name}')

    zero_suppress = True
    plot_many(got, output, zero_suppress)
    

def main():
    cli(obj=dict())

if '__main__' == __name__:
    main()
    
