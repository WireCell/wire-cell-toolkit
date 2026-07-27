#!/usr/bin/env -S uv run --script
# -*- python -*-
# /// script
# dependencies = ["click", "torch", "numpy", "matplotlib"]
# index-url.pytorch_cpu = "https://download.pytorch.org/whl/cpu"
# explicit.pytorch_cpu = true
# package.torch.source = "pytorch_cpu"
# ///
"""
Running a KernelConvolve job with a "debug_filename" option for it and its
kernels will produce a set of torch pickle files.  This script will plot them.

Example usage

$ wire-cell -l stderr -L debug  -A output=test-tdm-decon-%s.npz spng/test/test-tdm-decon.jsonnet
$ wcpy plot frame-image test-tdm-decon-adc.npz -o test-tdm-decon-adc.pdf --transform median
$ wcpy plot frame-image test-tdm-decon-sig.npz -o test-tdm-decon-sig.pdf
$ check-decon plot-raw *gauss.pkl

Other jobs that can dump pkl files:

- dnnroi-training.jsonnet -A dump=crossviews,wiener,dnnroi,stack,wthresh [etc]

"""

from pathlib import Path
import click
import torch
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as pcolors
import numpy

@click.group()
def cli():
    pass

def needs_log(tensor, title):
    if "_waveform" in title: return True
    if "field_response" in title: return True
    if "decon_kernel" in title and "amplitude" in title: return True
    if "_fr" in title and "phase" not in title: return True
    return False

def maybe_abs(tensor):
    if "complex" in str(tensor.dtype):
        return torch.abs(tensor)
    return tensor

def maybe_median(tensor, dim=-1):
    '''
    Subtract the median if the tensor holds numbers.
    '''
    if "bool" in str(tensor.dtype):
        return tensor
    return tensor - torch.median(tensor, dim=dim, keepdim=True).values
    

def count_plots(tensor):
    '''Return how many savefig calls do_plot will make for this tensor.'''
    ndim = tensor.dim()
    if ndim == 3:
        return sum(count_plots(t) for t in tensor)
    if "complex" in str(tensor.dtype):
        return 2
    if ndim in (1, 2):
        return 1
    return 0


class PngOutput:
    '''Drop-in for PdfPages that saves one PNG per savefig call.'''

    def __init__(self, stem, total):
        self.stem = stem
        self.ndigits = len(str(max(total - 1, 0)))
        self._count = 0

    def savefig(self, fig):
        filename = f"{self.stem}-{self._count:0{self.ndigits}d}.png"
        print(f'writing {filename}')
        fig.savefig(filename, dpi=150)
        self._count += 1

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass


def do_plot(tensor, title, pdf):
    ndim = tensor.dim()
    if ndim == 3:
        print(f'assuming 3D tensor of shape {tensor.shape} is batched on dim=0')
        for ind, one in enumerate(tensor):
            do_plot(one, f'{title} - (batch {ind})', pdf)
        return


    print(f'plotting {tensor.dtype} {tensor.shape} {title} sum={torch.sum(tensor)}')

    if "complex" in str(tensor.dtype):
        do_plot(torch.abs(tensor), '(amplitude) ' + title, pdf)
        do_plot(torch.angle(tensor), '(phase) ' + title, pdf)
        return

    axes = plt.gca()
    logit = needs_log(tensor, title)

    tensor = maybe_median(tensor)

    if ndim == 1:
        plt.plot(tensor)
        plt.title(title);
        if logit:
            axes.set_yscale('symlog')
    elif ndim == 2:
        norm = None
        if logit:
            pos = maybe_abs(tensor)
            avmin = torch.min(pos).item()
            avmax = torch.max(pos).item()
            linthresh = avmax*1e-8
            vlim = avmax*1e-3
            if "bool" in str(tensor.dtype):
                norm = None
            else:
                print(f'norm by {linthresh=} {vlim=}')
                norm = pcolors.SymLogNorm(linthresh=linthresh,
                                          vmin=-vlim, vmax=vlim)
        plt.imshow(tensor, aspect='auto', norm=norm, interpolation='none', cmap='rainbow')
        plt.colorbar()
        plt.title(title)
    else:
        print(f'{ndim}-dimension tensor not plotted');
        return
    plt.grid(True)
    pdf.savefig(plt.gcf())
    plt.close();


@cli.command("plot-raw")
@click.option("-f", "--fmt", type=click.Choice(["pdf", "png"]), default="pdf",
              show_default=True, help="Output format.")
@click.argument("debug_filenames", nargs=-1)
def cmd_plot_raw(fmt, debug_filenames):
    '''
    Plot the tensors as found from one or more debug files.

    For PDF format, a single .pdf file is written per input using PdfPages.
    For PNG format, one .png file per plot is written, named <stem>-NN.png
    where NN is a zero-padded integer with enough digits for all plots in
    that file.
    '''
    for debug_filename in debug_filenames:

        fname = Path(debug_filename).stem

        fp = torch.load(debug_filename, map_location='cpu', weights_only=False)

        if fmt == "pdf":
            output = fname + ".pdf"
            print(f'writing {output}')
            with PdfPages(output) as pdf:
                for name in fp:
                    do_plot(fp[name], f'{name} - {fname}', pdf)
        else:
            total = sum(count_plots(fp[name]) for name in fp)
            with PngOutput(fname, total) as png:
                for name in fp:
                    do_plot(fp[name], f'{name} - {fname}', png)
    
@cli.command("tonp")
@click.argument("debug_filenames", nargs=-1)
def cmd_tonp(debug_filenames):
    '''
    Convert .pkl files to WCT numpy .npz files.

    This makes a huge assumption that the torch tensor rows are in channel ID
    order and the channel IDs are sequential counts starting from zero.  In
    general this is NOT true.

    This should produce a .npz file that teepeesee can consume.
    '''

    for debug_filename in debug_filenames:

        fname = Path(debug_filename).stem

        output = fname + ".npz"
        print(f'writing {debug_filename} -> {output}')

        fp = torch.load(debug_filename, map_location='cpu', weights_only=False)
        tensors = list()
        for name, arr in fp.items():
            print(name, arr.shape, arr.dtype)
            if len(arr.shape) > 2:
                arr = torch.squeeze(arr)
            tensors.append(arr)
        tensors = torch.vstack(tuple(tensors))
        print(f'{tensors.shape=}')
        channels = torch.arange(tensors.shape[0]) # this makes huge assumptions!
        tickinfo = torch.tensor([0,500,0])
        arrays = {
            'frame_*_0': tensors.numpy(),
            'channels_*_0': channels.numpy(),
            'tickinfo_*_0': tickinfo.numpy()
        }
        numpy.savez_compressed(output, **arrays)
            


def main():
    cli(obj=dict())

if '__main__' == __name__:
    main()
