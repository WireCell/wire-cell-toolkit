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

"""

from pathlib import Path
import click
import torch
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as pcolors

@click.group()
def cli():
    pass

def needs_log(tensor, title):
    if "_waveform" in title: return True
    if "field_response" in title: return True
    if "decon_kernel" in title and "amplitude" in title: return True
    if "_fr" in title and "phase" not in title: return True
    return False

def do_plot(tensor, title, pdf):
    if "complex" in str(tensor.dtype):
        do_plot(torch.abs(tensor), '(amplitude) ' + title, pdf)
        do_plot(torch.angle(tensor), '(phase) ' + title, pdf)
        return

    axes = plt.gca()
    logit = needs_log(tensor, title)

    ndim = tensor.dim()
    if ndim == 1:
        plt.plot(tensor)
        plt.title(title);
        if logit:
            axes.set_yscale('symlog')
    elif ndim == 2:
        norm = None
        if logit:
            pos = torch.abs(tensor)
            avmin = torch.min(pos).item()
            avmax = torch.max(pos).item()
            linthresh = avmax*1e-8
            vlim = avmax*1e-3
            norm = pcolors.SymLogNorm(linthresh=linthresh, 
                                      vmin=-vlim, vmax=vlim)
        plt.imshow(tensor, aspect='auto', norm=norm, interpolation='none', cmap='rainbow')
        plt.colorbar()
        plt.title(title)
    else:
        print(f'{ndim}-dimension tensor not plotted');
        return

    pdf.savefig(plt.gcf())
    plt.close();


@cli.command("plot-raw")
@click.argument("debug_filenames", nargs=-1)
def cmd_plot_raw(debug_filenames):
    '''
    Plot the tensors as found from one or more debug files.

    For each input file, a same-name PDF file with .pdf extension is made in
    current directory.
    '''
    for debug_filename in debug_filenames:

        fname = Path(debug_filename).stem

        output = fname + ".pdf"
        print(f'writing {output}')

        fp = torch.load(debug_filename, map_location='cpu', weights_only=False)
        with PdfPages(output) as pdf:
            for name in fp:
                do_plot(fp[name], f'{name} - {fname}', pdf)
    
def main():
    cli(obj=dict())

if '__main__' == __name__:
    main()
