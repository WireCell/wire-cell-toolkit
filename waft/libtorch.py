import generic

def options(opt):
    generic._options(opt, "libtorch")

def configure(cfg):
    # warning, list of libs changes over version.

    libs = getattr(cfg.options, "with_libtorch_libs", None)
    with_cuda = generic.with_p(cfg.options, "cuda")
    if not libs:
        libs = ['torch', 'torch_cpu', 'c10']
        if with_cuda:
            libs += ['torch_cuda', 'c10_cuda']
        libs += ['c10']

    generic._configure(cfg, "libtorch", 
                       incs=["torch/script.h", "ATen/ATen.h"],
                       libs=libs,
                       mandatory=True)

    if with_cuda:               # fixme: why do we do this?
        setattr(cfg.env, 'LINKFLAGS_LIBTORCH',
                ['-Wl,--no-as-needed,-ltorch_cuda','-Wl,--as-needed'])
