import generic

from waflib import Task
from waflib.TaskGen import extension
from waflib.Tools import ccroot, c_preproc
from waflib.Configure import conf
from waflib.Logs import debug

import os

# from waf's playground
class cuda(Task.Task):
    run_str = '${NVCC} ${NVCCFLAGS} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F}${SRC} ${CXX_TGT_F} ${TGT}'
    color   = 'GREEN'
    ext_in  = ['.h']
    vars    = ['CCDEPS']
    scan    = c_preproc.scan
    shell   = False

@extension('.cu', '.cuda')
def c_hook(self, node):
    return self.create_compiled_task('cuda', node)

@extension('.cxx')
def cxx_hook(self, node):
    # override processing for one particular type of file
    if getattr(self, 'cuda', False):
        return self.create_compiled_task('cuda', node)
    else:
        return self.create_compiled_task('cxx', node)

def options(opt):
    generic._options(opt, "CUDA")

def configure(cfg):

    generic._configure(cfg, "CUDA", mandatory=False,
                       incs=["cuda.h"], libs=["cuda","cudart"], bins=["nvcc"])

    if not 'HAVE_CUDA' in cfg.env:
        return
    nvccflags = "-shared -Xcompiler -fPIC --gpu-architecture=all"
    nvccflags += os.environ.get("NVCCFLAGS","")
    cfg.env.NVCCFLAGS += nvccflags.strip().split()
    debug ("cuda: NVCCFLAGS = %s" % (' '.join(cfg.env.NVCCFLAGS)))
