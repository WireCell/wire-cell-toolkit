import generic

from waflib import Task
from waflib.TaskGen import extension
from waflib.Tools import ccroot, c_preproc
from waflib.Configure import conf
from waflib.Logs import debug

import os

# from waf's playground
class kokkos_gcc(Task.Task):
    run_str = '${CXX} ${KOKKOS_CXXFLAGS} ${CXXFLAGS} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F}${SRC} ${CXX_TGT_F} ${TGT}'
    color   = 'GREEN'
    ext_in  = ['.h']
    vars    = ['CCDEPS']
    scan    = c_preproc.scan
    shell   = False

class kokkos_cuda(Task.Task):
    run_str = '${NVCC} ${KOKKOS_NVCCFLAGS} ${NVCCFLAGS} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F}${SRC} ${CXX_TGT_F} ${TGT}'
    color   = 'GREEN'
    ext_in  = ['.h']
    vars    = ['CCDEPS']
    scan    = c_preproc.scan
    shell   = False

@extension('.kokkos')
def kokkos_hook(self, node):
    options = getattr(self.env, 'KOKKOS_OPTIONS', None)
    if 'cuda' in options:
        debug('kokkos: use nvcc on ' + str(node))
        return self.create_compiled_task('kokkos_cuda', node)
    else:
        debug('kokkos: use gcc on ' + str(node))
        return self.create_compiled_task('kokkos_gcc', node)

def options(opt):
    generic._options(opt, "KOKKOS")
    opt.add_option('--kokkos-options', type=str, help="cuda, ...")

def configure(cfg):
    generic._configure(cfg, "KOKKOS", mandatory=False,
                       incs=["Kokkos_Macros.hpp"], libs=["kokkoscore", "kokkoscontainers", "dl"], bins=["nvcc"])

    options = getattr(cfg.options, 'kokkos_options', None)
    setattr(cfg.env, 'KOKKOS_OPTIONS', options)
    options = getattr(cfg.env, 'KOKKOS_OPTIONS', None)
    cfg.start_msg("KOKKOS_OPTIONS:")
    cfg.end_msg(str(options))

    if not 'HAVE_KOKKOS' in cfg.env:
        return
    nvccflags = "-x cu -shared -Xcompiler -fPIC "
    # nvccflags += "--std=c++11 "
    nvccflags += "-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -expt-extended-lambda -arch=sm_75 -Xcompiler -fopenmp "
    nvccflags += os.environ.get("NVCCFLAGS","")
    cfg.env.KOKKOS_NVCCFLAGS += nvccflags.strip().split()
    cxxflags = " -x c++ "
    cfg.env.KOKKOS_CXXFLAGS += cxxflags.strip().split()
    cfg.start_msg("KOKKOS_NVCCFLAGS:")
    cfg.end_msg(str(cfg.env.KOKKOS_NVCCFLAGS))
    cfg.start_msg("KOKKOS_CXXFLAGS:")
    cfg.end_msg(str(cfg.env.KOKKOS_CXXFLAGS))
