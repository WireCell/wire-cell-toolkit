#!/usr/bin/env python

# Copyright 2015-2023 Brookhaven National Laboratory for the benefit
# of the Wire-Cell Team.
# 
# This file is part of the wire-cell-toolkit project and distributed
# according to the LICENSE file provided as also part of this project.

import os
import sys
import subprocess

# fixme: move into waft/
from waflib.Logs import debug, info

TOP = '.'
APPNAME = 'WireCell'

def determine_version():
    proc = subprocess.run(["git", "describe", "--tags"], capture_output=True)
    if proc.returncode:
        if os.path.exists("version.txt"):
            return open("version.txt", "r").readlines()[0].strip()
        raise FileNotFoundError("Wire-Cell Toolkit must either be built from a git clone or a version.txt must be provided in the source distribution")
    version = proc.stdout.decode().strip()
    proc = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True)
    branch = proc.stdout.decode().strip()
    if branch == "master" or branch[0].isdecimal():
        return version
    return f'{branch}-{version}'

VERSION = determine_version()


# Valid log level identifiers
log_levels = "trace debug info warn error critical off "
log_levels = (log_levels + log_levels.upper()).split()

# to avoid adding tooldir="waft" in all the load()'s
sys.path.insert(0, os.path.realpath("./waft"))

def options(opt):
    opt.load("wcb")

    # this used in cfg/wscript_build
    opt.add_option('--install-config', type=str, default="",
                   help="Install configuration files for given experiment")

    opt.add_option('--build-mode', type=str, default="",
                   help="Force the build mode (default, detect based on git branch)")

    # fixme: add to spdlog entry in wcb.py
    opt.add_option('--with-spdlog-static', type=str, default="yes",
                   help="Def is true, set to false if your spdlog is not compiled (not recomended)")
    opt.add_option('--with-spdlog-active-level',
                   default = "debug",
                   choices = log_levels,
                   help="The compiled minimum log level for SPDLOG_<LEVEL>() macros (def=info)")

    opt.add_option('--cxxstd', default='c++17',
                   help="Set the value for the compiler's --std= option, default 'c++17'")
    opt.add_option("--with-nvtx", action="store_true", default=False,
                   help="Use NVTX if available (default: False)")

def is_development():
    '''
    Return True if we are sitting on a "development" branch.

    A "development" branch is any that is not master or a numerical release
    branch.
    '''
    if VERSION == "master" or VERSION[0].isdecimal():
        return False
    return True


def configure(cfg):
    # Save to BuildConfig.h and env
    cfg.define("WIRECELL_VERSION", VERSION)
    cfg.env.VERSION = VERSION

    # see waft/wcb.py for how this is consumed.
    cfg.env.IS_DEVELOPMENT = is_development()
    if cfg.env.IS_DEVELOPMENT:
        info(f"configuring for DEVELOPMENT ({VERSION})")
        cfg.define("WIRECELL_DEVELOPMENT", VERSION)
    else:
        info(f"configuring for RELEASE ({VERSION})")
        cfg.define("WIRECELL_RELEASE", VERSION)

    # See https://github.com/WireCell/wire-cell-toolkit/issues/337
    if not cfg.options.libdir and cfg.env.LIBDIR.endswith("lib64"):
        cfg.env.LIBDIR = cfg.env.LIBDIR[:-2]
        debug(f'configure: forcing: {cfg.env.LIBDIR=} instead of lib64/, use explicit --libdir if you really want it')
    
    # Set to DEBUG to activate SPDLOG_DEBUG() macros or TRACE to activate both
    # those and SPDLOG_TRACE() levels.
    lu = cfg.options.with_spdlog_active_level.upper()
    cfg.define("SPDLOG_ACTIVE_LEVEL", 'SPDLOG_LEVEL_' + lu, quote=False)

    # See comments at top of Exceptions.h for context.
    cfg.load('compiler_cxx')
    cfg.check_cxx(lib='backtrace', use='backtrace',
                  uselib_store='BACKTRACE',
                 define_name = 'HAVE_BACKTRACE_LIB',
                 mandatory=False, fragment="""
#include <backtrace.h>
int main(int argc,const char *argv[])
{
    struct backtrace_state *state = backtrace_create_state(nullptr,false,nullptr,nullptr);
}
                 """)
    if cfg.is_defined('HAVE_BACKTRACE_LIB'):
        cfg.env.LDFLAGS += ['-lbacktrace']

    # fixme: this should go away when everyone is up to at least boost
    # 1.78.
    cfg.check_cxx(header_name="boost/core/span.hpp", use='boost',
                  define_name = 'HAVE_BOOST_CORE_SPAN_HPP',
                  mandatory=False)


    # cfg.env.CXXFLAGS += ['-Wpedantic', '-Werror']
    cfg.env.CXXFLAGS += ['-std='+cfg.options.cxxstd.lower()]
    
    cfg.env.CXXFLAGS += ['-std=c++17']
    cfg.env.CXXFLAGS += ['-DEIGEN_HAS_CXX11']


    if cfg.options.with_spdlog_static.lower() in ("yes","on","true"):
        cfg.env.CXXFLAGS += ['-DSPDLOG_COMPILED_LIB=1']
    
    if cfg.options.with_nvtx:
        cfg.env.CXXFLAGS += ['-DHAVE_NVTX=1']

    # in principle, this should be the only line here.  Any cruft
    # above that has accrued should be seen as a fixme: move to
    # wcb/waf-tools.

    cfg.load("wcb")

    cfg.env.CXXFLAGS += ['-I.']

    info("Configured version %s" % VERSION)


def build(bld):
    ### we used to be set sloppiness globally.  Now we use #pragma to
    ### selectively quell warnings.  See util/docs/pragma.org for some info.
    #bld.env.CXXFLAGS += '-Wall -Wpedantic'.split()
    bld.env.CXXFLAGS += '-Wall'.split()

    bld.load('wcb')

    
def dumpenv(bld):
    bld.load('wcb')


def packrepo(bld):
    bld.load('wcb')
