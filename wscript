#!/usr/bin/env python

# Copyright 2015-2023 Brookhaven National Laboratory for the benefit
# of the Wire-Cell Team.
# 
# This file is part of the wire-cell-toolkit project and distributed
# according to the LICENSE file provided as also part of this project.

import os

# fixme: move into waft/
from waflib.Build import BuildContext
from waflib.Logs import debug, info, error, warn

TOP = '.'
APPNAME = 'WireCell'
VERSION = os.popen("git describe --tags").read().strip()


def options(opt):
    opt.load("wcb")

    # this used in cfg/wscript_build
    opt.add_option('--install-config', type=str, default="",
                   help="Install configuration files for given experiment")

    # fixme: add to spdlog entry in wcb.py
    opt.add_option('--with-spdlog-static', type=str, default="yes",
                   help="Def is true, set to false if your spdlog is not compiled (not recomended)")

def configure(cfg):
    # get this into config.h
    cfg.define("WIRECELL_VERSION", VERSION)

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


    # fixme: should go into wcb.py
    cfg.find_program("jsonnet", var='JSONNET')

    # boost 1.59 uses auto_ptr and GCC 5 deprecates it vociferously.
    cfg.env.CXXFLAGS += ['-Wno-deprecated-declarations']
    cfg.env.CXXFLAGS += ['-Wall', '-Wno-unused-local-typedefs', '-Wno-unused-function']
    # cfg.env.CXXFLAGS += ['-Wpedantic', '-Werror']

    if cfg.options.with_spdlog_static.lower() in ("yes","on","true"):
        cfg.env.CXXFLAGS += ['-DSPDLOG_COMPILED_LIB=1']

    # in principle, this should be the only line here.  Any cruft
    # above that has accrued should be seen as a fixme: move to
    # wcb/waf-tools.
    cfg.load("wcb")

    cfg.env.CXXFLAGS += ['-I.']

    print("Configured version", VERSION)
#    print(cfg.env)

def build(bld):
    bld.load('wcb')



## fixme: move into waft
def dumpenv(bld):
    'print build environment'
    for key in bld.env:
        val = bld.env[key]
        if isinstance(val, list):
            val = ' '.join(val)
        if "-" in key:
            warn("replace '-' with '_' in: %s" % key)
            key = key.replace("-","_")
        print('%s="%s"' % (key, val))


class DumpenvContext(BuildContext):
    cmd = 'dumpenv'
    fun = 'dumpenv'

