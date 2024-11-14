#!/usr/bin/env python
'''
This is NOT a waf tool but generic functions to be called from a waf
tool, in particular by wcb.py.

There's probably a wafier way to do this.

The interpretation of options are very specific so don't change them
unless you really know all the use cases.  The rules are:

If package is optional:

    - omitting all --with-NAME* options will omit use the package

    - explicitly using --with-NAME=false (or "no" or "off") will omit
      use of the package.

If package is mandatory:

    - omitting all --with-NAME* options will use pkg-config to find
      the package.

    - explicitly using --with-NAME=false (or "no" or "off") will
      assert.

In either case:

    - explicitly using --with-NAME=true (or "yes" or "on") will use
      pkg-config to find the package.

    - using --with-NAME* with a path will attempt to locate the
      package without using pkg-config

Note, actually, pgk-config fails often to do its job.  Best to always
use explicit --with-NAME=[<bool>|<dir>].
'''

import sys
import os.path as osp
from waflib.Logs import debug

def boolish(b):
    '''
    Convert b to Bool if true'ish or false'ish string, else return b.
    '''
    if not isinstance(b, str):
        return b
    if b.lower() in 'yes on true'.split():
        return True
    if b.lower() in 'no off false'.split():
        return False
    return b


def with_p(opt, name):
    '''
    Return whether opt is with named package.
    '''
    return boolish(getattr(opt, f'with_{name}', None) \
                   or getattr(opt, f'with_{name}_lib', None) \
                   or getattr(opt, f'with_{name}_libs', None) \
                   or getattr(opt, f'with_{name}_include', None))


def do_options(opt, *pkgdesc, **pkgmap):
    for name, desc in list(pkgdesc) + list(pkgmap.items()):
        _options(opt, name,
                 desc.get("incs", None),
                 desc.get("libs", None))
        

def _options(opt, name, incs=None, libs=None):
    lower = name.lower()
    opt = opt.add_option_group('%s Options' % name)
    opt.add_option('--with-%s'%lower, type=str, default=None,
                   help="give %s installation location" % name)
    if incs is None or len(incs):
        opt.add_option('--with-%s-include'%lower, type=str, 
                       help="give %s include installation location"%name)
    if libs is None or len(libs):
        opt.add_option('--with-%s-lib'%lower, type=str, 
                       help="list %s lib installation locations"%name)
        opt.add_option('--with-%s-libs'%lower, type=str, 
                       help="list %s link libraries"%name)
    return

def do_configure(cfg, *pkgdesc, **pkgmap):
    for name, desc in list(pkgdesc) + list(pkgmap.items()):
        _configure(cfg, name, **desc)
    

def _configure(ctx, name, incs=(), libs=(), bins=(), pcname=None, mandatory=True, extuses=()):
    lower = name.lower()
    UPPER = name.upper()
    extuses = list(extuses)

    instdir = getattr(ctx.options, 'with_'+lower, None)
    incdir = getattr(ctx.options, 'with_%s_include'%lower, None)
    libdir = getattr(ctx.options, 'with_%s_lib'%lower, None)
    bindir = getattr(ctx.options, 'with_%s_bin'%lower, None)
    if libs:
        maybe = getattr(ctx.options, 'with_%s_libs'%lower, None)
        if maybe:
            libs = [s.strip() for s in maybe.split(",") if s.strip()]
    
    debug (f"configure: {name=} {instdir=} {incdir=} {libdir=} {mandatory=} {extuses=}")

    if mandatory:
        if instdir:
            if instdir.lower() in ['no','off','false']:
                raise ValueError(f'{name} is mandatory, "{instdir}" not allowed')

    else:                       # optional
        if instdir and instdir.lower() in ['no','off','false']:
            debug ("%s: skipping non mandatory %s, use --with-%s=[yes|<dir>] to force" % (lower, name, lower))
            return

    # rely on package config
    no_dirs = not any([instdir,incdir,libdir,bindir])
    have_instdir = (instdir and instdir.lower() in ['yes','on','true'])
    have_pcname = pcname is not None
    if have_pcname and (no_dirs or have_instdir):
        ctx.start_msg('Checking for %s in PKG_CONFIG_PATH' % name)
        args = "--cflags"
        if libs:                # things like eigen may not have libs
            args += " --libs"
        ctx.check_cfg(package=pcname,  uselib_store=UPPER,
                      args=args, mandatory=mandatory)
        if 'HAVE_'+UPPER in ctx.env:
            ctx.end_msg("located by pkg-config")
        else:
            ctx.end_msg("missing from pkg-config")
            return
    else:                       # do manual setting

        if incs:
            if not incdir and instdir:
                incdir = osp.join(instdir, 'include')
            if incdir:
                setattr(ctx.env, 'INCLUDES_'+UPPER, incdir.split(","))

        if libs:
            if not libdir and instdir:
                libdir = osp.join(instdir, 'lib')
            if libdir:
                setattr(ctx.env, 'LIBPATH_'+UPPER, libdir.split(","))

        if bins:
            if not bindir and instdir:
                bindir = osp.join(instdir, 'bin')
            if bindir:
                setattr(ctx.env, 'PATH_'+UPPER, [bindir])
            
    
    # now check, this does some extra work in the caseof pkg-config

    if libs:
        ctx.start_msg("Location for %s libs" % (name,))
        for tryl in libs:
            ctx.check_cxx(lib=tryl, define_name='HAVE_'+UPPER+'_LIB',
                          use=[UPPER] + extuses, uselib_store=UPPER, mandatory=mandatory)
        ctx.end_msg(str(getattr(ctx.env, 'LIBPATH_' + UPPER, None)))

        ctx.start_msg("Libs for %s" % UPPER)
        have_libs = getattr(ctx.env, 'LIB_' + UPPER, None)
        ctx.end_msg(str(have_libs))
        if ctx.is_defined('HAVE_'+UPPER+'_LIB'):
            ctx.env['HAVE_'+UPPER] = 1
            debug('%s: HAVE %s libs' % (lower, UPPER))
            if pcname:
                ctx.env.REQUIRES += [pcname]
        else:
            debug('%s: NO %s libs' % (lower, UPPER))
            

    if incs:
        ctx.start_msg("Location for %s headers" % name)
        for tryh in incs:
            ctx.check_cxx(header_name=tryh, define_name='HAVE_'+UPPER+'_INC',
                          use=[UPPER] + extuses, uselib_store=UPPER, mandatory=mandatory)
        have_incs = getattr(ctx.env, 'INCLUDES_' + UPPER, None)
        ctx.end_msg(str(have_incs))
        if ctx.is_defined('HAVE_'+UPPER+'_INC'):
            ctx.env['HAVE_'+UPPER] = 1
            debug('%s: HAVE %s includes' % (lower, UPPER))
            if pcname:
                ctx.env.REQUIRES += [pcname]
        else:
            debug('%s: NO %s includes' % (lower, UPPER))

    if bins:
        ctx.start_msg("Bins for %s" % name)
        found_bins = list()
        for tryb in bins:
            ctx.find_program(tryb, var=tryb.upper(), mandatory=mandatory)
            found_bins += ctx.env[tryb.upper()]
        ctx.end_msg(str(found_bins))

