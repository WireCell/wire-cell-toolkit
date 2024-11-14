# Aggregate all the waftools to make the main wscript a bit shorter.
# Note, this is specific to WC building

import generic
import os.path as osp
from waflib.Utils import to_list, subprocess
from waflib.Logs import debug, info, error, warn
from waflib.Build import BuildContext
from waflib.Configure import conf

# Use our own "wcsonnet" command instead of plain "jsonnet".
from smplpkgs import ValidationContext
ValidationContext.script_templates[".jsonnet"] = "${WCSONNET} ${SCRIPT}"

mydir = osp.dirname(__file__)

# These are packages descriptions which fit the generic functions.
# They will be checked in order so put any dependencies first.
package_descriptions = [

    # spdlog is "header only" but use library version for faster recompilation
    # wire-cell-util and ZIO both use this
    # Need to build with -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    ('spdlog',   dict(incs=['spdlog/spdlog.h'], libs=['spdlog'], pcname='spdlog')),

    ('ZLib',     dict(incs=['zlib.h'], libs=['z'], pcname='zlib')),
    ('FFTW',     dict(incs=['fftw3.h'], libs=['fftw3f'], pcname='fftw3f')),
    ('FFTWThreads', dict(libs=['fftw3f_threads'], pcname='fftw3f', mandatory=False)),
    ('JsonCpp',  dict(incs=["json/json.h"], libs=['jsoncpp'], pcname='jsoncpp')),

    ('Eigen',    dict(incs=["Eigen/Dense"], pcname='eigen3')),

    ('GLPK',    dict(incs=["glpk.h"], libs=["glpk"], mandatory=False)),

    # for faster parsing, consider:
    # ./wcb configure --with-jsonnet-libs=gojsonnet 
    # see https://github.com/WireCell/wire-cell-toolkit/issues/342
    ('Jsonnet',  dict(incs=["libjsonnet.h"], libs=['jsonnet'])),
    ('TBB',      dict(incs=["tbb/parallel_for.h"], libs=['tbb'], pcname='tbb', mandatory=False)),
    ('HDF5',     dict(incs=["hdf5.h"], libs=['hdf5'], pcname='hdf5', mandatory=False)),

    ('ZMQ',      dict(incs=["zmq.h"], libs=['zmq'], pcname='libzmq', mandatory=False)),
    ('CZMQ',     dict(incs=["czmq.h"], libs=['czmq'], pcname='libczmq', mandatory=False)),
    ('ZYRE',     dict(incs=["zyre.h"], libs=['zyre'], pcname='libzyre', mandatory=False)),
    ('ZIO',      dict(incs=["zio/node.hpp"], libs=['zio'], pcname='libzio', mandatory=False,
                      extuses=("ZYRE","CZMQ","ZMQ"))),

    # Note, this list may be modified (appended) in wscript files.
    # The list here represents the minimum wire-cell-toolkit requires.
]


def options(opt):

    # from here
    opt.load('boost')
    opt.load('smplpkgs')
    opt.load('rootsys')
    opt.load('libtorch')
    opt.load('cuda')
    opt.load('kokkos')
    #opt.load('protobuf')

    generic.do_options(opt, *package_descriptions)

    opt.add_option('--build-debug', default='-O2 -ggdb3',
                   help="Build with debug symbols")


def find_submodules(ctx):
    sms = list()
    for wb in ctx.path.ant_glob("**/wscript_build"):
        sms.append(wb.parent.name)
    sms.sort()
    return sms


def configure(cfg):
    info ('Compile options: %s' % cfg.options.build_debug)
    cfg.load('boost')
    cfg.load('smplpkgs')

    generic.do_configure(cfg, *package_descriptions)

    def with_p(name):
        return generic.with_p(cfg.options, name)

    debug(f'wcb: {cfg.options=}')
    if with_p('libtorch'):
        cfg.load('libtorch')

    if with_p('cuda'):
        cfg.load('cuda')

    if with_p('kokkos'):
        cfg.load('kokkos')

    if with_p('root'):
        cfg.load('rootsys')


    # Check for stuff not found in the wcb-generic way
    #
    # Record having a lib as a 'HAVE'.
    def haveit(one):
        one=one.upper()
        if 'LIB_'+one in cfg.env:
            cfg.env['HAVE_'+one] = 1
            debug('HAVE %s libs'%one)
        else:
            info('NO %s libs'%one)

    cfg.check_boost(lib='system filesystem graph thread program_options iostreams regex')
    haveit('boost')

    cfg.check(header_name="dlfcn.h", uselib_store='DYNAMO', lib=['dl'], mandatory=True)
    haveit('dynamo')

    cfg.check(features='cxx cxxprogram', lib=['pthread'], uselib_store='PTHREAD')
    haveit('pthread')


    cfg.env.CXXFLAGS += to_list(cfg.options.build_debug)
    cfg.env.CXXFLAGS += ['-DEIGEN_FFTW_DEFAULT=1']

    cfg.env.LIB += ['z']
    
    submodules = find_submodules(cfg)
    debug(f'wcb: {submodules=}')

    # Remove WCT packages if they an optional dependency wasn't found
    for pkg,ext in [
            ("root","HAVE_ROOTSYS"),
            ("tbb","HAVE_TBB LIB_FFTWTHREADS"),
            ("patrec", "HAVE_GLPK"),
            ("cuda","HAVE_CUDA"),
            ("hio", "INCLUDES_HDF5"),
            ("pytorch", "LIB_LIBTORCH"),
            ("zio", "LIB_ZIO LIB_ZYRE LIB_CZMQ LIB_ZMQ")
    ]:
        exts = to_list(ext)
        for have in exts:
            if have in cfg.env or have in cfg.env.define_key:
                continue
            if pkg in submodules:
                info('Removing package "%s" due to lack of external dependency "%s"'%(pkg,have))
                submodules.remove(pkg)

    cfg.env.SUBDIRS = submodules
    info ('Configured for submodules: %s' % (', '.join(submodules), ))
    bch = 'WireCellUtil/BuildConfig.h' 
    cfg.env['BuildConfig'] = bch
    info ('Writing build config: %s' % bch) 
    cfg.write_config_header(bch)

    # Define env vars for wire-cell-python CLI's.
    # Extend this list manually as more are developed!
    # ls -l wire-cell-python/wirecell/*/__main__.py
    for one in to_list("aux gen img ls4gan pgraph plot pytorch resp sigproc test util validate"):
        cmd = 'wirecell-' + one
        var = 'WC' + one.upper()
        cfg.find_program(cmd, var=var, mandatory=False)

    cfg.find_program("pandoc", var="PANDOC", mandatory=False)
    cfg.find_program("emacs", var="EMACS", mandatory=False)

    cfg.env.REQUIRES = list(set(cfg.env.REQUIRES))

    debug("dump: " + str(cfg.env))


# helper used in various tools
@conf
def cycle_group(bld, gname):
    if gname in bld.group_names:
        bld.set_group(gname)
    else:
        bld.add_group(gname)


def build(bld):
    bld.load('smplpkgs')

    bch = bld.path.find_or_declare(bld.env.BuildConfig)
    bld.install_files('${PREFIX}/include/' + bch.parent.name, bch)

    subdirs = bld.env.SUBDIRS
    debug ('wcb: subdirs %s' % (', '.join(subdirs), ))
    bld.recurse(subdirs)

    if hasattr(bld, "smplpkg_graph"):
        # fixme: this writes directly.  Should make it a task, including
        # running graphviz to produce PNG/PDF
        debug ("wcb: writing wct-depos.dot")
        bld.path.make_node("wct-deps.dot").write(str(bld.smplpkg_graph))

    # fixme: this replicates a chunk of code in smplpkgs
    org = bld.path.find_resource("README.org")
    if org:
        bld_dir = bld.root.find_dir(bld.out_dir)
        bld.cycle_group("documentation")
        for ext in ('html', 'pdf'):
            feat = 'org2'+ext
            if feat in bld.env.DOCS:
                out = org.parent.find_or_declare(org.name.replace(".org","."+ext))
                bld(name=org.parent.name + '-' + out.name.replace(".","-"),
                    features=feat, source=org, target=out)
                if ext == 'html':
                    bld.install_files(bld.env.DOCS_INSTALL_PATH, out,
                                      cwd=bld_dir, relative_trick=True)

    debug(f'wcb: {bld.smplpkg_names=}')
    link_libs = 'WireCellAux WireCellIface WireCellUtil'.split()
    debug(f'wcb: {link_libs=}')
    # Produce a pkg-config .pc file
    bld(source='wire-cell-toolkit.pc.in',
        name="pkg-config-file",
        LLIBS = ' '.join([f'-l{n}' for n in link_libs]),
        REQUIRES = ' '.join(bld.env.REQUIRES),
        install_path = '${LIBDIR}/pkgconfig/')


    # Produce a libtool .la file.  This needs one for each lib.
    # bld(source='libWireCellXxxx.la.in', 
    #     install_path = '${LIBDIR}')    

def dumpenv(bld):
    'dump build environment to stdout'

    # This is an ugly hack as the info needed is not in this build context.
    for app in "wire-cell wcsonnet wcwires".split():
        bld.env[app.upper().replace("-","_")] = bld.out_dir + "/apps/" + app

    for key in bld.env:
        val = bld.env[key]
        if isinstance(val, list):
            val = ' '.join(val)
        if "-" in key:
            warn("replace '-' with '_' in: %s" % key)
            key = key.replace("-","_")
        info('wcb: %s="%s"' % (key, val))


class DumpenvContext(BuildContext):
    cmd = 'dumpenv'
    fun = 'dumpenv'
        
def packrepo(bld):
    '''
    Pack existing test data repo to archive files.
    '''

    indir=bld.path.find_dir("build/tests/input")
    if indir.exists():
        cmd = "tar -C %s -cf input.tar input" % (indir.parent.abspath())
        info(cmd)
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            warn("archive failed")
    else:
        warn("no input files found")
    
    rels = bld.env.TEST_DATA_VERSIONS

    hdir = bld.path.find_dir("build/tests/history")
    hdirs = hdir.ant_glob("*")

    for vdir in hdir.children:
        if rels and vdir not in rels:
            info(f"skip {vdir}, not in releases")
            continue
        cmd = "tar -C %s -cf history-%s.tar history/%s" % (
            hdir.parent.abspath(), vdir, vdir)
        info(cmd)
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            warn("archive failed")
        
class PackrepoenvContext(BuildContext):
    cmd = 'packrepo'
    fun = 'packrepo'
