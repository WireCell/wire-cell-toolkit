import os
import os.path as osp
import waflib
from waflib.Utils import to_list
from waflib.Configure import conf
from waflib.Logs import debug

_tooldir = osp.dirname(osp.abspath(__file__))


def options(opt):
    opt = opt.add_option_group('ROOT Options')
    opt.add_option('--with-root', default=None,
                   help="give ROOT installation location")
    return

@conf
def check_root(cfg, mandatory=False):
    instdir = cfg.options.with_root

    if instdir and instdir.lower() in ['no','off','false']:
        if mandatory:
            raise RuntimeError("ROOT is mandatory but disabled via command line")
        debug ("root: optional ROOT dependency disabled by command line")
        return

    cfg.env.CXXFLAGS += ['-fPIC']

    path_list = list()
    for topdir in [getattr(cfg.options, 'with_root', None), os.getenv('ROOTSYS', None)]:
        if topdir:
            path_list.append(osp.join(topdir, "bin"))
        
    kwargs = dict(path_list=path_list)

    cfg.find_program('root-config', var='ROOT_CONFIG', mandatory=mandatory, **kwargs)
    if not 'ROOT_CONFIG' in cfg.env:
        if mandatory:
            raise RuntimeError("root-config not found but ROOT required")
        debug ("root: skipping non mandatory ROOT, use --with-root to force")
        return

    cfg.check_cfg(path=cfg.env['ROOT_CONFIG'], uselib_store='ROOTSYS',
                  args = '--cflags --libs --ldflags', package='', mandatory=mandatory)


    cfg.find_program('rootcling', var='ROOTCLING', path_list=path_list, mandatory=mandatory)
    # cfg.find_program('rootcint', var='ROOTCINT', path_list=path_list, mandatory=mandatory)
    # cfg.find_program('rlibmap', var='RLIBMAP', path_list=path_list, mandatory=False)

    cfg.check_cxx(header_name="Rtypes.h", use='ROOTSYS',
                  mandatory=mandatory)


    return

def configure(cfg):
    cfg.check_root()

@conf
def gen_rootcling_dict(bld, name, linkdef, headers = '', includes = '', use=''):
    '''
    rootcling -f NAMEDict.cxx -rml libNAME.so -rmf libNAME.rootmap myHeader1.h myHeader2.h ... LinkDef.h
    '''
    use = to_list(use) + ['ROOTSYS']
    includes = to_list(includes)
    # for u in use:
    #     more = bld.env['INCLUDES_'+u]
    #     debug('root: USE(%s)=%s: %s' % (name, u, more))
    #     includes += more

    # make into -I...
    incs = list()
    for inc in includes:
        debug(f'root: INC({name}): {inc}')
        if inc.startswith('/'):
            newinc = '-I%s' % inc
        else:
            newinc = '-I%s' % bld.path.find_dir(inc).abspath()
        if not newinc in incs:
            incs.append(newinc)
    incs = ' '.join(incs)
    debug('root: INCS(%s): %s' % (name, str(incs)))
    
    dict_src = name + 'Dict.cxx'
    # dict_lib = 'lib' + name + 'Dict.so' # what for Mac OS X?
    # dict_map = 'lib' + name + 'Dict.rootmap'
    # dict_pcm =         name + 'Dict_rdict.pcm'
    dict_lib = 'lib' + name + '.so' # what for Mac OS X?
    dict_map = 'lib' + name + '.rootmap'
    dict_pcm =         name + 'Dict_rdict.pcm'

    if type(linkdef) == type(""):
        linkdef = bld.path.find_resource(linkdef)
    source_nodes = to_list(headers) + [linkdef]
    sources = ' '.join([x.abspath() for x in source_nodes])
    debug(f'root: SOURCES({name}): {sources}')
    rule = '${ROOTCLING} -f ${TGT[0].abspath()} -rml %s -rmf ${TGT[1].abspath()} %s %s' % (dict_lib, incs, sources)
    bld(source = source_nodes,
        target = [dict_src, dict_map, dict_pcm],
        rule=rule, use = use)

    # bld.shlib(source = dict_src,
    #           target = name+'Dict',
    #           includes = includes,
    #           use = use + [name])

    bld.install_files('${PREFIX}/lib/', dict_map)
    bld.install_files('${PREFIX}/lib/', dict_pcm)


@conf
def gen_rootcint_dict(bld, name, linkdef, headers = '', includes=''):
    '''Generate a rootcint dictionary, compile it to a shared lib,
    produce its rootmap file and install it all.
    '''
    headers = to_list(headers)
    incs = ['-I%s' % bld.path.find_dir(x).abspath() for x in to_list(includes)]
    incs = ' '.join(incs)
    
    dict_src = name + 'Dict.cxx'

    bld(source = headers + [linkdef],
        target = dict_src,
        rule='${ROOTCINT} -f ${TGT} -c %s ${SRC}' % incs)

    bld.shlib(source = dict_src,
              target = name+'Dict',
              includes = includes,
              use = 'ROOTSYS')

    rootmap = 'lib%sDict.rootmap'%name
    bld(source = [linkdef], target=rootmap,
        rule='${RLIBMAP} -o ${TGT} -l lib%sDict.so -d lib%s.so -c ${SRC}' % (name, name))
    
    bld.install_files('${PREFIX}/lib/', rootmap)
                      
