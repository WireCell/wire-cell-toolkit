#!/usr/bin/env python

import os


TOP = '.'
APPNAME = 'WireCell'

def options(opt):
    opt.load("wcb")

    # this used in cfg/wscript_build
    opt.add_option('--install-config', type=str, default="",
                   help="Install configuration files for given experiment")


def configure(cfg):
    cfg.load("wcb")

    # fixme: should go into wcb.py
    cfg.find_program("jsonnet", var='JSONNET')

    # not meant for users: re-generate code
    cfg.find_program("moo", var='MOO', mandatory=False)

    # boost 1.59 uses auto_ptr and GCC 5 deprecates it vociferously.
    cfg.env.CXXFLAGS += ['-Wno-deprecated-declarations']
    cfg.env.CXXFLAGS += ['-Wall', '-Wno-unused-local-typedefs', '-Wno-unused-function']
    # cfg.env.CXXFLAGS += ['-Wpedantic', '-Werror']


from collections import defaultdict

def build(bld):

    if "MOO" in bld.env:
        for schema in bld.srcnode.ant_glob("*/schema/*.schema"):
            if schema.name.startswith("build"):
                continue
            do_schema(bld, schema)
        bld.add_group()         # make sure this is done first

    bld.load('wcb')

def do_schema(bld, schema):

    schema_path = schema.abspath()
    print("generating schema:",schema_path)
    Pkg, Comp = os.path.splitext(os.path.basename(schema_path))[0].split("_")

    params=dict(Pkg=Pkg, pkg=Pkg.lower(), Comp=Comp, top=bld.srcnode.abspath())
    inc_cfg=bld.srcnode.make_node('{pkg}/inc/WireCell{Pkg}/Cfg/{Comp}'.format(**params))


    src = [schema]
    if Pkg != "Util":           # likely a dependency
        src.append(bld.path.find_resource("util/schema/Util_Base.schema"))

    out = inc_cfg.make_node("Structs.hpp")
    bld(rule="""${{MOO}} -g /lang:ocpp.jsonnet \
    -M {top}/cfg -M {top}/util/schema \
    -A path="WireCell{Pkg}.Cfg.{Comp}" \
    -A os="${{SRC[0]}}" \
    render omodel.jsonnet ostructs.hpp.j2 \
    > ${{TGT[0]}}""".format(**params),
        source=src, target=[out])

    out = inc_cfg.make_node("Nljs.hpp")
    bld(rule="""${{MOO}} -g /lang:ocpp.jsonnet \
    -M {top}/cfg -M {top}/util/schema \
    -A path="WireCell{Pkg}.Cfg.{Comp}" \
    -A os="${{SRC[0]}}" \
    render omodel.jsonnet onljs.hpp.j2 \
    > ${{TGT[0]}}""".format(**params),
        source=src, target=[out])


    out='cfg/schema/{pkg}/{Comp}.jsonnet'.format(**params)
    out = bld.srcnode.make_node(out)
    bld(rule="""${{MOO}} \
    -T {top}/util/schema -M {top}/cfg -M {top}/util/schema \
    -A path="WireCell{Pkg}.Cfg.{Comp}" \
    -A os="${{SRC[0]}}" \
    render omodel.jsonnet wct-cfg-ctor.jsonnet.j2 \
    > ${{TGT[0]}}""".format(**params),
        source=src, target=[out])

