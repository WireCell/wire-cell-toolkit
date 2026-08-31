# -*- python -*-
'''Generate CMake package-config files for wire-cell-toolkit.

Although WCT builds with waf, this tool lets downstream projects that build
with CMake depend on an installed WCT using the conventional

    find_package(WireCellToolkit REQUIRED)
    target_link_libraries(myapp PRIVATE WireCell::Util WireCell::Aux ...)

The files are generated at `./wcb install` time from information already
gathered by the waf build and are installed into the conventional location
${LIBDIR}/cmake/WireCellToolkit/:

 - WireCellToolkitConfig.cmake         entry point for find_package()
 - WireCellToolkitConfigVersion.cmake  version compatibility check

Data sources (all populated by waft/smplpkgs.py and waft/wcb.py):

 - bld.smplpkg_libs    names of WCT packages that build a shared library
 - bld.smplpkg_graph   inter-library "use" dependency graph
 - bld.env             VERSION, cxxshlib_PATTERN, per-dependency results

The config files re-discover WCT's own external dependencies (Boost, Eigen,
...) so a downstream consumer inherits them transitively and need not
find_package(Boost) etc. itself.  Mandatory dependencies (needed by every
build) are found with find_dependency() and hard-fail if missing.  Optional
dependencies are found quietly and each WCT imported target is guarded so it is
only defined when the dependencies it (transitively) needs were found -- a
consumer that uses only the core libraries therefore does not need ROOT, HDF5,
etc. just because the local build happened to include them.

See https://github.com/WireCell/wire-cell-toolkit/issues/484
'''

import re
from waflib.Logs import debug, info, warn
from waflib.Utils import to_list

PACKAGE = "WireCellToolkit"
NAMESPACE = "WireCell"

# Prefixes stripped from a WCT library name to form its CMake target name,
# e.g. WireCellUtil -> WireCell::Util, WCPQuickhull -> WireCell::Quickhull.
LIB_PREFIXES = ("WireCell", "WCP")

# The "public" libraries gathered into the convenience WireCell::WireCell
# aggregate target (mirrors the libs exported by the pkg-config .pc file).
# All are mandatory so the aggregate is never guarded.
AGGREGATE_LIBS = ("WireCellAux", "WireCellIface", "WireCellUtil")


def _cmake(package, libs, components=(), includes=(), found=None, required=False):
    'A dependency re-found with a native CMake package (config or find-module).'
    return dict(kind="cmake", package=package, components=list(components),
                libs=list(libs), includes=list(includes),
                found=found or (package + "_FOUND"), required=required)


def _pc(tag, pcname, required=False):
    'A dependency re-found via pkg-config (CMake PkgConfig module).'
    var = "WCT_PC_" + tag
    return dict(kind="pkgconfig", var=var, pcname=pcname,
                libs=["PkgConfig::" + var], includes=[],
                found=var + "_FOUND", required=required)


def _findlib(tag, lib, header=None, required=False):
    'A dependency with neither a CMake package nor a usable .pc.'
    libvar = "WCT_LIB_" + tag
    includes = []
    incvar = None
    if header:
        incvar = "WCT_INC_" + tag
        includes = ["${%s}" % incvar]
    return dict(kind="findlib", libvar=libvar, lib=lib, header=header,
                incvar=incvar, libs=["${%s}" % libvar], includes=includes,
                found=libvar, required=required)


# Map a WCT external "use" token (upper-cased) to how a downstream CMake
# project re-discovers it.  Tokens absent here fall back to find_library().
#
# Mandatory tokens (required=True) are needed by every build and hard-fail if
# missing; the rest are optional and gate the WCT targets that use them.  The
# pcname values mirror waft/wcb.py package_descriptions.
DEPMAP = {
    # --- mandatory: needed by WireCellUtil and the core chain ---
    "BOOST": _cmake("Boost",
                    ["Boost::filesystem", "Boost::graph", "Boost::thread",
                     "Boost::program_options", "Boost::iostreams", "Boost::regex"],
                    components=["filesystem", "graph", "thread",
                                "program_options", "iostreams", "regex"],
                    required=True),
    "EIGEN":   _cmake("Eigen3", ["Eigen3::Eigen"], found="Eigen3_FOUND", required=True),
    "SPDLOG":  _cmake("spdlog", ["spdlog::spdlog"], found="spdlog_FOUND", required=True),
    "ZLIB":    _cmake("ZLIB", ["ZLIB::ZLIB"], required=True),
    "BZIP2":   _cmake("BZip2", ["BZip2::BZip2"], found="BZip2_FOUND", required=True),
    "FFTW":    _pc("FFTW", "fftw3f", required=True),
    "JSONCPP": _pc("JSONCPP", "jsoncpp", required=True),
    "JSONNET": _findlib("JSONNET", "jsonnet", "libjsonnet.h", required=True),
    "DYNAMO":  dict(kind="builtin", libs=["${CMAKE_DL_LIBS}"], includes=[],
                    found=None, required=True),

    # --- optional ---
    "TBB":      _cmake("TBB", ["TBB::tbb"], found="TBB_FOUND"),
    "PROTOBUF": _cmake("Protobuf", ["protobuf::libprotobuf"], found="Protobuf_FOUND"),
    "GRPC":     _cmake("gRPC", ["gRPC::grpc++", "gRPC::grpc", "gRPC::gpr"], found="gRPC_FOUND"),
    "PYTHON":   _cmake("Python3", ["Python3::Python"],
                       components=["Development.Embed"], found="Python3_FOUND"),
    "ROOTSYS":  _cmake("ROOT", ["${ROOT_LIBRARIES}"], includes=["${ROOT_INCLUDE_DIRS}"],
                       found="ROOT_FOUND"),
    "LIBTORCH": _cmake("Torch", ["${TORCH_LIBRARIES}"], includes=["${TORCH_INCLUDE_DIRS}"],
                       found="Torch_FOUND"),
    "CUDA":     _cmake("CUDAToolkit", ["CUDA::cudart"], found="CUDAToolkit_FOUND"),
    "PTHREAD":  _cmake("Threads", ["Threads::Threads"], found="Threads_FOUND"),
    # HDF5: CMake's FindHDF5 is unreliable; WCT resolves it via the hdf5 .pc.
    "HDF5":     _pc("HDF5", "hdf5"),
    "ZMQ":      _pc("ZMQ", "libzmq"),
    "CZMQ":     _pc("CZMQ", "libczmq"),
    "ZYRE":     _pc("ZYRE", "libzyre"),
    "ZIO":      _pc("ZIO", "libzio"),
    "GLPK":        _findlib("GLPK", "glpk", "glpk.h"),
    "FFTWTHREADS": _findlib("FFTWTHREADS", "fftw3f_threads"),
    "TRITON":      _findlib("TRITON", "grpcclient", "grpc_client.h"),
}


def _target_name(libname):
    'Map a WCT library name to its namespaced CMake target name.'
    for pre in LIB_PREFIXES:
        if libname.startswith(pre):
            return "%s::%s" % (NAMESPACE, libname[len(pre):])
    return "%s::%s" % (NAMESPACE, libname)


def _is_internal(token):
    return token.startswith(LIB_PREFIXES)


def _lib_deps(graph, pkg):
    'Return the "use" (library-interface) dependency tokens of a package.'
    deps = []
    for (a, b), attrs in graph._edges.items():
        if a == pkg and 'lib' in attrs:
            deps.append(b)
    deps.sort()
    return deps


def _toposort(libs, graph):
    'Order libraries so internal dependencies precede their dependents.'
    libset = set(libs)
    ordered, seen = [], set()

    def visit(name):
        if name in seen or name not in libset:
            return
        seen.add(name)
        for dep in _lib_deps(graph, name):
            if dep in libset:
                visit(dep)
        ordered.append(name)

    for name in sorted(libs):
        visit(name)
    return ordered


def _spec_for(token, unknown):
    'Return (UPPER token, normalized spec) for an external dependency token.'
    up = token.upper()
    spec = DEPMAP.get(up)
    if spec is None:
        unknown.add(token)
        spec = _findlib(re.sub(r'[^A-Za-z0-9_]', '_', up), token.lower())
    return up, spec


def _hints(paths, prefix_subdir):
    'Render a find_library/find_path HINTS clause from discovered paths.'
    hints = ['"%s"' % p for p in paths]
    hints.append('"${WireCellToolkit_PREFIX}/%s"' % prefix_subdir)
    return " HINTS " + " ".join(hints)


def _resolve(up, spec, env):
    '''Resolve a dependency spec against the configured build environment.

    Returns dict(find, libs, includes, found, needs_pc).  find_library /
    find_path fallbacks take their library names and search hints from what
    waf actually discovered (env LIB_*/LIBPATH_*/INCLUDES_*), so they reflect
    this build (e.g. gojsonnet vs jsonnet) and resolve in the install prefix.
    '''
    kind = spec["kind"]
    required = spec.get("required", False)

    if kind == "builtin":
        return dict(find=[], libs=spec["libs"], includes=spec["includes"],
                    found=None, needs_pc=False)

    if kind == "cmake":
        args = spec["package"]
        if spec["components"]:
            args += " COMPONENTS " + " ".join(spec["components"])
        find = ["find_dependency(%s)" % args] if required \
            else ["find_package(%s QUIET)" % args]
        return dict(find=find, libs=spec["libs"], includes=spec["includes"],
                    found=spec["found"], needs_pc=False)

    if kind == "pkgconfig":
        req = "REQUIRED " if required else ""
        find = ["pkg_check_modules(%s %sIMPORTED_TARGET %s)"
                % (spec["var"], req, spec["pcname"])]
        return dict(find=find, libs=spec["libs"], includes=spec["includes"],
                    found=spec["found"], needs_pc=True)

    # findlib: derive names and hints from the configured build.
    tag = re.sub(r'[^A-Za-z0-9_]', '_', up)
    names = to_list(getattr(env, 'LIB_' + up, None) or spec.get("lib") or "")
    libhints = to_list(getattr(env, 'LIBPATH_' + up, None) or [])
    inchints = to_list(getattr(env, 'INCLUDES_' + up, None) or [])

    find, libs, found = [], [], None
    if len(names) <= 1:
        var = "WCT_LIB_" + tag
        find.append("find_library(%s NAMES %s%s)"
                    % (var, names[0] if names else tag.lower(),
                       _hints(libhints, "lib")))
        libs = ["${%s}" % var]
        found = var
    else:
        for i, nm in enumerate(names):
            var = "WCT_LIB_%s_%d" % (tag, i)
            find.append("find_library(%s NAMES %s%s)" % (var, nm, _hints(libhints, "lib")))
            libs.append("${%s}" % var)
            if i == 0:
                found = var

    includes = []
    if spec.get("header"):
        incvar = "WCT_INC_" + tag
        find.append("find_path(%s NAMES %s%s)"
                    % (incvar, spec["header"], _hints(inchints, "include")))
        includes = ["${%s}" % incvar]

    return dict(find=find, libs=libs, includes=includes, found=found, needs_pc=False)


def _version_triplet(version):
    'Extract a numeric X.Y.Z from a git-describe style version string.'
    m = re.search(r'(\d+)\.(\d+)\.(\d+)', version or "")
    if m:
        return "%s.%s.%s" % m.groups()
    m = re.search(r'(\d+)\.(\d+)', version or "")
    if m:
        return "%s.%s.0" % m.groups()
    return "0.0.0"


def _config_text(bld):
    'Build the text of WireCellToolkitConfig.cmake.'
    graph = bld.smplpkg_graph
    libs = _toposort(getattr(bld, "smplpkg_libs", []), graph)
    libset = set(libs)
    pattern = bld.env.cxxshlib_PATTERN or "lib%s.so"
    unknown = set()

    # Resolve each external token once against the build environment.
    resolved = {}

    def resolve(token):
        up, spec = _spec_for(token, unknown)
        if up not in resolved:
            resolved[up] = _resolve(up, spec, bld.env)
        return up, resolved[up]

    # Collect the unique external-dependency find commands in the order
    # dependencies are first encountered, and remember whether any need
    # pkg-config.
    find_lines, seen_tokens = [], set()
    needs_pkgconfig = False
    for name in libs:
        for dep in _lib_deps(graph, name):
            if _is_internal(dep):
                continue
            up, info = resolve(dep)
            if info["needs_pc"]:
                needs_pkgconfig = True
            if up in seen_tokens:
                continue
            seen_tokens.add(up)
            find_lines += info["find"]

    # For each library, accumulate the set of "found" guard variables it (and
    # its internal deps, transitively) need.  Topological order guarantees a
    # dependency's guards are known before its dependents'.
    guard = {}

    L = []
    L.append("# Generated by waf (waft/cmake.py) for wire-cell-toolkit -- do not edit.")
    L.append("# Lets CMake projects build against an installed wire-cell-toolkit:")
    L.append("#   find_package(WireCellToolkit REQUIRED)")
    L.append("#   target_link_libraries(myapp PRIVATE WireCell::Util WireCell::Aux)")
    L.append("")
    L.append("cmake_minimum_required(VERSION 3.16)")
    L.append("include(CMakeFindDependencyMacro)")
    L.append("")
    # The file installs to <prefix>/lib/cmake/WireCellToolkit so the install
    # prefix is three directories up: keep the package relocatable.
    L.append('get_filename_component(_wct_cmakedir "${CMAKE_CURRENT_LIST_DIR}" REALPATH)')
    L.append('get_filename_component(WireCellToolkit_PREFIX "${_wct_cmakedir}/../../.." REALPATH)')
    L.append('set(WireCellToolkit_INCLUDE_DIRS "${WireCellToolkit_PREFIX}/include")')
    L.append('set(WireCellToolkit_LIBRARY_DIR "${WireCellToolkit_PREFIX}/lib")')
    L.append("")
    if needs_pkgconfig:
        L.append("find_package(PkgConfig REQUIRED)")
    if find_lines:
        L.append("# Re-discover WCT's external dependencies so consumers inherit them.")
        L.append("# Mandatory ones hard-fail; optional ones gate the targets below.")
        L += find_lines
    L.append("")

    # Define one imported target per built library, in dependency order.
    for name in libs:
        tgt = _target_name(name)
        link, incs, conds = [], [], set()
        for dep in _lib_deps(graph, name):
            if _is_internal(dep):
                if dep in libset:
                    link.append(_target_name(dep))
                    conds |= guard.get(dep, set())
                else:
                    debug('cmake: %s skips unbuilt internal dep %s' % (name, dep))
                continue
            _up, info = resolve(dep)
            link += info["libs"]
            incs += info["includes"]
            _spec = DEPMAP.get(_up)
            if not (_spec and _spec.get("required")) and info["found"]:
                conds.add(info["found"])
        guard[name] = conds

        indent = "  " if conds else ""
        if conds:
            L.append("if(%s)" % " AND ".join(sorted(conds)))
        all_incs = ['${WireCellToolkit_INCLUDE_DIRS}'] + incs
        L.append('%sif(NOT TARGET %s)' % (indent, tgt))
        L.append('%s  add_library(%s SHARED IMPORTED)' % (indent, tgt))
        L.append('%s  set_target_properties(%s PROPERTIES' % (indent, tgt))
        L.append('%s    IMPORTED_LOCATION "${WireCellToolkit_LIBRARY_DIR}/%s"'
                 % (indent, pattern % name))
        L.append('%s    INTERFACE_INCLUDE_DIRECTORIES "%s"' % (indent, ";".join(all_incs)))
        L.append('%s    INTERFACE_COMPILE_FEATURES "cxx_std_17"' % indent)
        if link:
            L.append('%s    INTERFACE_LINK_LIBRARIES "%s"' % (indent, ";".join(link)))
        L.append('%s    )' % indent)
        L.append('%sendif()' % indent)
        if conds:
            L.append("endif()")
        L.append("")

    # Convenience aggregate target over the (mandatory) public libraries.
    agg = [_target_name(n) for n in AGGREGATE_LIBS if n in libset]
    if agg:
        L.append('if(NOT TARGET %s::%s)' % (NAMESPACE, NAMESPACE))
        L.append('  add_library(%s::%s INTERFACE IMPORTED)' % (NAMESPACE, NAMESPACE))
        L.append('  set_target_properties(%s::%s PROPERTIES' % (NAMESPACE, NAMESPACE))
        L.append('    INTERFACE_LINK_LIBRARIES "%s")' % ";".join(agg))
        L.append('endif()')
        L.append("")

    L.append('set(WireCellToolkit_LIBRARIES %s::%s)' % (NAMESPACE, NAMESPACE))
    L.append('set(WireCellToolkit_FOUND TRUE)')
    L.append("")

    if unknown:
        warn('cmake: no dependency mapping for %s, used find_library fallback'
             % ', '.join(sorted(unknown)))

    return "\n".join(L)


def _version_text(bld):
    'Build the text of WireCellToolkitConfigVersion.cmake.'
    full = bld.env.VERSION or ""
    triplet = _version_triplet(full)
    L = []
    L.append("# Generated by waf (waft/cmake.py) for wire-cell-toolkit -- do not edit.")
    L.append('set(PACKAGE_VERSION "%s")' % triplet)
    L.append('set(WireCellToolkit_VERSION_STRING "%s")' % full)
    L.append("")
    L.append("if(PACKAGE_VERSION VERSION_LESS PACKAGE_FIND_VERSION)")
    L.append("  set(PACKAGE_VERSION_COMPATIBLE FALSE)")
    L.append("else()")
    L.append("  set(PACKAGE_VERSION_COMPATIBLE TRUE)")
    L.append("  if(PACKAGE_FIND_VERSION STREQUAL PACKAGE_VERSION)")
    L.append("    set(PACKAGE_VERSION_EXACT TRUE)")
    L.append("  endif()")
    L.append("endif()")
    return "\n".join(L)


def write_cmake_config(bld):
    '''Register tasks to generate and install the CMake package-config files.

    Call from a build() after bld.recurse() so the smplpkg graph is complete.
    '''
    if not hasattr(bld, "smplpkg_graph") or not getattr(bld, "smplpkg_libs", None):
        debug("cmake: no smplpkg graph/libs, skipping CMake config generation")
        return

    cfgnode = bld.path.find_or_declare(PACKAGE + "Config.cmake")
    vernode = bld.path.find_or_declare(PACKAGE + "ConfigVersion.cmake")

    # Write the files eagerly during build-graph construction.  A waf rule task
    # would be signed by the rule function's bytecode (identical regardless of
    # the generated text), so it would not regenerate when the content changes.
    cfgnode.parent.mkdir()
    cfgnode.write(_config_text(bld))
    vernode.write(_version_text(bld))

    bld.install_files('${LIBDIR}/cmake/' + PACKAGE, [cfgnode, vernode])
    debug("cmake: will install %s package config to ${LIBDIR}/cmake/%s"
          % (PACKAGE, PACKAGE))
