# Building Wire-Cell Toolkit with CMake

This directory holds a modern-CMake build for the Wire-Cell Toolkit that runs
**in parallel to the waf (`./wcb`) build**.  It does not replace waf and makes
no changes to `wscript`, `*/wscript_build` or `waft/*`.

The commands below assume you work from a directory that contains:

```
.
├── toolkit/   # the WCT source (this checkout; holds the top-level CMakeLists.txt)
├── build/     # out-of-tree CMake build products (created for you)
└── local/     # install prefix of the external dependencies (Boost, ROOT, ...)
```

`local/` is the prefix where WCT's dependencies are installed; it is passed to
CMake as `CMAKE_PREFIX_PATH` so `find_package`/`pkg-config` locate them.  Keep
the CMake `build/` directory distinct from any `toolkit/build/` used by waf.

> Throughout, `"$PWD/local"` is the dependency prefix.  Use an absolute path
> (CMake resolves relative `CMAKE_PREFIX_PATH` against each project, not your
> shell's CWD).

## 1. Configure

### Unix Makefiles (default generator)

```bash
cmake -S toolkit -B build -DCMAKE_PREFIX_PATH="$PWD/local"
```

### Ninja

```bash
cmake -S toolkit -B build -G Ninja -DCMAKE_PREFIX_PATH="$PWD/local"
```

The configure step prints the detected version/build-mode, the C++ standard,
and which dependencies and packages were found, e.g.:

```
-- Wire-Cell Toolkit 0.35.0-337-g0447d488 [triplet 0.35.0]
--   build mode      : release (strict warnings)
-- WCT dependencies found  : SPDLOG;BOOST;FFTW;EIGEN;...
-- WCT packages configured: util;iface;aux;gen;...;tbb;hio
```

## 2. Build

```bash
cmake --build build -j            # generator-agnostic (Make or Ninja)
# or, with the Ninja build tree:
ninja -C build
```

Build a single target (library, app, or test):

```bash
cmake --build build -j --target WireCellUtil
cmake --build build -j --target wire-cell
```

Shared libraries are collected under `build/lib/`; applications
(`wire-cell`, `wcsonnet`, `wcwires`, ...) under their package build dirs.

## 3. Common options

Pass these as `-D<NAME>=<VALUE>` at configure time.

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_PREFIX_PATH` | – | Where to find dependencies (point at `local/`). |
| `CMAKE_INSTALL_PREFIX` | `/usr/local` | Install destination (§6). |
| `CMAKE_BUILD_TYPE` | `RelWithDebInfo` | `Debug`/`Release`/`RelWithDebInfo`/`MinSizeRel` (≈ waf `-O2 -ggdb3`). |
| `WCT_BUILD_MODE` | `` (auto) | Force `development` (relaxed warnings) or `release` (strict, `-Werror`). Auto-detected from the git branch otherwise. |
| `WCT_CXX_STANDARD` | `17` | C++ standard (`17`/`20`/`23`); mirrors `--cxxstd`. |
| `WCT_SPDLOG_ACTIVE_LEVEL` | `trace` | Compiled-in minimum spdlog level (`trace`…`off`). |
| `WCT_SPDLOG_STATIC` | `ON` | Define `SPDLOG_COMPILED_LIB`. |
| `WCT_WITH_TESTS` | `OFF` | Build and register the test suite (§5); mirrors `--tests`. |
| `WCT_INSTALL_CONFIG` | `` | Comma list of experiments (or `all`) whose `cfg/pgrapher/experiment/*` configs to install; mirrors `--install-config`. |

## 4. Dependencies

Each external dependency `<TOKEN>` (`SPDLOG`, `BOOST`, `FFTW`, `EIGEN`,
`JSONCPP`, `JSONNET`, `ZLIB`, `BZIP2`, `TBB`, `HDF5`, `ROOTSYS`, ...) has a
family of `WITH_*` cache options that mirror waf's
`--with-NAME[/-include/-lib/-libs]`:

| Option | waf equivalent | Meaning |
|---|---|---|
| `WITH_<TOKEN>` | `--with-name` | `` (auto) / `no` / `yes` / an **install-prefix path** |
| `WITH_<TOKEN>_INCLUDE` | `--with-name-include` | include dir(s) (comma list), overrides `<prefix>/include` |
| `WITH_<TOKEN>_LIB` | `--with-name-lib` | library dir(s) (comma list), overrides `<prefix>/lib` |
| `WITH_<TOKEN>_LIBS` | `--with-name-libs` | exact library name(s) (comma list), overrides the built-in default |

`WITH_<TOKEN>` alone behaves like waf: empty means auto (mandatory deps are
required, ordinary optional deps are probed quietly); `no` disables; `yes`
requires; a path is used as an install-prefix search hint and requires the dep.

### Discovering every `WITH_*` option

The options are declared during configuration, so configure once, then list
them (with their one-line help) from the cache:

```bash
cmake -S toolkit -B build -DCMAKE_PREFIX_PATH="$PWD/local"
cmake -LH build | grep -B1 '^WITH_'
```

`cmake -LAH build` shows all cache variables (including advanced ones); or use
the interactive `ccmake build` / `cmake-gui build`.

### Boost

Boost is treated **uniformly**, not special-cased — there is no separate set of
`--boost-*` flags as in waf.  Use `WITH_BOOST=/prefix` (or
`WITH_BOOST_INCLUDE` / `WITH_BOOST_LIB`) to point at a non-standard Boost.
Modern Boost ships `BoostConfig.cmake`, which is preferred automatically, so the
multithreaded/ABI selection waf did with `--boost-mt` is handled by CMake.

### Fine-grained include/lib locations and exact library names

`WITH_<TOKEN>_INCLUDE` and `WITH_<TOKEN>_LIB` add explicit search dirs (fed to
`find_package`, `find_library`/`find_path`, and pkg-config), and
`WITH_<TOKEN>_LIBS` overrides the library name(s) to link.  For example, to link
the Go implementation of Jsonnet (`libgojsonnet.so`) instead of the C one — the
CMake equivalent of waf's `--with-jsonnet-libs=gojsonnet`:

```bash
cmake -S toolkit -B build -DCMAKE_PREFIX_PATH="$PWD/local" \
      -DWITH_JSONNET_LIBS=gojsonnet
# or point at a specific tree and library:
cmake -S toolkit -B build \
      -DWITH_JSONNET=/opt/jsonnet \
      -DWITH_JSONNET_INCLUDE=/opt/jsonnet/include \
      -DWITH_JSONNET_LIB=/opt/jsonnet/lib \
      -DWITH_JSONNET_LIBS=gojsonnet
```

### Opt-in dependencies

`WITH_ROOTSYS`, `WITH_CUDA`, `WITH_LIBTORCH`, `WITH_KOKKOS` are *never*
auto-detected (matching waf's `with_p()` gating); they are used only when their
`WITH_*` is non-empty.  Their packages (`root`, `cuda`, `pytorch`, ...) are
gated on the corresponding dependency being found.

```bash
# Build with ROOT and CUDA enabled (deps live under local/):
cmake -S toolkit -B build -G Ninja \
      -DCMAKE_PREFIX_PATH="$PWD/local" \
      -DWITH_ROOTSYS=yes -DWITH_CUDA=yes
```

A package whose optional dependency is absent is skipped, and configure says so:

```
-- WCT package 'zio' disabled (missing dependency: ZMQ;CZMQ;ZYRE;ZIO)
```

## 5. Tests

Tests are opt-in.  Configure with `-DWCT_WITH_TESTS=ON`, build, then use CTest:

```bash
cmake -S toolkit -B build -G Ninja -DCMAKE_PREFIX_PATH="$PWD/local" -DWCT_WITH_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Tests carry **labels** that map to the waf test groups; select with `-L`:

```bash
ctest --test-dir build -L atomic       # the test_*/atomic* unit tests
ctest --test-dir build -L doctest      # the aggregated wcdoctest-<pkg> runners
ctest --test-dir build -L script       # interpreted (.py/.sh/.bats/.jsonnet) tests
ctest --test-dir build -L history      # history-group tests (need test data, below)
ctest --test-dir build -R WireCellUtil # by name regex
ctest --test-dir build -N              # list without running
```

`check_*.cxx` programs are built but not run (waf's "check" group).  Run a
compiled test under a wrapper (mirrors `--testcmd`):

```bash
cmake -S toolkit -B build -DWCT_WITH_TESTS=ON \
      -DWCT_TEST_LAUNCHER="valgrind;--error-exitcode=1" ...
```

### Test data

History/report and some atomic tests read a downloaded data repository.  Fetch
it into `build/tests/` (mirrors waf's data-repo handling):

```bash
cmake --build build --target test-data        # download + unpack
cmake --build build --target pack-test-data   # re-archive (≈ ./wcb packrepo)
```

Override the source/versions at configure time with `WCT_TEST_DATA_URL` and
`WCT_TEST_DATA_VERSIONS`.  A CTest fixture `WCTTestData` is also available for
on-demand provisioning of data-dependent tests.

## 6. Install

```bash
cmake -S toolkit -B build -G Ninja \
      -DCMAKE_PREFIX_PATH="$PWD/local" \
      -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build -j
cmake --install build
```

This installs libraries to `<prefix>/lib`, headers to `<prefix>/include`
(including the generated `WireCellUtil/BuildConfig.h`), apps to `<prefix>/bin`,
Jsonnet configuration to `<prefix>/share/wirecell`, the pkg-config file
`<prefix>/lib/pkgconfig/wire-cell-toolkit.pc`, and the CMake package config
under `<prefix>/lib/cmake/WireCellToolkit/`.

## 7. Using WCT from a downstream project

CMake `find_package` (targets are namespaced `WireCell::`, with an aggregate
`WireCell::WireCell`):

```cmake
find_package(WireCellToolkit REQUIRED)
target_link_libraries(myapp PRIVATE WireCell::Util WireCell::Aux)
```

```bash
cmake -S mydir -B mybuild -DCMAKE_PREFIX_PATH="$PWD/install;$PWD/local"
```

(The WCT install prefix and the dependency prefix are both on
`CMAKE_PREFIX_PATH` so the config can re-find WCT's external dependencies.)

Or via pkg-config:

```bash
PKG_CONFIG_PATH="$PWD/install/lib/pkgconfig:$PWD/local/lib/pkgconfig" \
  pkg-config --cflags --libs wire-cell-toolkit
```

## 8. Reconfigure / clean

```bash
rm -rf build                       # full clean: just delete the build tree
cmake --build build --target clean # remove build products, keep the cache
```

Source files are discovered by globbing (`CONFIGURE_DEPENDS`), so adding or
removing a `.cxx`/`.h` is picked up on the next build without editing any
CMake file.  Changing a package's dependency tokens does require editing its
`CMakeLists.txt` (the one piece of per-package metadata duplicated with the
waf `wscript_build`).

## 9. Maintaining the CMake build alongside waf

The CMake build is **additive** and runs in parallel with waf; keeping both is
well under 2× the effort of one because almost everything that changes often is
either shared or convention-driven.  Only two kinds of metadata are duplicated:

1. **Per-package dependency tokens.** A package's `use=`/`app_use=`/`test_use=`
   in `*/wscript_build` is mirrored by `wct_package(... USE ... APP_USE ...
   TEST_USE ...)` in `*/CMakeLists.txt`.  Adding an inter-package edge or an
   external dependency to a package means editing both one-liners.

2. **The external-dependency recipe.** waf's `waft/wcb.py:package_descriptions`
   + `waft/generic.py` are mirrored by `cmake/WCTDependencies.cmake` (discovery)
   and the small re-find table in `cmake/WCTExport.cmake` (so downstream
   consumers inherit the dep).  Adding a *new* external dependency touches both
   build systems here.

Everything else is **not** duplicated: adding/removing a source, header, app or
test file is picked up by globbing in both builds; the C++ code, the install
layout, `BuildConfig.h` contents, and the optional-dependency gating rules have
a single conceptual definition each.

**Adding a package** — create the directory with `src/`, `inc/<Name>/`,
`test/`, add a one-line `wscript_build` and a one-line `CMakeLists.txt`
(`wct_package(<Name> USE ...)`), and add its directory to the package list in
the top-level `CMakeLists.txt` (or the optional-gating table if it depends on an
optional external).

**Adding an external dependency** — add its `package_descriptions` entry (waf)
and a discovery block in `cmake/WCTDependencies.cmake` (a `WCT::<TOKEN>` target
+ `WCT_HAVE_<TOKEN>` flag); if consumers need it transitively, add a one-line
case to `_wct_refind` in `cmake/WCTExport.cmake`; if it should gate a package,
add a row to the optional table in the top-level `CMakeLists.txt`.

**Guarding against drift** — the parity harness under `cmake/test/` builds and
installs both ways and diffs the results (libraries, apps, headers,
`BuildConfig.h`, `.pc`, and the CMake target set).  Run it before a release, or
wire it into CI (`cmake/test/ci-parity.yml`):

```bash
WCT_DEPS=$PWD/local bash toolkit/cmake/test/parity.sh
```

See `cmake/test/README.md` for the harness details and the expected
version/environment-probe caveats.
