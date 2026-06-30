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

Each external dependency has a `WITH_<TOKEN>` option that mirrors waf's
`--with-<name>`:

* `WITH_<TOKEN>=` (empty) — auto: mandatory deps are required, ordinary
  optional deps are probed quietly.
* `WITH_<TOKEN>=no` — do not use (error if the dependency is mandatory).
* `WITH_<TOKEN>=yes` — require it.
* `WITH_<TOKEN>=/some/prefix` — use that install prefix as a search hint, then
  require it.

**Opt-in dependencies** — `WITH_ROOTSYS`, `WITH_CUDA`, `WITH_LIBTORCH`,
`WITH_KOKKOS` — are *never* auto-detected (matching waf's `with_p()` gating);
they are used only when their `WITH_*` is non-empty.  Their packages
(`root`, `cuda`, `pytorch`, ...) are gated on the corresponding dependency
being found.

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
