# CMake ⇄ waf parity harness

Files here validate that the CMake build (see `../README.md`) stays at feature
parity with the waf (`./wcb`) build.  Nothing in this directory is processed by
the CMake build itself: the top-level `CMakeLists.txt` only `add_subdirectory()`s
package directories and only `include()`s named modules from `cmake/`, so
`cmake/test/` is inert to CMake.

## `parity.sh`

Configures, builds and installs WCT with **both** build systems into separate
prefixes and diffs the installed results:

| Check | Severity | What it compares |
|---|---|---|
| `libraries` | hard | `lib/libWireCell*.so`, `lib/libWCP*.so` basenames |
| `applications` | hard | WCT executables in `bin/` |
| `header-packages` | soft | `include/WireCell*`, `include/WCP*` dirs |
| `BuildConfig.h` | hard | `WireCellUtil/BuildConfig.h`, normalized (version-independent) |
| `pc-Libs` / `pc-Requires` | hard | `wire-cell-toolkit.pc` `Libs:` and `Requires:` set |
| `cmake-targets` | hard | the `WireCell::*` imported targets the package config defines |
| test outcomes | info | `ctest` vs `./wcb --tests` (only with `WCT_RUN_TESTS=1`) |

Exit status is nonzero if any *hard* check differs; per-check diffs are written
under the work dir.

### Run it

Full run (builds both sides fresh into clean prefixes — the intended mode):

```bash
WCT_DEPS=$PWD/local \
WCT_CMAKE_ARGS="-DWITH_ROOTSYS=yes" \
WCT_WAF_ARGS="--with-root=$PWD/local" \
  bash toolkit/cmake/test/parity.sh
```

Compare-only (reuse existing installs, e.g. against a waf install prefix):

```bash
WCT_WAF_PREFIX=$PWD/local \
WCT_CMAKE_PREFIX=/path/to/cmake-install \
  bash toolkit/cmake/test/parity.sh
```

Ninja and tests:

```bash
WCT_DEPS=$PWD/local WCT_GENERATOR=Ninja WCT_RUN_TESTS=1 \
  bash toolkit/cmake/test/parity.sh
```

All knobs (`WCT_SRC`, `WCT_DEPS`, `WCT_WORK`, `WCT_GENERATOR`, `WCT_CMAKE_ARGS`,
`WCT_WAF_ARGS`, `WCT_JOBS`, `WCT_RUN_TESTS`, `WCT_WAF_PREFIX`,
`WCT_CMAKE_PREFIX`) are documented in the script header.

### Interpreting results

For a faithful comparison, build **both** sides on the **same machine** with the
**same dependency prefix and the same optional set** (e.g. both with ROOT, or
neither).  Two differences are expected and not build-system defects:

* **Version string** — `WIRECELL_VERSION`/`_DEVELOPMENT`/`_RELEASE` reflect the
  checkout; the `BuildConfig.h` check normalizes these away.
* **Environment probes** — macros like `HAVE_BACKTRACE_LIB` depend on whether a
  header/library is present in the toolchain at build time.  Comparing a fresh
  build against a waf tree built in a *different* environment can therefore flag
  `HAVE_BACKTRACE_LIB`; a same-machine run of both agrees.

## `ci-parity.yml`

A GitHub Actions workflow template that builds WCT with CMake (configure →
build → ctest) and with waf, then runs `parity.sh`.  It assumes a container
image that already provides the WCT dependencies (Boost, Eigen, spdlog, FFTW,
JsonCpp, Jsonnet, ROOT, ...).  To activate it, copy it to
`.github/workflows/` and set the `container.image` (and any `WCT_*` args) to
match your dependency image:

```bash
mkdir -p .github/workflows
cp toolkit/cmake/test/ci-parity.yml .github/workflows/
```
