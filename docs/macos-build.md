# Building WCT natively on macOS

macOS is not an officially supported platform (CI runs on Linux), but the
toolkit builds and runs natively on macOS with Homebrew-provided
dependencies.  This recipe was validated on macOS 15 (Apple Silicon) with
Apple clang 21 and Homebrew Boost 1.90.

The differences from the standard Linux build (see the top-level
[README](../README.org)) come down to: one dependency built from source
(C++ libjsonnet), a few extra compiler/linker flags, explicit dependency
paths at configure time, and a runtime environment variable for plugin
loading.

## Prerequisites

Xcode Command Line Tools and cmake, plus:

```sh
brew install boost eigen spdlog jsoncpp fftw
# fmt is pulled in as a spdlog dependency
# tbb, root, nlohmann-json are optional
```

### C++ libjsonnet (built from source)

Homebrew's `jsonnet` package is go-jsonnet: a command-line binary with no
linkable library.  WCT links the C++ `libjsonnet`, so build it once:

```sh
git clone --depth 1 https://github.com/google/jsonnet ~/opt/jsonnet-src
cd ~/opt/jsonnet-src
cmake -S . -B build -DBUILD_SHARED_BINARIES=ON -DBUILD_TESTS=OFF \
      -DBUILD_JSONNET=OFF -DBUILD_JSONNETFMT=OFF \
      -DCMAKE_INSTALL_PREFIX=$HOME/opt/jsonnet \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j4 && cmake --install build
```

This installs `~/opt/jsonnet/include/libjsonnet.h` and
`~/opt/jsonnet/lib/libjsonnet.dylib`.

## Configure

Two macOS-specific concerns are handled through environment variables,
exported before running configure (they are baked into the configure
cache, so only re-run configure to change them):

- Homebrew's spdlog is compiled against an *external* fmt, so fmt's
  include and link paths must be supplied.
- Boost.Stacktrace needs `BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED` since
  macOS has no glibc.

```sh
export CXXFLAGS="-I$(brew --prefix fmt)/include -DBOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED"
export LINKFLAGS="-L$(brew --prefix fmt)/lib -lfmt"
export LDFLAGS="-L$(brew --prefix fmt)/lib"
```

Homebrew kegs are not in the default header/library search paths, so each
dependency is given explicitly:

```sh
./wcb configure \
  --prefix=$HOME/opt/wct \
  --with-root=no --with-tbb=no \
  --boost-includes=$(brew --prefix boost)/include \
  --boost-libs=$(brew --prefix boost)/lib \
  --with-eigen-include=$(brew --prefix eigen)/include/eigen3 \
  --with-spdlog=$(brew --prefix spdlog) \
  --with-jsoncpp=$(brew --prefix jsoncpp) \
  --with-jsonnet=$HOME/opt/jsonnet \
  --with-fftw=$(brew --prefix fftw)
```

If you have ROOT or TBB installed via brew the corresponding packages can
be enabled in the usual way.

## Build and install

```sh
./wcb -j8
./wcb install
```

## Runtime environment

```sh
export PATH=$HOME/opt/wct/bin:$PATH
export DYLD_FALLBACK_LIBRARY_PATH=$HOME/opt/wct/lib
```

The `DYLD_FALLBACK_LIBRARY_PATH` is needed because plugin libraries are
`dlopen()`ed by bare name (eg `libWireCellGen.dylib`) and macOS provides
no `ld.so.conf`-style search path configuration.

Check that the build works, including plugin/factory loading:

```sh
wire-cell --version
wire-cell -p WireCellGen --help
```

### SIP caveat

System Integrity Protection strips `DYLD_*` variables from the
environment whenever a process is launched via a SIP-protected binary
(`/usr/bin/perl`, `/usr/bin/timeout`, system shells in some contexts).
If plugin loading mysteriously fails under a test harness or wrapper
script, run `wire-cell` directly.
