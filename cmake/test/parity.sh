#!/usr/bin/env bash
#
# parity.sh -- waf vs CMake build parity harness for Wire-Cell Toolkit.
# (epic wct-ike, task wct-ike.13)
#
# Configures, builds and installs WCT with BOTH build systems into separate
# prefixes and diffs the results:
#   - installed WCT file lists (libraries, headers, apps, share, pkgconfig, cmake)
#   - the generated WireCellUtil/BuildConfig.h  (normalized: version-independent)
#   - wire-cell-toolkit.pc                       (Libs + Requires set)
#   - the installed CMake package config          (the WireCell:: target set)
#   - optionally, the test outcomes              (WCT_RUN_TESTS=1)
#
# Living under cmake/test/, this file is never processed by the CMake build
# (nothing add_subdirectory()s cmake/, and include() only pulls named modules
# from cmake/).  Run it by hand or from CI (see ci-parity.yml).
#
# Configuration (all via environment):
#   WCT_SRC          toolkit source dir            [default: this repo]
#   WCT_DEPS         dependency prefix (=CMAKE_PREFIX_PATH / waf --with hints)
#   WCT_WORK         scratch dir                   [default: mktemp -d]
#   WCT_GENERATOR    CMake generator              [default: "Unix Makefiles"; try "Ninja"]
#   WCT_CMAKE_ARGS   extra cmake -D args           (e.g. "-DWITH_ROOTSYS=yes")
#   WCT_WAF_ARGS     extra ./wcb configure args    (e.g. "--with-root=$PREFIX")
#   WCT_JOBS         parallel build jobs           [default: nproc]
#   WCT_RUN_TESTS    1 to also compare test runs   [default: 0]
#
# Compare-only mode (skip building one/both side): point at existing installs
#   WCT_WAF_PREFIX   use this existing waf install prefix (skip waf build)
#   WCT_CMAKE_PREFIX use this existing cmake install prefix (skip cmake build)
#
# Exits nonzero if any *hard* check (libraries, apps, BuildConfig, .pc, cmake
# targets) differs.  Header/share differences are reported as INFO.

set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=${WCT_SRC:-$(cd "$HERE/../.." && pwd)}
DEPS=${WCT_DEPS:-}
WORK=${WCT_WORK:-$(mktemp -d)}
GEN=${WCT_GENERATOR:-Unix Makefiles}
CMAKE_ARGS=${WCT_CMAKE_ARGS:-}
WAF_ARGS=${WCT_WAF_ARGS:-}
JOBS=${WCT_JOBS:-$(nproc 2>/dev/null || echo 4)}
RUN_TESTS=${WCT_RUN_TESTS:-0}
WAF_PREFIX=${WCT_WAF_PREFIX:-}
CMAKE_PREFIX=${WCT_CMAKE_PREFIX:-}

mkdir -p "$WORK"
FAILS=0
say()  { printf '%s\n' "$*"; }
hdr()  { printf '\n=== %s ===\n' "$*"; }
pass() { printf '  [PASS] %s\n' "$*"; }
fail() { printf '  [FAIL] %s\n' "$*"; FAILS=$((FAILS+1)); }
info() { printf '  [INFO] %s\n' "$*"; }

# ---------------------------------------------------------------------------
# Build/install each side (unless a prefix was provided).
# ---------------------------------------------------------------------------
build_cmake() {
  [ -n "$CMAKE_PREFIX" ] && { say "cmake: using existing $CMAKE_PREFIX"; return; }
  local b="$WORK/cmake-build" i="$WORK/cmake-install"
  hdr "CMake build+install -> $i"
  local pfx=(); [ -n "$DEPS" ] && pfx=(-DCMAKE_PREFIX_PATH="$DEPS")
  cmake -S "$SRC" -B "$b" -G "$GEN" "${pfx[@]}" \
        -DCMAKE_INSTALL_PREFIX="$i" $CMAKE_ARGS || { fail "cmake configure"; return 1; }
  cmake --build "$b" -j "$JOBS" || { fail "cmake build"; return 1; }
  cmake --install "$b" >/dev/null || { fail "cmake install"; return 1; }
  CMAKE_PREFIX="$i"; CMAKE_BUILD="$b"
}
build_waf() {
  [ -n "$WAF_PREFIX" ] && { say "waf: using existing $WAF_PREFIX"; return; }
  local i="$WORK/waf-install"
  hdr "waf build+install -> $i"
  ( cd "$SRC" && ./wcb configure --prefix="$i" $WAF_ARGS \
      && ./wcb -j "$JOBS" && ./wcb install ) || { fail "waf build/install"; return 1; }
  WAF_PREFIX="$i"
}

# ---------------------------------------------------------------------------
# Comparison helpers.  All operate on install prefixes.
# ---------------------------------------------------------------------------
# WCT-owned library basenames under lib/.
libs_of()  { ( cd "$1/lib" 2>/dev/null && ls libWireCell*.so libWCP*.so 2>/dev/null ) | sort -u; }
# WCT-owned app names under bin/ (WCT apps only, not dependency executables).
apps_of()  { ( cd "$1/bin" 2>/dev/null && ls 2>/dev/null ) \
             | grep -E '^(wire-cell|wcsonnet|wcwires|wcuboone)' | sort -u; }
# WCT-owned public header dirs under include/.
hdrs_of()  { ( cd "$1/include" 2>/dev/null && ls -d WireCell* WCP* 2>/dev/null ) | sort -u; }
# Normalized BuildConfig.h: keep HAVE_*/SPDLOG_ACTIVE_LEVEL, drop version macros.
normbc()   { grep -E '^#define (HAVE_|SPDLOG_ACTIVE)' "$1" 2>/dev/null | grep -vE 'WIRECELL_' | sort; }
pc_reqs()  { grep -E '^Requires:' "$1" 2>/dev/null | sed 's/Requires://' | tr ' ,' '\n\n' | grep . | sort; }
pc_libs()  { grep -E '^Libs:'     "$1" 2>/dev/null | tr ' ' '\n' | grep -E '^-l' | sort; }
# The set of WireCell:: imported targets *defined* by the installed CMake
# package.  waf writes one monolithic config under lib/cmake/WireCellToolkit/;
# the CMake build splits definitions across Config + WireCellTargets*.cmake
# under lib/cmake/WireCell/.  Scan whichever package dir exists (glob matches
# both names) for add_library(WireCell::...) definitions.
cfg_tgts() { grep -rhoE 'add_library\(WireCell::[A-Za-z]+' "$1"/lib/cmake/WireCell*/ 2>/dev/null \
             | sed 's/add_library(//' | sort -u; }

diffcheck() { # name  fileA  fileB  hard|soft
  local name="$1" a="$2" b="$3" sev="$4"
  if diff -u "$a" "$b" >"$WORK/$name.diff" 2>&1; then
    pass "$name"
  else
    if [ "$sev" = hard ]; then fail "$name (see $WORK/$name.diff)"; else info "$name differs (see $WORK/$name.diff)"; fi
    sed 's/^/      /' "$WORK/$name.diff" | head -30
  fi
}

compare() {
  local W="$WAF_PREFIX" C="$CMAKE_PREFIX"
  hdr "Comparing  waf=$W  cmake=$C"

  libs_of  "$W" >"$WORK/libs.waf";  libs_of  "$C" >"$WORK/libs.cmake"
  diffcheck libraries "$WORK/libs.waf" "$WORK/libs.cmake" hard

  apps_of  "$W" >"$WORK/apps.waf";  apps_of  "$C" >"$WORK/apps.cmake"
  diffcheck applications "$WORK/apps.waf" "$WORK/apps.cmake" hard

  hdrs_of  "$W" >"$WORK/hdrs.waf";  hdrs_of  "$C" >"$WORK/hdrs.cmake"
  diffcheck header-packages "$WORK/hdrs.waf" "$WORK/hdrs.cmake" soft

  normbc "$W/include/WireCellUtil/BuildConfig.h" >"$WORK/bc.waf"
  normbc "$C/include/WireCellUtil/BuildConfig.h" >"$WORK/bc.cmake"
  diffcheck BuildConfig.h "$WORK/bc.waf" "$WORK/bc.cmake" hard

  local pcw="$W/lib/pkgconfig/wire-cell-toolkit.pc" pcc="$C/lib/pkgconfig/wire-cell-toolkit.pc"
  pc_libs "$pcw" >"$WORK/pclibs.waf"; pc_libs "$pcc" >"$WORK/pclibs.cmake"
  diffcheck pc-Libs "$WORK/pclibs.waf" "$WORK/pclibs.cmake" hard
  pc_reqs "$pcw" >"$WORK/pcreq.waf"; pc_reqs "$pcc" >"$WORK/pcreq.cmake"
  diffcheck pc-Requires "$WORK/pcreq.waf" "$WORK/pcreq.cmake" hard

  cfg_tgts "$W" >"$WORK/tgts.waf"; cfg_tgts "$C" >"$WORK/tgts.cmake"
  diffcheck cmake-targets "$WORK/tgts.waf" "$WORK/tgts.cmake" hard
}

compare_tests() {
  [ "$RUN_TESTS" = 1 ] || return 0
  hdr "Test outcomes (informational)"
  if [ -n "${CMAKE_BUILD:-}" ]; then
    ( cd "$CMAKE_BUILD" && ctest -j "$JOBS" ) >"$WORK/ctest.log" 2>&1
    info "cmake: $(grep -E '% tests passed' "$WORK/ctest.log" | tail -1)"
  fi
  ( cd "$SRC" && ./wcb --tests -j "$JOBS" ) >"$WORK/waftest.log" 2>&1
  info "waf: $(grep -E 'tests that (pass|fail)' "$WORK/waftest.log" | tr '\n' ' ')"
  info "logs: $WORK/ctest.log  $WORK/waftest.log"
}

# ---------------------------------------------------------------------------
main() {
  say "WCT parity harness"
  say "  src=$SRC  deps=${DEPS:-<none>}  work=$WORK  generator=$GEN"
  build_cmake || true
  build_waf   || true
  [ -n "$WAF_PREFIX" ] && [ -n "$CMAKE_PREFIX" ] || { say "ERROR: need both install prefixes"; exit 2; }
  compare
  compare_tests
  hdr "Summary"
  if [ "$FAILS" -eq 0 ]; then say "  PARITY OK (all hard checks passed)"; exit 0
  else say "  PARITY FAILURES: $FAILS (diffs under $WORK)"; exit 1; fi
}
main "$@"
