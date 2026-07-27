# Merging `origin/master` into `apply-pointcloud` — 2026-07-27

Merge of 479 upstream commits (`origin/master` `66b288a7`, 2026-07-23) into
`apply-pointcloud` (`d1bf5782`), merge base `61900e03` (2026-07-10), with
byte-identity validation of the SBND, PDHD and PDVD reconstruction chains.

- merge commit `5b248b93`
- follow-up `bc8869bc` — `pgraph: drop the LIBTORCH link dependency master added`

**Result: no reconstruction output changes.** Every stage that could be gated
deterministically is byte-identical across the merge. Two gates could not be
decided by identity because the underlying stage is nondeterministic *before*
the merge; both were resolved by same-binary repeat experiments and are
documented below as pre-existing issues, not merge regressions.

## Repro

```bash
cd /nfs/data/1/xqian/toolkit-dev/toolkit

# 0. conflict set, without mutating anything
git merge-tree --write-tree --name-only apply-pointcloud origin/master
#   -> tree ddc6975e..., CONFLICT (content): .gitignore   (that is the whole list)

# the gate scripts live in wcp-porting-img (commit f2b3217):
AB=../wcp-porting-img/abtest

# 1. config gate -- compile every live job against both cfg/ trees
#    (pre-merge this used a cfg/ tree extracted from the merge-tree result:
#     git archive ddc6975e cfg | tar -x -C /home/xqian/tmp/mergecfg/
#     then CFGROOT=/home/xqian/tmp/mergecfg/cfg $AB/compile_all_cfg.sh ...)
$AB/compile_all_cfg.sh /home/xqian/tmp/cfgbase          # before
$AB/compile_all_cfg.sh /home/xqian/tmp/cfgpost          # after
$AB/cmp_cfg.sh /home/xqian/tmp/cfgbase /home/xqian/tmp/cfgpost

# 2. build + unit tests.  NB --notests only skips *running* tests
#    (waft/wcb_unit_test.py:184); the doctest binaries, including master's new
#    all-package `wcdoctest`, are still BUILT by wcbuild, so they are fresh.
./wcb build --notests -p && ./wcb install --notests -p
for p in clus match flash img sigproc aux util iface pgraph gen; do
    ./build/$p/wcdoctest-$p; done

# 3. reco arms (pre-merge first -- baselines expire on rebuild)
$AB/mergegate_run.sh pre        # PDHD+PDVD NF+SP -> img -> clus
$AB/simgate_run.sh pre          # PDHD+PDVD gen/ simulation
$AB/sbndgate_run.sh pre 1 2     # SBND img -> clus -> Q/L
(cd ../wcp-porting-img/sbnd/sbnd_xin && ./run_perf54_nusel.sh pre m66)
(cd ../wcp-porting-img/qlport/scripts && ./sweep_5384.sh m66pre 4)
# ... merge, rebuild, then the same with 'post' / m66post, and compare:
$AB/mergegate_cmp.sh pre post
$AB/fixedinput_gate.sh          # img+clus with SP input held fixed
(cd ../wcp-porting-img/sbnd/sbnd_xin && python3 p54_ab_report.py \
     --base-tag m66pre --opt-tag m66post)
(cd ../wcp-porting-img/qlport/scripts && ./ab_check.sh m66post m66pre)
```

Pre-merge binaries are preserved at `/home/xqian/tmp/premerge-bin/`
(`COMMIT.txt` = `d1bf5782`) so any arm can be re-measured against the old code
without a checkout and rebuild.

## What the merge actually contains

479 commits, dominated by material that cannot reach our chains: a new `spng/`
package (410 files), a complete CMake build (additive — waf is still the
build), and a parallel `cfg/omni/` config tree. Nothing under `cfg/pgrapher`
imports `cfg/omni`, and all our cross-directory jsonnet imports are fully
path-qualified.

The merge is textually trivial: **`.gitignore` was the only conflict**, and all
121 files this branch had changed survive byte-identically. Master modifies
exactly one file inside our packages — `clus/src/SimpleClusGeomHelper.cxx`
(`format(` → `String::format(`, `b29ff85b`) — which auto-merged.

`.gitignore` was resolved as a union: master's 125-line sorted file plus our
`!cfg/pgrapher/experiment/sbnd/sbnd_track_fitting.json` in sorted position
(126 lines), matching how upstream resolved the same conflict in `6ad16499`.

## Gate results

| # | stage | scope | result |
|---|---|---|---|
| 1 | compiled config | 16 jobs, SBND+PDHD+PDVD, reco + sim | **PASS** — identical after `del(._pnode)`; same element counts, component order and Pgrapher `edges` |
| 2 | unit tests | 10 packages | **PASS** — 0 failures (clus 565/565, util 42567 incl. master's new doctests) |
| 3 | simulation (`gen/`) | PDHD + PDVD, track + noise | **PASS** — 8/8 archives identical |
| 4 | SBND PR / taggers | 30 events | **PASS** — 121/121 identical (mabc-pr.zip, pctree, TSV, tracking-stm.root) |
| 5 | SBND img → clus → Q/L | 2 events | **PASS** — 26/26 identical |
| 6 | PDHD+PDVD img + clus | 6 events, SP input held fixed | **PASS** — 178/178 identical |
| 7 | PDHD+PDVD NF+SP | 6 events, 32 archives | first run **29/32** — the 3 differences are pre-existing nondeterminism (§ below); a later full run (`landed`) came out **32/32** |
| 8 | uBooNE | 35 events | Bee zips **35/35 identical**; tagger logs differ, pre-existing nondeterminism (§ below) |
| 9 | full chain, landed binary | 6 events, NF+SP → img → clus | **PASS** — 32/32 + 64/64 + 114/114 identical to the pre-merge arm |

Gate 9 was run last, with the shipped (post-revert) binary, to leave `work/`
coherent (see "State left on disk"). It reproduced the pre-merge output exactly
at *every* stage including NF+SP. That is a stronger end-to-end result than
gate 7, but note what it does and does not show: NF+SP is demonstrably
nondeterministic (§ below), so a 32/32 NF+SP result is one favourable draw, not
proof of determinism. The load-bearing evidence for imaging and clustering
remains gate 6, where the SP input is held fixed by construction.

Gate 6 is the one that matters for imaging and clustering. Because NF+SP is
nondeterministic (gate 7), a naive pre-vs-post comparison of img/clus mixes two
effects. `fixedinput_gate.sh` removes the variable by re-running imaging and
clustering with the post-merge binary on the *exact* SP frames the pre arm
used: **178/178 identical**, i.e. `img/` and `clus/` are unchanged by the merge.

### The two config changes that looked dangerous, and why they are inert

Master rewrites `cfg/pgraph.jsonnet` so `pnode()` injects `_pnode:{nin,nout}`
into every inode, and reimplements `wc.unique_list` (`foldl` → `std.set`+sort),
which feeds `pgraph.edges()` and `pgraph.uses()`. A reordering of `edges[]`
would change Pgrapher's node construction order and therefore output.

Measured, not argued — every live job compiled against both `cfg/` trees with
`wcsonnet` (the production compiler, so this is not subject to the Go-vs-C++
`std.sort` stability question):

- component **order is identical** on all 16 jobs, including the 390-element
  PDVD NF+SP graph, and the `edges` arrays match element for element;
- `_pnode` appears only as a **sibling of `data`** at the config-array element
  root — `jq '[.[]|.data|..|objects|select(has("_pnode"))]|length'` is 0 on
  every job — so it never reaches a component's `configure()`.
  `ConfigManager` reads only `type`/`name`/`data`
  (`util/src/ConfigManager.cxx:41-42,93-95`);
- no new duplicate `type:name` entries appear (the 26/30 duplicates in the
  PDHD/PDVD NF+SP configs are present identically on both arms — pre-existing).

Consequence for future gates: **"compiled config byte-identical" no longer
holds literally.** Compare with `del(._pnode)` (`cmp_cfg.sh` does this).

## `pgraph` libtorch dependency — measured and reverted (`bc8869bc`)

Master's `9cc3c18c` adds `LIBTORCH` to `WireCellPgraph`'s `use` list for
optional NVTX annotations. `--with-nvtx` is off by default so the macros expand
to nothing and pgraph references no torch symbol, but `LINKFLAGS_LIBTORCH`
(`-Wl,--no-as-needed,-ltorch_cuda`) forces a `DT_NEEDED` on `libtorch_cuda`
into `libWireCellPgraph.so`. Every wire-cell job loads that library.

Three-way measurement on the 30-event SBND PR scan:

| arm | `pgraph/wscript_build` | peak RSS | wall |
|---|---|---:|---:|
| `m66pre` | pre-merge | 515 MB | 96 s |
| `m66post` | master's line (`use='WireCellAux LIBTORCH'`) | **812 MB** | 102 s |
| `m66nolt` | reverted (`use='WireCellAux'`) | **516 MB** | 93 s |

`p54_ab_report.py m66post vs m66nolt` = **121/121 identical**, so the
dependency costs ~300 MB per process and ~9% wall while changing no output. The
same ~+300 MB was visible on every PDHD/PDVD NF+SP, imaging and clustering
stage, which matters at the 6-way batch parallelism this box uses. Reverted;
re-add if NVTX profiling is ever wanted.

**Identity chain for the shipped binary.** Gate 4 above compares `m66pre` with
`m66post`, but the binary that ships is `m66nolt` (with the revert). The chain
is `m66pre ≡ m66post` (121/121) and `m66post ≡ m66nolt` (121/121), therefore
`m66pre ≡ m66nolt`. Gates 3, 5 and 6 were additionally re-run against the final
reverted binary and reproduced exactly: 8/8, 26/26, 178/178.

## Pre-existing problems found while validating — NOT caused by the merge

### PDHD NF+SP is nondeterministic run-to-run

`pdhd 027305 evt 0` produced **three distinct content hashes for anode0 in
three identical runs of the PRE-merge binary**, under `setarch x86_64 -R`:

```
run1  132eff40606fa005   run2  0b28333998869f2b   run3  3f3ee6b75d4fe6cc
```

anode2 is bistable (stable across 3 pre-merge runs, but the post-merge binary
produced both `5114f0a8…` and the pre-merge `c9cc3d64…` on two runs). PDVD
`039252 evt 5` anode6 behaves the same way sporadically.

This is why gate 7 cannot be a byte-identity gate today, and it means **any
PDHD/PDVD A/B gate that regenerates SP frames is unreliable** — the doc-54
style harnesses that symlink fixed SP/imaging inputs are not affected. Not
investigated further here; it predates this merge. Related:
`project_wct_run_nondeterminism` (FFTW planner cache) and M4.

### uBooNE `ab_check.sh` gate 2 does not pass against itself

`ab_check.sh` compares `wire-cell-uboone-tagger-compare` logs between two arms.
Running the **same binary** twice (`m66post` vs `m66post2`) gives
`TAGGER: identical=2 diff=33` — the identical signature to the pre-vs-post
comparison. The differing arrays are permutations of the same values, i.e. a
container-ordering effect in the uBooNE tagger output. Gate 1 of that script
(Bee zip content hashes) is deterministic and passed 35/35 in both comparisons.
Treat gate 2's verdict as advisory until the ordering is fixed.

### Stale-install shadowing broke the first build (variant of M1)

The first `./wcb build` failed linking `wcdoctest-hio` with `undefined
reference to WireCell::Hio::show_errors(bool)` — even though the symbol *is*
defined in `hio/src/HIO.cxx:231` and exported by the freshly built
`build/hio/libWireCellHio.so`. Cause: `-L/nfs/.../local/lib` precedes `-Lhio`
on the link line, so the **stale installed** `libWireCellHio.so` (Jul 24,
pre-merge, without the symbol) shadowed the new one. Fixed by
`./wcb install --notests -p --targets=WireCellHio` first, after which the full
build and install succeed. This is the stale-library trap biting at *link*
time; worth knowing whenever master adds a symbol to an existing library.

## Report-only (upstream defects carried in, not fixed here)

- Master's bare `data` ignore line defeats `!sigproc/test/data/*.json` — git
  cannot re-include a file under an excluded directory. Tracked files are
  unaffected; new files under `sigproc/test/data/` would be invisible.
- `iface/src/IChannel.cxx:17-25` — the new `IChannel::global_order()` uses `&`
  where `|` is clearly intended, so it returns 0 for every channel. No callers
  today; do not adopt it.
- `clus/CMakeLists.txt` and `match/CMakeLists.txt` arrive stale against this
  branch's clus sources, and master's CMake `WCT_PACKAGES` omits `flash/`
  entirely. waf is the build, so neither blocks anything.
- `cfg/pgrapher/experiment/pdhd/dnnroi.jsonnet` arrives from master, is
  byte-identical to `cfg/omni/pdhd/dnnroi.jsonnet`, and nothing imports it. Our
  PDHD/PDVD DNN chains import `dnnroi_pp.jsonnet`, whose topology differs —
  a shadowing trap for a future edit.
- The `pgrapher/common/fileio.jsonnet` rot (M12) is unchanged. Note that the
  three `wct-sim-check.jsonnet` files still *compile*: their `io.` uses are
  commented out, and jsonnet imports are lazy, so the missing import is never
  resolved.

## State left on disk

- `<det>/work/<run>_<evt>/` for the six `abtest/events.txt` events was
  regenerated end-to-end (NF+SP → imaging → clustering) with the **landed**
  binary as arm `landed`, so the tree is self-consistent: verified 6/6 events'
  SP frames match the `landed` snapshot by content hash. This matters because
  the validation itself wrote several different SP vintages there — including,
  at one point, frames restored from the *pre-merge* arm by
  `fixedinput_gate.sh` — and NF+SP is nondeterministic, so "just re-run it"
  does not restore a known state (M11).
- Snapshots of every arm are under `/home/xqian/tmp/mergegate/`
  (`pre`, `post`, `post2`, `landed`, `sim-*`, `sbnd-*`, `fixedinput*`), and the
  pre-merge binaries under `/home/xqian/tmp/premerge-bin/` (`COMMIT.txt` =
  `d1bf5782`). These are scratch and will not survive a cleanup; the durable
  record is this doc plus the `m66*` tags
  (`sbnd_xin/work-*-m66{pre,post,nolt}`, `qlport/scripts/sweep/m66*`).
- **`abtest/snap/premerge-20260727/` is not a usable snapshot.** It contains
  only failure logs from an aborted first attempt (imaging skipped because the
  work dirs had no SP frames yet) and was superseded by
  `/home/xqian/tmp/mergegate/`. Do not feed it to `ab_compare.sh`. It was left
  in place rather than deleted because `abtest/snap/` is the protected record
  tree (M13) — the owner should remove it.

## Deliberately not done

`./wcb configure` was **not** re-run. There is no `autoconfig`, so
`./wcb build` reuses `build/c4che/_cache.py`, whose `SUBDIRS` has no `spng` —
the new package is therefore not compiled and the spdlog configure-flag change
is not re-evaluated. Reconfiguring (which forces a full rebuild and compiles
~80 spng sources against libtorch headers) is a separate follow-up.
