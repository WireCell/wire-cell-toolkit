# TaggerCheckTGM: `apply-pointcloud` vs the `tgm` branch

Comparison of two independent ports of the WCP prototype
`WCPPID::ToyFiducial::check_tgm` into WCT's `TaggerCheckTGM` component.

- **A ("ours")** — `apply-pointcloud`, commit `79edf84b`
  (*clus/cfg: TaggerCheckTGM -- through-going-muon tagger port*), 314 lines.
- **B ("colleague's")** — `origin/tgm`
  (<https://github.com/WireCell/wire-cell-toolkit/tree/tgm>), head `449d6f16`;
  the tagger arrives in `dd856c8d` and is amended by `040d32bf`
  (SCE true space) and `b373c0bb` (Bee filter). 325 lines.

Both branches share the merge base `86c39e6c "wcp tgm logic"`, which added the
algorithm walkthrough `clus/docs/tgm/check_tgm.html`. The two ports were
written independently from that walkthrough and from the prototype.

**This document is descriptive only. No code was changed.**

## Repro

```bash
cd /nfs/data/1/xqian/toolkit-dev/toolkit
git fetch origin tgm
git merge-base apply-pointcloud origin/tgm            # 86c39e6c
git log --oneline apply-pointcloud..origin/tgm        # 6 commits
git show origin/tgm:clus/src/TaggerCheckTGM.cxx       # implementation B
sed -n '1331,1602p' prototype_base/pid/src/Cosmic_tagger.h   # the prototype
git diff apply-pointcloud origin/tgm -- \
    cfg/pgrapher/common/clus.jsonnet \
    cfg/pgrapher/experiment/sbnd/clus.jsonnet \
    clus/src/MultiAlgBlobClustering.cxx \
    clus/inc/WireCellClus/MultiAlgBlobClustering.h
```

Prototype reference throughout: `prototype_base/pid/src/Cosmic_tagger.h`,
`check_tgm` at lines 1331–1602.

Note on reading the branch diff: `origin/tgm` forked before a large amount of
unrelated `apply-pointcloud` work (cathode-connect knobs, `steiner
require_beam_flash`, the whole `clus_pr` stage, …). Those show up as
"deletions" in the diff and are branch-age artifacts, **not** TGM differences.
Only the files named above are relevant here.

## Verdict

**Not functionally identical.** Fed the same cluster and the same fiducial
volume, the two will usually agree, but they differ in *which clusters they are
even asked about*, in *how the fiducial margins are applied*, in *whether an
in-beam-window cluster may be tagged*, and in one prototype-transcription index.

What is genuinely shared — both are recognisable ports of the same function:

- pairwise loop over the groups of `cluster.get_extreme_wcps()`;
- per-group "first point that falls outside the FV" scan producing
  `flag_p1_inside` / `p1_index`;
- **CASE A** (both groups exit): 3 interior sample points at ¼, ½, ¾ of the
  exit-to-exit chord; `out_vec_wcps.size()==2` shortcut; the `flag_check_again`
  waypoint re-check through the other extreme groups; the
  `temp_length > 0.45 * length_limit` guard;
- **CASE B** (at least one group looks inside): the
  `|90° − ∠(dir_test, dir_main)| > 75° || (i==0 && k==1)` gate, the 30 cm
  `vhough_transform` direction flipped by −1, the two early-outs
  (`∠(dir, dir_test) > 120°` for p1, `< 60°` for p2), the prolonged-signal
  U/V/W angle test (`<10°`, `<10°`, `<5°`) feeding
  `FiducialUtils::check_signal_processing`, and the
  `|90° − ∠(dir, dir_main)| > 60°` gate feeding
  `FiducialUtils::check_dead_volume`;
- neither ports `check_neutrino_candidate()` (the prototype's Dijkstra
  path-topology neutrino veto) — see the flash-type row in the table below.

## 1. Features in the colleague's (B) that ours (A) does not have

| # | Feature | Where | Notes |
|---|---------|-------|-------|
| 1 | **Debug per-point charge array** | `TaggerCheckTGM.cxx` `mark_debug_points()`; config `debug`, `debug_charge_array` (`"tgm_charge"`), `debug_charge_pcname` (`"tgm_debug"`), `debug_endpoint_charge` (10000), `debug_body_charge` (100) | On each tagged cluster, writes an `npoints`-long `double` pcarray: the two main-axis extreme points get the endpoint charge, everything else the body charge. Default OFF. |
| 2 | **`MultiAlgBlobClustering` `charge_array` / `charge_pcname`** on `bee_points_sets` | `clus/inc/WireCellClus/MultiAlgBlobClustering.h`, `clus/src/MultiAlgBlobClustering.cxx` (`fill_bee_points_from_cluster` gains two args) | When `charge_array` is non-empty, only clusters carrying that cluster-level array are dumped, and per-point Bee charge is read from the array **by global point index** instead of being computed as the per-plane charge mean. Backs feature 1, but is a general-purpose facility. Absent from A entirely. |
| 3 | **Tags *all* live clusters** | `visit()`: `for (auto* cluster : grouping.children())` | A only visits clusters carrying `Flags::main_cluster`. B's choice matches a bulk cosmic-tagging use case; A's matches the prototype, which is called per flash-TPC bundle on that bundle's main cluster (`Cosmic_tagger.h:174`). Practical effect: on the same event B will tag many more clusters. |
| 4 | **SCE-true-space pipeline** | `cfg/pgrapher/experiment/sbnd/clus.jsonnet` (`sce_field` = `SCEFieldTH3/sbnd_dualmap` on CVMFS `SCEoffsets_SBND_E500_dualmap_CV_voxelTH3.root`, wired into `DetectorVolumes` per-APA metadata key `sce_field`; `use_sce` arg on `clus_all_apa`; `switch_scope(correction_name='SCECorrection')` replacing `T0Correction`) | `SCECorrection::forward` = T0 **and** SCE (`clus/inc/WireCellClus/SCECorrection.h`), producing scope `x_sce/y_sce/z_sce`. The whole downstream pipeline — clustering, TGM, Bee — then runs in true space, where the fiducial boundary really is a box. **This is a cfg/pipeline feature, not something inside `TaggerCheckTGM.cxx`**: the component simply box-tests whatever scope points it is handed and is SCE-agnostic. |
| 5 | **One fiducial for everything** | `clus_all_apa` builds `BoxFiducial/all-overall-fv` from `dvm(...).overall` with the `FV_*_margin` values **already subtracted from the bounds**, passes it as `clustering_methods(fiducial=fv_box)`, so `MakeFiducialUtils` — and therefore both the FV test and the dead/SP walks — use the same box. | In A the FV test and the dead/SP walks use *different* fiducials (§2.3). |
| 6 | **Literal prototype wire directions** | `U_dir = (0, cos60, sin60)`, `V_dir = (0, cos60, −sin60)`, `W_dir = (0,1,0)`, hard coded, as in prototype lines 1339–1341 | A computes these from geometry instead (§2.2). B is the more literal transcription; A is the more portable one. |

Also present on `origin/tgm` but outside the scope of this comparison: a truth
labeler with SCE forward map and trackid blob-scalar labels (`fc725327`,
`40aaa672`, `449d6f16`).

## 2. Features in ours (A) that the colleague's (B) does not have

Listed for completeness — this is the opposite direction from the question
asked.

### 2.1 Beam-window protection (largest behavioral divergence)

The prototype gates its tags on `main_flash->get_type()==2` (beam flash): for a
beam flash it requires `!check_neutrino_candidate(...) && temp_length >
0.45*length_limit` before tagging; otherwise it tags unconditionally
(prototype 1416–1437, 1564–1583).

- **A** replaces the flash-type test with a window on `cluster_t0`
  (`beam_window_low`/`beam_window_high`, both 0 by default = gate disabled) and
  is **conservative**: an in-window cluster is never tagged through the
  protected branches, because `check_neutrino_candidate` is unported.
- **B** treats every flash as type ≠ 2 and therefore **always** tags in those
  branches.

Both are documented choices in their respective file headers. They are opposite
choices, and on beam-coincident clusters they give opposite answers.

### 2.2 Geometric per-(apa, face) wire directions

A calls `Facade::compute_wireplane_params(cluster.grouping()->wpids(), m_dv, …)`
and, for each endpoint, resolves the plane via `m_dv->contained_by(p)` to pick
the U/V/W directions of that (apa, face). B hard-codes the uBooNE ±60° /
vertical values.

For SBND the two agree, but not quite trivially. Dumping
`sbnd-wires-larsoft-v1.json.bz2` gives, in `(y, z)` unit components:

| | anode 0 (TPC0) | anode 1 (TPC1) | prototype / B |
|---|---|---|---|
| U | (+0.5, +0.866) = +60° | (+0.5, −0.866) = **−60°** | (+0.5, +0.866) = +60° |
| V | (+0.5, −0.866) = −60° | (+0.5, +0.866) = **+60°** | (+0.5, −0.866) = −60° |
| W | (+1, 0) = 0° | (+1, 0) = 0° | (0, 1, 0) = 0° |

TPC1's U and V are mirrored relative to TPC0, so B's hard-coded pair is
literally wrong there — but harmlessly: the prolonged test is
`(angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)`, and U and V share the same
10° threshold, so swapping them cannot change the result. The test is also
invariant under a sign flip of a wire direction, since `dir_1.angle(wire_dir)`
is unsigned and only `sin(angle)` is used. Net effect on SBND: none. On a
detector with different or asymmetric wire angles, A would follow the geometry
and B would not.

Repro for the table:
```bash
python3 -c "
import bz2,json,math
d=json.load(bz2.open('\$WIRECELL_DATA/sbnd-wires-larsoft-v1.json.bz2'))
st=d['Store']; pts=[p['Point'] for p in st['points']]; ws=[w['Wire'] for w in st['wires']]
for f in st['faces']:
    for pi in f['Face']['planes']:
        p=st['planes'][pi]['Plane']; w=ws[p['wires'][len(p['wires'])//2]]
        t=pts[w['tail']]; h=pts[w['head']]
        dy=h['y']-t['y']; dz=h['z']-t['z']; n=math.hypot(dy,dz)
        print(p['ident'], round(dy/n,4), round(dz/n,4), round(math.degrees(math.atan2(dz,dy)),2))
"
```

A comment in A notes that the `|dir.x|` construction makes the drift-angle test
sign-agnostic, so a single `(1,0,0)` drift axis serves both SBND drift
directions — B relies on the same property implicitly.

### 2.3 Where the fiducial margins enter, and which fiducial each test uses

A's `inside_fv()` is a line-for-line copy of
`FiducialUtils::inside_fiducial_volume` (`clus/src/FiducialUtils.cxx:72–113`) —
the same 6-probe outward-shift walk with the same
`[x_lo,x_hi,y_lo,y_hi,z_lo,z_hi]` / `[x,y,z]` / uniform tolerance-vector
convention — applied to a *separately configured* `IFiducial`
(`NeedFiducial`, config key `fiducial`). B calls
`fid.inside_fiducial_volume(p)` with the default empty tolerance, which reduces
to `m_sd.fiducial->contained(p)`.

So the difference is **not** "A has a tolerance and B does not". It is:

| | margins enter via | FV in/out test uses | `check_signal_processing` / `check_dead_volume` use |
|---|---|---|---|
| **A** | `fv_tolerance` config vector, applied at probe time | `sbnd_pr_fv` — an *unshrunk* `BoxFiducial` spanning both TPCs | the grouping's `FiducialUtils`, built with `fiducial=dv` = `DetectorVolumes`' union of per-face sensitive volumes (which **excludes** the CPA slab, \|x\| < 0.45 cm) |
| **B** | pre-subtracted into the box bounds (`FV_xmin + FV_xmin_margin`, …) | `all-overall-fv` — a *shrunk* `BoxFiducial` from `dvm.overall` | the same `all-overall-fv` box |

A's split is deliberate (its header says the dead-region / SP checks keep the
per-(apa,face) dead-channel logic), but it does mean A's two kinds of FV test
disagree about the CPA slab, whereas B's agree by construction.

### 2.4 Knobs

`length_limit_frac` (default 0.45, the prototype constant) and `enable_case_b`
(default true; when false, CASE B is skipped entirely). B hard-codes 0.45 and
always runs CASE B.

### 2.5 Prototype-faithful `flag_check_again` second sub-loop

Prototype line 1457:

```cpp
Point p3(out_vec_wcps.at(kkk).at(0).x
         + (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(0).x)/4.*(kk+1), …);
```

The mixed endpoints (`at(kkk).at(0)` as the base, but `at(k).at(p2_index)` minus
`at(i).at(0)` as the step) look like a prototype typo, but both ports chose to
keep the shape.

- **A**: `out_vec_wcps[kkk][0] + (pe2 - out_vec_wcps[i][0]) * ((kk+1)/4.)` —
  matches the prototype exactly.
- **B**: `em + (e2 - g0) * (kk+1)/4` with `g0 = out_vec_wcps[0][0]` — uses
  group **0** where the prototype uses group **i**.

They differ only when `i != 0`, i.e. when there are ≥3 extreme-point groups and
the outer pair does not start at group 0. In that case B samples a differently
directed chord and can reach a different `flag_check_again` verdict.

### 2.6 Robustness / integration

- A wraps `check_tgm` in `try`/`catch` and logs a warning on exception;
  B does not.
- A skips clusters already carrying `Flags::TGM` (idempotent re-entry).
- A added a skip in `TaggerCheckSTM` (`clus/src/TaggerCheckSTM.cxx:126–130`) so
  a main already flagged TGM is not re-examined by the STM tagger; B has no STM
  interaction.
- A has a `WCT_TGM_DEBUG` environment-variable debug log printing the endpoint
  pair, the mid-point verdict and the chord/limit lengths. B's debug facility is
  the Bee charge array (§1.1) instead.

## 3. Line-level comparison against the prototype

`P` = `prototype_base/pid/src/Cosmic_tagger.h`.

| P lines | Prototype | A (`apply-pointcloud`) | B (`origin/tgm`) | Same? |
|---|---|---|---|---|
| 1334 | `main_cluster1->get_extreme_wcps()` | same; returns false if <2 groups or any group empty | same | yes |
| 1337–1341 | `drift_dir(1,0,0)`; hard-coded U/V/W | drift `(1,0,0)`; U/V/W from `compute_wireplane_params` per (apa,face) | identical hard-coded values | **no** (§2.2) |
| 1343–1345 | `length_limit` = \|group0[0] − group1[0]\| | same | same | yes |
| 1350–1374 | scan each group for first outside point → `flag_pN_inside`, `pN_index` | same | same | yes |
| 1355 / 1369 | `inside_fiducial_volume(p, offset_x)` | own `IFiducial` + `fv_tolerance` probes; `offset_x = 0` (post-`switch_scope`) | `FiducialUtils::inside_fiducial_volume(p)`, margins in the box; scope already T0(+SCE) corrected | **no** (§2.3) |
| 1406–1411 | CASE A: 3 interior samples on the exit-to-exit chord | same | same | yes |
| 1416–1437 | `flag_check` true: type 2 → `!check_neutrino_candidate && len>0.45·limit`; else → tag | in beam window → **never tag**; else tag | always tag | **no** (§2.1) |
| 1440–1441 | `flag_check` false and exactly 2 groups → tag | same | same | yes |
| 1447–1454 | `flag_check_again` first sub-loop: 4 samples from `p1_index` toward `group[kkk][0]` | same | same | yes |
| 1456–1461 | second sub-loop: base `group[kkk][0]`, step `(group[k][p2_index] − group[i][0])/4` | **same (uses `i`)** | uses `group[0][0]` instead of `group[i][0]` | **no** (§2.5) |
| 1463–1476 | `!flag_check_again`: `!check_neutrino_candidate && len > 0.45·limit` → tag (no flash-type gate here) | in beam window → skip; else `len > length_limit_frac·limit` → tag | `len > 0.45·limit` → tag | equivalent apart from §2.1/§2.4 |
| 1480–1484 | `dir_main` = PCA axis 0; `dir_test` = group[i][0] − group[k][0] | same (`dir_main` hoisted out of the loops) | same (recomputed per pair) | yes |
| 1491 | gate `\|90° − ∠(dir_test,dir_main)\| > 75° \|\| (i==0 && k==1)` | same, expressed as an early `continue` | same, expressed as an `if` block | yes |
| 1496–1499 | p1: 30 cm `VHoughTrans`, ×(−1), `continue` if ∠ > 120° | same (`skip_pair` flag then `continue`) | same (direct `continue`) | yes |
| 1504–1526 | p1 prolonged test → `check_signal_processing(p1, dir, 1 cm)` | same, per-(apa,face) wire dirs | same, hard-coded wire dirs | see §2.2 |
| 1528–1529 | p1 `\|90° − ∠(dir,dir_main)\| > 60°` → `check_dead_volume` | same | same | yes |
| 1533–1538 | p2: Hough, ×(−1), `continue` if ∠ < 60° | same | same | yes |
| 1541–1561 | p2 prolonged + dead-volume, mirroring p1 | same | same | see §2.2 |
| 1564–1583 | both pushed outside: type 2 → neutrino-candidate + 0.45 guard; else tag | in beam window → skip; else tag | tag | **no** (§2.1) |
| 1519 | `angle4` computed and unused ("not added for now", XQ 7/11/2018) | not ported | not ported | yes |
| 1601 | return false | same | same | yes |

Structural note: the prototype's `continue` statements at 1499 and 1538
continue the **`k`** loop, discarding any work already done for `p1`. Both ports
reproduce that — A via a `skip_pair` flag checked after each block, B via a
direct `continue` at the same nesting depth. The two spellings are equivalent.

## 4. Config and pipeline surface

### `cfg/pgrapher/common/clus.jsonnet::tagger_check_tgm`

**A**
```jsonnet
tagger_check_tgm(name="", fiducial="", fv_tolerance=[], beam_window_low=0,
                 beam_window_high=0, length_limit_frac=0.45, enable_case_b=true)
  data: { grouping, fv_tolerance, beam_window_low, beam_window_high,
          length_limit_frac, enable_case_b } + dv + pcts
        + (if fiducial == "" then {} else { fiducial: fiducial })
```

**B**
```jsonnet
tagger_check_tgm(name="", debug=false, debug_charge_array="tgm_charge",
                 debug_charge_pcname="tgm_debug",
                 debug_endpoint_charge=10000.0, debug_body_charge=100.0)
  data: { grouping, debug, debug_charge_array, debug_charge_pcname,
          debug_endpoint_charge, debug_body_charge } + dv + pcts
```

The two key sets are disjoint apart from `grouping`, `detector_volumes` and
`pc_transforms` — a merge would have to union them.

### SBND placement

**A** — inside `clus_pr`, the pattern-recognition job that reloads the persisted
post-Q/L point-cloud tree (`tensor_outname` tarball) through a
`TensorFileSource`. Visitor order:
`switch_scope → steiner → fiducialutils → tagger_check_tgm → tagger_check_stm →
tagger_check_neutrino`, selected by name from `pipeline_names`. The tagger gets
`fiducial=sbnd_pr_fv` (both-TPC box, `x ∈ [−201.05, 201.05] cm`,
`y ∈ [−199.312, 199.312] cm`, `z ∈ [0.85, 500.15] cm`) with
`fv_tolerance = [−2, −2, −2.5, −2.5, −3, −3] cm`, plus the job's beam window.
Scope is `x_t0cor` (+ optional `y_cor/z_cor` per-TPC offsets).

**B** — inside `clus_all_apa`, the in-line all-APA clustering job, appended
after the clustering passes: `… → examine_bundles → fiducialutils →
tagger_check_tgm(debug) → MABC`. `MultiAlgBlobClustering` gains a `"tgm"` entry
in `bee_points_sets` (with `charge_array`, and `filter: 1` after `b373c0bb`, to
clip to the active volume like clustering does). Scope is `x_sce/y_sce/z_sce`
when `use_sce=true`, else `x_t0cor`.

So the two taggers do not even run at the same point in the same job: A runs in
a separate downstream PR job on a reloaded tree; B runs in-line at the tail of
all-APA clustering.

## 5. Smaller observations

- **B's `NeedPCTS` is vestigial.** `TaggerCheckTGM` on `origin/tgm` derives from
  `Clus::NeedPCTS` and calls `NeedPCTS::configure(config)`, but nothing in the
  file uses the transform set — the SCE application moved out of the component
  in `040d32bf`. The `pc_transforms` config key is now inert for this component.
- **B has a stale header comment.** The `.cxx` header says the pipeline inserts
  `switch_scope(SCECorrection, backward, t0=0)` materializing
  `x_sce = backward(x_t0cor)`. The cfg comment and `SCECorrection.h` both
  describe the *forward* direction from raw: `x_t0 = x_raw − dirx·t0·v_drift`,
  then `x_sce = x_t0 + displacement`. The cfg and header agree with each other;
  only the `.cxx` comment is out of date. Documentation only — no behavior
  claim rests on it.
- **Loggers differ**: A uses `clus.NeutrinoPattern` (shared with the other
  taggers) at INFO per cluster; B uses a dedicated `clus.TGM` at DEBUG with a
  single per-event summary.
- **Neither ports `check_neutrino_candidate`.** It remains the one substantive
  piece of the prototype function that is missing from both.

## 6. Open questions (for a later decision — no recommendation here)

1. **Which clusters should be tagged** — every live cluster (B) or only
   `main_cluster`-flagged ones (A)? These serve different downstream consumers
   (bulk cosmic rejection vs the per-bundle neutrino-selection chain).
2. **Should the beam-window protection survive a merge?** A refuses to tag
   in-window clusters; B tags them. Resolving this probably means deciding
   whether `check_neutrino_candidate` gets ported.
3. **Should the SCE true-space scope be adopted more widely?** It is
   independent of the tagger, and it changes the meaning of the box FV for
   every downstream consumer, not just TGM.
4. **Should the MABC `charge_array` / `charge_pcname` Bee facility be adopted?**
   It is generally useful for visual debugging of any per-cluster verdict and
   is not coupled to TGM.
5. **`flag_check_again` index** (§2.5): confirm whether the prototype's `at(i)`
   is intentional or a typo before either implementation is changed —
   `clus/docs/porting/porting_dictionary.md` does not currently list it (M15).
6. **Fiducial consistency** (§2.3): should the FV test and the dead/SP walks
   share one fiducial (B) or stay split (A)?
