# PDVD: dead-blob placement bugs in the Bee display — diagnosis & fix

**Status:** FIXED & verified (run039324 evts 0–4). **Date:** 2026-06-07.

Three independent bugs surfaced by the PDVD run039324 dead-area display. Two are
in the **wire-cell-bee3 viewer** (geometry table); one is in the **toolkit**
(dead-blob corner reconstruction). The live/active clustering was always correct.

## Symptoms

1. **Anodes 1, 3, 4, 6**: dead-area slab drawn on the **cathode side** instead of
   the anode plane.
2. **Anode 2**: dead blobs spill ~4 cm **past the anode's Y boundary** (the
   thing visible in the user's screenshot).

## Root causes & fixes

### Bug 1 — viewer `driftDir` (symptom 1). `wire-cell-bee3`
`experiment.js` base `driftDir(i) = ((i%2)-0.5)*-2` alternates every index, but
PDVD is grouped: anodes 0–3 bottom-drift, 4–7 top-drift. The dead slab is placed
at `anode_x = center − driftDir·halfx`, so the alternating sign is wrong for
exactly **{1,3,4,6}**. Fix: override `driftDir` in `ProtoDUNEVD` → `+1` for 0–3,
`−1` for 4–7. (Committed to `wire-cell-bee3` `main`, `7d4b5b7`.)

### Bug 2 — viewer TPC YZ box mapping (the empty/mis-placed box). `wire-cell-bee3`
The `ProtoDUNEVD` TPC table had rows **1↔2 and 5↔6 swapped** relative to the
LArSoft anode numbering (verified against `protodunevd-wires-larsoft-v3`: anode 2
is Y[0.6,336]/Z[0,149], but row 2 was Y[-337,0]/Z[150,299]). The per-anode box
outline is drawn from this table, so anode 2's box landed in the wrong quadrant.
Fix: swap the two row pairs so box `i` matches anode `i`. (Same bee3 commit.)
**Note:** the box YZ table does NOT move the dead blobs themselves (they render at
their true world Y,Z minus a global union offset) — it only draws the outline.

### Bug 3 — toolkit dead-blob corner reconstruction (symptom 2 spillover). toolkit
The dead blob itself is built correctly at imaging time (raw blob corners are
clipped to the wire extent, Y≥0.6; the live 3-view blobs never overflow). The
spillover is injected on **reload**:

- The dead-area bee patch (`MultiAlgBlobClustering::fill_bee_patches_from_cluster`)
  uses the per-blob `"corner"` point cloud, built by
  `Aux::make_corner_dataset` from `iblob->shape().corners()`.
- For JSON cluster files (PDVD), `ClusterLoader::to_blob`
  (`aux/src/ClusterHelpersLoader.cxx`) rebuilds the RayGrid `Blob` from the wire
  strips via `add()` **without reconstructing the two boundary layers** (see the
  `// fixme: how to set bounding box?`). For a **2-view** (dead) blob — 2 real
  planes + 1 full-width dummy plane — the missing boundary layers let the
  re-derived corners run to the dummy plane's physical edge, ~4 cm past the wires
  on one edge of anode 2 (deterministic; ~0 elsewhere because their dead channels
  don't reach an edge in the unconstrained direction).
- The **correct** corners are already saved in the JSON (`jblob["corners"]`,
  cartesian x,y,z) but were discarded on load.

**Fix (toggle, default OFF → byte-identical for every existing config):**
- `ClusterFileSource` gains `restore_corners` (default false). When true, the
  JSON loader carries the saved corners onto the `SimpleBlob`
  (`set_stored_corners`).
- `make_corner_dataset` uses the stored corners (via `dynamic_cast` to
  `Aux::SimpleBlob`) when present, else falls back to the reloaded shape.
- This touches **only** the dead-area `"corner"` PC (its sole consumer is the
  dead bee patch; the live consumer in `Facade_Blob` is commented out). The
  RayGrid `shape().corners()` is intentionally left as-is — fixing the
  boundary-layer reconstruction would change geometry for all consumers and is
  out of scope.
- Enabled in `pdvd/wct-clustering.jsonnet` (`restore_corners: true`).
- The numpy cluster-file path (`ClusterArrays::to_cluster`) is **also covered**:
  the same `restore_corners` flag carries the stored (y,z) corners from the `b`
  array (count at col `(sigu_col+1)+4+6 = 14`, then y,z pairs) onto the blob
  (x set to 0 — the dead patch consumes only y,z).

## Cross-detector status (does each chain have the bug, and is it fixed?)
The overflow only occurs for **2-view dead blobs** (2 real planes + 1 full-width
dummy) when the loader drops the RayGrid boundary layers.  Whether a chain is
affected depends on its load path **and** whether its dead blobs are 2-view.

| chain | load path | dead-blob views (tested evt) | overflow? | `restore_corners` |
|---|---|---|---|---|
| **PDVD** run039324 | JSON (`ClusterLoader`) | 2-view at anode-2 edge | **YES, −3.9 cm** | enabled, **fixes it** |
| **PDHD** 027380 evt0 | JSON (`ClusterLoader`) | all 680 blobs **3-view** (max plane width 3–8 wires) | no (10 µm rounding only) | enabled, protective + byte-identical |
| **SBND** evt11 (`wct-clustering`/`perevt`) | numpy (`ClusterArrays`) | 3-view; numpy `{0,1}`+nudge reconstruction already faithful | no | enabled, **faithful no-op** (sentinel-proven the readback fires) |

Notes:
- PDHD/SBND don't exhibit PDVD's spill in the tested events because their dead
  blobs are 3-view (fully constrained — boundary layers redundant).  The flag is
  enabled on the same proven code path so a 2-view edge blob, if it ever appears,
  is corrected; live output is byte-identical in all cases.
- **SBND `wct-clus-matching-standalone.jsonnet` (the main `run_clust_QL_evt.sh`
  BEE chain) images the dead view in-graph** (not file-loaded), so its dead blobs
  use fresh imaging-time corners and never had the bug — `restore_corners` is not
  set there (its `ClusterFileSource` loads only active).
- The numpy readback was verified live by a sentinel: injecting `y+12345` shifted
  the SBND dead corners by exactly +1234.5 cm, then reverted.
- Configs enabling it: `pdvd/wct-clustering.jsonnet`, `pdhd/wct-clustering.jsonnet`,
  `sbnd_xin/wct-clustering.jsonnet`, `sbnd_xin/wct-clus-matching-perevt.jsonnet`.

Files: `aux/inc/WireCellAux/SimpleBlob.h`, `aux/inc/WireCellAux/ClusterHelpers.h`,
`aux/src/ClusterHelpersLoader.cxx` (JSON), `aux/inc/WireCellAux/ClusterArrays.h`
+ `aux/src/ClusterArrays.cxx` (numpy), `aux/src/SamplingHelpers.cxx`,
`sio/inc/WireCellSio/ClusterFileSource.h`, `sio/src/ClusterFileSource.cxx`,
and the `*/wct-clustering.jsonnet` (+ SBND `wct-clus-matching-perevt.jsonnet`)
configs that set `restore_corners: true`.

## Verification (run039324, 5 events)
- Anode 2 dead-area Y_min: **−3.3 → 0.60 cm** (exactly the wire edge); all 8
  anodes show 0 overflow.
- Dead patches non-empty, polygon counts unchanged (≈40–100/anode).
- Live/active layers (`clustering-global`, `clustering-apaN-faceM`) **byte-identical**
  old-vs-new with the flag ON — the restored corners are inert for live blobs.

## Which fix resolves which symptom (coordination)
The toolkit fix is **data-side**; the two bee3 fixes are **viewer-side** (already
on bee3 `main`). For a regenerated mabc/bee zip:
- new zip + **old** viewer → spillover fixed, drift/box still wrong;
- new zip + **bee3 `main`** viewer → all three fixed.
