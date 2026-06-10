# PDHD: dead blobs missing from Bee links — diagnosis & fix

**Status:** FIXED & verified (run027409 evts 0–4). **Date:** 2026-06-07.
**Reference good example:** SBND (`sbnd_xin`) chain, which renders the dead layer correctly.

## TL;DR of the fix (four changes, all landed)

1. `gen/src/Reframer.cxx` + `Reframer.h` — new opt-in `keep_masks` (default
   false → byte-identical for every existing config). The PDHD imaging Reframer
   was silently dropping the channel-mask-map, so the dead/masked fork saw zero
   bad channels. `pdhd/wct-img-all.jsonnet` sets `keep_masks: true`.
2. `pdhd/wct-img-all.jsonnet` — override `det.volumes` to the **data**
   geometry (one APA face `null`/non-sensitive) while keeping simparams' 500 ns
   tick. simparams alone marks BOTH faces sensitive, which made the dead fork
   tile phantom blobs on the driftless face.
3. `img/src/GridTiling.cxx` — emit no blobs on a non-sensitive face
   (`m_face->sensitive().empty()`). Byte-identical wherever faces are sensitive
   (every existing config); suppresses the phantom driftless-face dead blobs.
4. `clus/src/clustering_deghost.cxx` — null-guard in `find_max_cluster`
   (defensive; the `get_closest_2d_point_info` empty-index sentinel `-1.0`
   would otherwise be counted as a "close" null cluster and `at()`-throw).
5. `pdhd/clus.jsonnet` — `dead_area_version: 2` on the MABC nodes so the dead
   slab lands on the correct PD anode face (tpc=apa, not -1).

Verified: dead area now appears only on each APA's sensitive face
(apa0→face0, apa1→face1, apa2→face0, apa3→face1); clustering emits 0 empty 2D
indices and runs clean; Bee link carries live `clustering-global` + dead
`channel-deadarea` for all 5 events.

## Symptom

The PDHD Bee link (run027409, 5 events) shows live charge blobs but **no dead
("channel-deadarea") layer**. SBND Bee uploads show it.

## Root cause (verified `file:line` + runtime logs)

The dead layer is **empty at the imaging source**, for two independent reasons.

### 1. The Bee link was built by the imaging-only path, which has no dead capability

`pdhd/run_bee_img_evt.sh` converts clusters with `wirecell-img bee-blobs`
(`wire-cell-python wirecell/img/__main__.py:474`). That tool only ever writes
`type="cluster"` from blob points — it has **no dead-area code path at all**.
So a link produced this way can never contain a dead layer, by construction.

The dead layer is instead produced by the **clustering** node
`MultiAlgBlobClustering` (`clus/src/MultiAlgBlobClustering.cxx`):
`save_deadarea` → `fill_bee_patches_from_grouping("dead")` (:1961) →
`channel-deadarea-apa%d-face%d.json` (:309). PDHD `pdhd/clus.jsonnet:223`
already sets `save_deadarea: true`.

### 2. PDHD imaging produces an EMPTY dead grouping → MABC writes no dead area

`MultiAlgBlobClustering` only flushes a patch if `patches.size() > 0`
(`clus/src/MultiAlgBlobClustering.cxx:449-455`). PDHD's are empty because the
masked/dead imaging fork emits no blobs:

- All `work/<run>_<evt>/clusters-apa-apa*-ms-masked.tar.gz` are **155 bytes**
  (one trivial graph), across every event and anode. Active tarballs are MB-sized.
- Imaging log (`work/027409_0/wct_img_027409_0.log`):
  `<MaskSlice:slicing-apa0-ms-masked_0> ... cmm:[ bad:0 ]` and `cmm.size(): 0`.

  The masked 2-view fork (`pdhd/img.jsonnet:173`,
  `masked_planes=[[0,1],[1,2],[0,2]]`) needs the bad/dead channel map to know
  where to tile. It receives **zero** bad channels.

- The SP archive **does** contain a non-empty mask: `chanmask_bad_40896.npy`
  has shape `(34, 3)` (34 fully-dead channel ranges). It is simply **not being
  loaded into the frame's channel-mask-map** by `FrameFileSource` in
  `pdhd/wct-img-all.jsonnet`.

- **Not DNN-ROI specific.** Imaging the *standard* SP frames
  (`input_data/.../protodunehd-sp-frames`) gives the identical `cmm bad:0` and a
  155-byte masked tarball. So this is a general PDHD imaging gap, independent of
  the DNN-ROI/L1SP chain.

### SBND contrast (why it works)

Same imaging skeleton, but SBND's imaging log shows
`cmm:[ bad:92 ]` / `cmm.size(): 92` reaching the masked slicer, its
`icluster-apa*-masked.npz` are ~900 KB (real dead blobs), and its
`work/ql_*/mabc-*.zip` contain `channel-deadarea-apa{0,1}-face0.json`. SBND
frames also carry `chanmask_bad` AND the CMM actually populates (`bad:92`).
The SBND frame additionally shows a frame tag (`tags:[ "orig1" ]`) where the
PDHD frame shows `tags:[ ]` — a candidate lead for why the mask associates on
SBND but not PDHD.

### 3. Secondary, independent gap: `dead_area_version` not set in PDHD

SBND `cfg/.../sbnd/clus.jsonnet:241,350` sets `dead_area_version: 2`
(wire-cell-bee3 v2 wrapper with `tpc=apa`). PDHD `clus.jsonnet` omits it →
defaults to v1 with `tpc=-1` (`MultiAlgBlobClustering.cxx:313`), so even once
dead blobs exist, the dead slab would not land on the correct PD anode face.
Required for PD multi-face geometry, but **not sufficient** while the dead
grouping is empty.

## Open question (needs decision / PD-domain input)

Why does PDHD's `FrameFileSource` deliver `cmm bad:0` while SBND's delivers
`bad:92`, when both archives carry `chanmask_bad`? Candidates:
1. PDHD SP frame sink writes `chanmask_bad` but not in a form `FrameFileSource`
   re-associates (frame tag / ident keying — note PDHD `tags:[ ]` vs SBND
   `tags:[ "orig1" ]`).
2. `FrameFileSource` in `pdhd/wct-img-all.jsonnet` needs explicit config to load
   the channel-mask-map.
3. PDHD masked-fork output is partly expected to be sparse (34 scattered
   single-plane-dead ranges may yield few 2-view-dead crossings), but `bad:0`
   means the fork never even sees the dead channels — so this is a real gap,
   not just sparsity.

## UPDATE — the crash, and the deeper geometry bug it exposed

Enabling the dead blobs (Reframer `keep_masks`) made `run_clus_evt.sh` crash with
`std::out_of_range: unordered_map::at` in `clustering_deghost.cxx:217`
(`ClusteringDeghost::visit` → per-APA deghost). Backtrace key frame:
`find_max_cluster` lambda doing `map_cluster_index.at(nullptr)`.

### Immediate crash mechanism
`DynamicPointCloud::get_closest_2d_point_info` returns `{-1.0, nullptr, -1}`
when the per-(plane,face,apa) KD index is empty (`DynamicPointCloud.cxx:310`).
In `clustering_deghost.cxx`, the caller gates on `distance <= dis_cut/3`, and
**-1.0 passes that test**, so the `nullptr` cluster gets counted into
`map_cluster_num`; `find_max_cluster` then `at(nullptr)`-throws.
Guarded (`if (!c) continue;`, `clustering_deghost.cxx:216`) — byte-identical for
any input without a null key (everything that worked before).

### Root geometry bug (the user's "something not wired right")
The empty index that produced the null is **`a1f0pW` — apa1, face0, W plane —
which is the NON-SENSITIVE face**. Each PDHD APA has exactly one non-sensitive
face (`AnodePlane.cxx:159-163`, "anode N face M is not sensitive": no cathode →
empty sensitive volume). Confirmed per event:

| APA | non-sensitive face | live index | dead area emitted on |
|---|---|---|---|
| apa0 | face1 | a0f1 empty | face0 only ✓ |
| apa1 | face0 | a1f0 empty | face0 **+** face1 ✗ |
| apa2 | face1 | a2f1 empty | face0 **+** face1 ✗ |
| apa3 | face0 | a3f0 empty | face0 **+** face1 ✗ |

The **active** imaging correctly emits nothing on the non-sensitive face (no
charge → empty live index). The **masked/dead** fork tiles dead channels purely
geometrically, so it emits dead blobs on BOTH faces — including the
non-sensitive one, which has no drift volume. Those phantom dead blobs are
(a) wrong in the Bee dead layer and (b) exactly what makes the deghost project a
live point onto the empty non-sensitive-face index → null → crash.

### Fix applied (chosen gate: GridTiling sensitivity)
Suppress dead/masked blobs on non-sensitive faces at `GridTiling`: it now emits
nothing when `m_face->sensitive()` is empty/degenerate
(`AnodeFace::sensitive()`, `gen/inc/WireCellGen/AnodeFace.h:32`). This is
byte-identical wherever the tiled faces are sensitive (every existing config);
it changes output only where it fixes this bug.

Crucially, the gate is inert unless imaging *knows* the face is non-sensitive.
`pdhd/wct-img-all.jsonnet` previously used `simparams.jsonnet`, which marks BOTH
APA faces sensitive ("cryostat side included"). The fix keeps simparams (for its
500 ns tick — data is resampled 512→500 before imaging) but overrides
`det.volumes` with the **data** geometry (`params.jsonnet`, one face `null`), so
each APA's driftless face is correctly non-sensitive and GridTiling skips it.

The `clustering_deghost` null-guard stays regardless (defensive; the sentinel
`-1.0`-as-"close" is a latent trap independent of PD). Note it is crash-safety
only — with the guard the phantom point would still be mis-counted into the
deghost thresholds, so removing the phantom blobs at the imaging source (above)
is the real fix, not the guard.

## Why this is "wired right" now (answers "shouldn't the code skip the back face?")

Live/active imaging *does* skip the driftless face naturally: no charge → no
activity → no live blobs (the live index is empty, always). The dead/masked fork
does **not** skip naturally — it tiles dead *channels* geometrically, independent
of charge, so it needs to be *told* a face is driftless. The fix tells it
(data `det.volumes`) and makes the tiler honor it (GridTiling sensitivity check).
