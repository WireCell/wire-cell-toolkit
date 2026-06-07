# Audit: toolkit wire geometry & channel→wire mapping vs official ProtoDUNE (HD/VD)

**Status:** investigation only — *no code/config/wire-file changes made*.
**Date:** 2026-06-07.

> **Update (same day):** the initial draft of this doc concluded, from git
> provenance alone, that the PDHD wire file was "strongly-evidenced stale." That
> was **revised after empirically regenerating the wire files from the official
> GDMLs** (dunecore/dunereco cloned at the matching release tag). The regeneration
> shows **PDHD is current and reproducible**, and surfaces a **concrete PDVD U/V
> channel-assignment discrepancy** instead. The provenance argument was misleading;
> the empirical results below supersede it.

## 1. Scope, method, and what was actually done

Compares the WireCell **toolkit** wire-geometry files for ProtoDUNE-HD (PDHD) and
ProtoDUNE-VD (PDVD) against the **official** DUNE geometry / channel maps, with
focus on the **PDHD "3-wire shift"** (known issue, reportedly fixed upstream) and
a suspected **stale/incorrect PDVD** mapping.

**Repos used (cloned into `../` at tag `v10_20_08d00`, matching `duneprototypes`):**

| Repo | Path | Provides |
|---|---|---|
| `duneprototypes` (`protodunecode` symlink) | `../duneprototypes` | official channel-map text + mapmakers |
| `dunecore` | `../dunecore` | official GDML geometry (`…/Geometry/gdml/`) |
| `dunereco` | `../dunereco` | official wire-cell reference wires (`…/DUNEWireCell/`) |
| generator | `../wire-cell-python/wirecell/util/wires/gdml.py` | `wirecell-util gdml-to-wires` |

**Method (read-only):** decompressed every wire JSON and tallied
anodes/planes/channels/positions; read the official channel-map text + mapmakers;
**regenerated** the wire files from the current GDMLs with `wirecell-util
gdml-to-wires` and diffed channel↔position per face-plane; ran the upstream
integration tests. No wire file, config, or code was modified.

---

## 2. Disambiguation: there are TWO different "±3 shifts"

Easy to conflate; they live on **different axes**.

### (A) Official "visible-wire" ±3 channel shift — channel ↔ wire-position
A *cyclic permutation of offline channels within each U/V plane* so the offline
channel numbering follows where wires **emerge from under the head boards**
(the visible/active portion, matching the offline geometry's wire endpoints)
rather than their **soldered anchor points**.

- Intro in offline channel map **v3** — `duneprototypes` `6bbccf25` (2022-06-29):
  *"V3 PDHD channel map with +-3 channel shifts in U and V to account for
  anchor-point vs active area shifts."*
- Carried to WIBEth as `…_visiblewires_v1.txt` (`76aadccd`, 2024-08-05).
- **Sign-flipped 2025-06-30** — `c8f43809` / PR #89: *"flip the sign of the +/- 3
  wire shift … fix the cathode-crossing inefficiency and a sign mistake."*
- `+3` on inverted (North/Lower) APAs, `-3` on upright (South/Upper); **U/V only,
  W exempt**. Verified directly in the text maps: **v5 = (v4 + 3) mod 800** within
  each U/V plane block (e.g. v4 planechan-0 offlchans 119,159,199,… → v5
  122,162,202,…). This lives **entirely in the offline-channel labeling
  (LArSoft channel map)** — it does *not* move wire positions.

### (B) Toolkit `coh_group_shift` — channel ↔ FEMB electronics grouping
`cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet` (`9e151a19`, 2026-04-28): pure
offline-channel arithmetic `n*2560 + std.mod(40*u + shift + j, 800)` deciding
which 40 offline channels share a coherent-noise group; `+3` on anodes 0&2 /
`-3` on 1&3, U/V only; from a run-027409 FEMB-edge correlation audit.

Both are magnitude 3, U/V-only, APA-orientation sign-flipped — but (A) is
channel↔**wire-position** and (B) is channel↔**FEMB grouping**.

---

## 3. PDHD — wire file is CURRENT and reproducible (no geometry-file shift)

**Config:** `cfg/pgrapher/experiment/pdhd/params.jsonnet:167`
`wires: "protodunehd-wires-larsoft-v1.json.bz2"`.

**Finding 1 — identical to the official release.** The toolkit file is
**byte-identical** (md5 `ef190f4b…`, and identical decompressed) to the official
`dunereco/dunereco/DUNEWireCell/pdhd/protodunehd-wires-larsoft-v1.json.bz2` at the
current release tag **`v10_20_08d00`**. The 2022 git date in `wire-cell-data` is
misleading: the *content* equals the current DUNE release. **It is not stale vs
upstream.**

**Finding 2 — exactly reproducible from the current GDML.** Regenerating with
```bash
wirecell-util gdml-to-wires -d protodunehd_v8 -o regen.json.bz2 \
   ../dunecore/dunecore/Geometry/gdml/protodunehd_v8_refactored.gdml
```
yields the **same channel↔position mapping in all 24 face-planes** — `0/8`
differing in U, `0/8` in V, `0/8` in W (seg-0 wires, ordered by pitch). Structure
identical (4 anodes / 8 faces / 24 planes / 22 208 wires / 44 416 points).

**Finding 3 — the channel numbering is purely geometric.** `build_hd_channel_map`
(`gdml.py`) allocates channels **sequentially by pitch position** over the GDML's
**active-wire** endpoints (face1 then face0, with the documented U/V
descending/ascending pitch convention), explicitly "to reproduce the reference
file convention." The ±3 *visible-wire* shift (§2A) is **not** encoded in this
geometric file; it is a LArSoft channel-**map** convention applied at decode time.

**Conclusion (PDHD).** The wire-geometry file is **current and self-consistent**:
content = current DUNE release, and reproducible from the current `protodunehd_v8`
GDML. Because its channels are pitch-ordered from the **active-wire** GDML
endpoints — the same "emergence" geometry the **visible-wire** convention targets,
and the channel map was **sign-flipped in 2025-06-30 to match that geometry** —
the file is consistent with the corrected (post-2025) convention. **There is no
stale "3-wire shift" baked into the wire-geometry file.**

> **One residual, honestly bounded gap.** This proves the file is consistent with
> the GDML geometry, *not* that a given data file's per-trace channel labels match
> it. The ±3 visible-wire shift lives in LArSoft's channel map; confirming the
> *data path* is byte-aligned needs a LArSoft channel-map dump
> (`PD2HDChannelMapService` at runtime), which requires the LArSoft runtime (not
> available here). All on-disk evidence points to consistency.

The toolkit `coh_group_shift` **(B)** therefore remains a **separate**, empirically
real FEMB-electronics grouping correction for the specific run-027409 data — it is
**not** evidence of a wire-geometry-file bug. (If anything, (B) is the place a
residual data-vs-convention ±3 would show up, and it was patched in config.)

---

## 4. PDVD — geometry is right, but the U/V channel↔wire assignment is UNVERIFIED and differs from the current converter

**Config:** `cfg/pgrapher/experiment/protodunevd/params.jsonnet:180`
`wires: "protodunevd-wires-larsoft-v3.json.bz2"`.

**Finding 1 — identical to the official release, and channel-count-correct.**
Byte-identical (md5 `217c2f66…`) to
`dunereco/.../protodunevd/protodunevd-wires-larsoft-v3.json.bz2` at `v10_20_08d00`.
v3 has **8 anodes (idents 0–7)** and **12 288 channels (0–12 287)** — matching the
official scheme (`vd/ChannelMap/mapmakers/addcrpnum.cxx`). The old v1 (16 anodes /
12 304 ch) is a superseded numbering. **v3 is the latest, count-correct file.**

**Finding 2 — geometry (positions) matches the current v4 *and* v5 GDML.** The v4
and v5 GDML produce **identical** channel↔position mappings (`0/48` face-planes
differ between `regen-v4` and `regen-v5`); v5 only adds more wire *segments*
(25 344 vs 23 792 points) without changing the seg-0 channel↔position. The upstream
`test_gdml_integration_v4` passes (19/19), validating Z-plane wire **counts and
endpoint coordinates** (0.1 mm tol) against the shipped v3. **So the v4-vs-v5
"which geometry" question is moot for the channel↔wire map** — they agree.

**Finding 3 — the shipped v3's U/V channel↔wire assignment does NOT match the
current converter, and upstream does not test it.** Diffing the shipped v3 against
`regen-v4`/`regen-v5` (same converter that reproduces PDHD perfectly):

| plane | face-planes differing (of 16) |
|---|---|
| U | 12 |
| V | **16** |
| W | 12 |

The wire **positions and segment histograms are identical**; what differs is the
**global channel number assigned to the same physical wire** (e.g. anode-0 face-0:
V starts at offline 1239 in shipped vs 1055 in regen; W at 2196 vs 1244). The
upstream `test_gdml_integration_v4` **passes anyway because it deliberately does
not check U/V channel assignment** — its own comment: *"The U/V planes of the
positive-X (seg=1) face can swap relative to the reference due to different
drift-order conventions, so we only assert the Z-plane count and the total per
face."* It checks only **Z-plane geometry**, never the U/V channel↔wire pairing.

**Conclusion (PDVD).** The v3 **geometry** is correct and current, but its **U/V
channel-to-wire assignment is on a convention that the current wire-cell converter
does not reproduce and that no upstream test validates.** This is the concrete
locus matching the user's suspicion. It is *not* a ±3-style shift and *not* a
position error — it is a **U/V drift-order / channel-block convention** difference
between the (older) tool that produced the shipped v3 and the current converter.
**Which one matches the LArSoft VD offline channels (`PD2VDTPCChannelMap_v2`) is
unconfirmed** and is the thing to settle.

---

## 5. Why HD reproduces but VD doesn't (and the cross-cutting NF hypothesis)

- **HD** channels reproduce exactly because `build_hd_channel_map` is the same tool
  (or direct ancestor) that made the HD reference, with a fully pinned pitch
  convention. **VD**'s `assign_vd_channels` has an acknowledged U/V drift-order
  ambiguity, and the shipped v3 was evidently produced by a different/older
  channel-assignment pass — hence the mismatch.
- **Cross-cutting hypothesis (flagged, not asserted):** the 2026-04
  `coh_group_shift` (B) is empirical evidence that the run-027409 *data* needed a
  ±3 offline-channel correction in U/V. Given §3, the HD **wire geometry** is not
  the cause; the more likely story is a **data/channel-map-version** mismatch
  (which offline-channel convention the data was decoded with) rather than a
  geometry-file bug. Confirming this needs a LArSoft channel-map dump.

---

## 6. Recommended next steps (no change made here)

1. **PDHD — treat as good.** No wire-file change indicated. If a residual data-vs-
   convention ±3 is suspected, dump `PD2HDChannelMapService` for run-027409's
   release and compare its offlchan↔(plane,wire) to the wire file's
   channel↔position — that, not the geometry file, is where (A) would bite. Keep
   `coh_group_shift` (B) as the electronics-grouping fix it is.
2. **PDVD — resolve the U/V channel convention.** Compare the shipped v3 and the
   converter's `regen` U/V channel↔wire against the LArSoft VD channel map
   (`vd/ChannelMap/PD2VDTPCChannelMap_v2.txt`, 2025-07) to decide which is correct.
   If the converter is correct, regenerate and adopt
   `wirecell-util gdml-to-wires -d protodunevd_v5 …` (v5 = current geometry); if the
   shipped v3 is correct, the converter's `assign_vd_channels` U/V drift-order needs
   fixing and the upstream test should be extended to assert U/V channel↔wire.
3. Any wire-file swap is a **geometry change** — gate behind the usual
   bit-identical-by-default discipline and rebaseline imaging/clustering snapshots.

---

## 7. Empirical results table (this audit)

| check | PDHD | PDVD |
|---|---|---|
| toolkit file == dunereco release `v10_20_08d00`? | **yes (byte-identical)** | **yes (byte-identical)** |
| reproducible from current GDML (channel↔position)? | **yes, 0/24 face-planes differ** | **no — U 12/16, V 16/16, W 12/16 differ** |
| geometry (positions) matches current GDML? | yes | yes (Z-plane, 0.1 mm, test passes) |
| upstream test covers U/V channel↔wire? | (HD test) — n/a here | **no (Z-plane only, by design)** |
| channel numbering source | pure geometric pitch order (active-wire GDML) | differs between shipped v3 and converter |
| verdict | **current & consistent; no geometry-file shift** | **geometry OK; U/V channel↔wire convention unverified / differs** |

---

## 8. Reference appendix

**Toolkit** — `cfg/pgrapher/experiment/pdhd/params.jsonnet:167` (PDHD wires);
`cfg/pgrapher/experiment/protodunevd/params.jsonnet:180` (PDVD wires);
`cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet` (`coh_group_shift`, B);
`wire-cell-data/protodune{hd,vd}-wires-larsoft-*.json.bz2`.

**Official (tag `v10_20_08d00`)** —
`dunecore/dunecore/Geometry/gdml/protodunehd_v8_refactored.gdml`,
`…/protodunevd_v4_refactored.gdml`, `…/protodunevd_v5_ggd.gdml`;
`dunereco/dunereco/DUNEWireCell/{pdhd,protodunevd}/…-wires-larsoft-*.json.bz2`
(byte-identical to the toolkit copies);
`duneprototypes/.../Protodune/hd/ChannelMap/` (README, `…_visiblewires_v1.txt`,
v4/v5, mapmakers — the ±3 logic), `…/vd/ChannelMap/` (`PD2VDTPCChannelMap_v2`,
`addcrpnum.cxx`).

**Generator** — `wire-cell-python/wirecell/util/wires/gdml.py`
(`gdml-to-wires`; configs `protodunehd_v8`, `protodunevd_v4`, `protodunevd_v5`;
`build_hd_channel_map` pitch-order allocation; `assign_vd_channels` U/V drift-order
note); tests `…/util/test/test_gdml_integration_{hd,v4,v5}.py`.

**Key channel-map commits (duneprototypes)** — `6bbccf25` 2022-06-29 (intro ±3);
`76aadccd` 2024-08-05 (WIBEth visiblewires); `c8f43809` 2025-06-30 (**sign-flip**,
PR #89); `2f71e7df`/`4e92f76c` 2025-07 (VD top-CRP / v2 maps); `a67411a4`
2025-11-28 (VD TDE map).
