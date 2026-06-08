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

## 4b. PDVD follow-up: geometry vs channel-mapping, and the v4 test (2026-06-07)

Driven by a real symptom: **PDVD *data* imaging shows gaps** — fired U/V/W wires
that should cross don't, so blobs break up. That is the signature of either a
channel→wire mapping error or a wire-geometry inconsistency. We separated the two.

### Geometry: v3 is built from v4 GDML; v5 moved the bottom CRPs
- The toolkit `…-v3.json.bz2` reproduces the **v4** GDML to ~µm (matched wire-by-wire
  by exact position: 6 µm match, 7.2 mm to next-nearest → unambiguous bijection).
- The **v5** GDML (`protodunevd_v5_ggd.gdml`) repositions wires: **bottom CRPs
  (anodes 0–3) shift ~5.5 mm in W and ~2.5 mm in U/V; top CRPs (4–7) unchanged.**
  Within each CRP this is **near-rigid** (best-fit translation, residual ~0.2 mm).
  Wire **angles** (U 30°, V 150°, W 90°) and counts/wrapping are identical v3↔v5.
- Method caveat learned here: **rank/perp-matching is unreliable for wrapped
  induction.** It invented a 2.3 mm PDHD-U "shift" that **channel-identity matching
  showed was 8 µm**. Only channel-matched (or exact-position-matched) comparisons
  are trusted below. (PDHD is separately confirmed **bit-consistent with v8 GDML** —
  all 22208 wires within 8 µm.)

### The v4 file (exact, no approximation)
Built `wire-cell-data/protodunevd-wires-larsoft-v4.json.bz2` = **v3's channel
assignment + the exact v5 GDML wire positions**. Construction is exact: each v3
wire is matched to its v4-geometry counterpart by exact position (6 µm, unique),
then given the corresponding v5 endpoints (regen-v4↔regen-v5 correspond 1:1 by
ident/channel/segment). Verified: channel/ident/segment identical to v3; positions
equal v5 to 0.00000 mm and differ from v3 by up to 5.700 mm; channel↔plane still
matches the official map 100%.

### Result: geometry is NOT the gap cause
Re-imaged + re-clustered run 39324 evts 0–4 with v4. **Imaging is essentially
identical to v3** — blob counts agree to within 0–5 per anode (e.g. anode0
2681→2678, anode4 13572→13572). A near-rigid per-CRP shift *moves* each CRP's image
without filling gaps, exactly as predicted. **So updating to v5 geometry does not
close the PDVD gaps.** Bee link (v4, evts 0–4):
`https://www.phy.bnl.gov/twister/bee/set/3fa7822a-52d0-4143-8267-922e33571fbd`.

### Where that leaves it
The gap cause is **not** the wire positions (v4/v5) — it is the **within-plane
channel→wire mapping**. What is established: v3's channel→**plane** assignment and
channel count match `PD2VDTPCChannelMap_v2` **100%**, and the v1→v2 update did not
touch channel→wire (electronics-only). What remains **unverified** is the
within-plane channel→wire ordering, which **cannot be settled from files alone**
(wrapped induction + 2-D collection defeat position-rank methods; the converter
renumbers channels so there's nothing to channel-match against). The definitive
test is a **LArSoft `PD2VDChannelMapService` wire-dump** (`channel tpc plane wire
x y z`) — with it, a channel-matched comparison gives a yes/no in seconds (as it
did for PDHD). The production config remains on **v3** (v4 changes nothing); v4 is
left in `wire-cell-data` for further testing.

---

## 5. Why HD reproduces but VD doesn't

- **HD** channels reproduce exactly because `build_hd_channel_map` is the same tool
  (or direct ancestor) that made the HD reference, with a fully pinned pitch
  convention. **VD**'s `assign_vd_channels` has an acknowledged U/V drift-order
  ambiguity, and the shipped v3 was evidently produced by a different/older
  channel-assignment pass — hence the mismatch.

---

## 5b. The PDHD `coh_group_shift` ±3, and sim↔data consistency (resolved)

Follow-up question: is the channel↔wire mapping consistent between **simulation**
and **real data**, and is the ±3 in the coherent-noise removal / FEMB-saturation
tagger a wire-mapping bug or just a FEMB-grouping confusion?

**The key distinction — two independent lookups from an offline channel:**
- **channel↔wire** (geometry) — where the signal sits in 3D. Lives in the
  wire-geometry file; used by SP / imaging / reco.
- **channel↔FEMB** (electronics) — which 64 channels share a cold-electronics
  motherboard → coherent noise & saturation. **Not in the wire file**; it is the
  `n*2560 + 40*u` block arithmetic in `chndb-base.jsonnet`.
The ±3 lives **entirely in the second relation.** A correct wire mapping can still
need a ±3 in the grouping.

**Q1 — channel↔wire consistent sim↔data? YES.**
- WCT **sim and reco share the same wire file** (`params.jsonnet:167`;
  `simparams.jsonnet` inherits `files:` via `super`) — identical by construction.
- It matches the real data empirically: **run-027409 imaging closes richly**
  (6712 / 7647 / 3035 / 6287 blob points in APA0–3). A 3-wire U/V offset
  (opposite sign per APA) would destroy U∩V∩W closure and collapse those counts.

**Q2 — is the ±3 just a FEMB-grouping confusion? YES (confirmed from the map).**
For APA0 U-plane (`PD2HDChannelMap_WIBEth_*`):

| map | physical FEMB → offline channels |
|---|---|
| **electronics** (`…_electronics_v1`) | clean `[40m … 40m+39]` (FEMB10 = 0–39, FEMB20 = 400–439) |
| **visible-wire** (`…_visiblewires_v1`, post-flip) | same FEMBs, offline labels **−3** (FEMB9 = 37–76, FEMB10 wraps 797…36) |

A physical FEMB is a fixed set of 40 wires. WCT's `40*u+j` grouping assumes the
**electronics** block layout, but the data's offline numbering follows the
**visible-wire** convention — shifting the blocks by ±3 in U/V (W has no wrap → no
shift, exactly as in `chndb-base.jsonnet`). Hence `40u, 40u+1, 40u+2` carry the
*previous* FEMB's common mode — the exact run-027409 symptom. The FEMB-saturation
tagger uses the same shifted groups (`femb-negpulse-groups-shifted_v2.jsonnet`),
so it is the same effect. **A grouping issue, not a wire-mapping issue.**

**Q3 — channel↔wire mismatch sim↔data? Not operative for this data — but a
convention flag.** The APA-dependent sign (`+3` on upright anodes 0&2, `−3` on
upside-down 1&3 — a colleague's finding) is the official visible-wire structure
(mapmaker `udown = (icrate==1||icrate==3)`; WCT anode0 = `APA_P02SU` = crate 0 =
upright). Verified against git:

| convention | upright APA0 U | FEMB10 offline (e.g.) |
|---|---|---|
| **pre**-flip (before `c8f43809`, 2025-06-30) | electronics **+3** | 3–42 |
| **post**-flip (current `v10_20_08d00`) | electronics **−3** | wraps 797…36 |

`coh_group_shift` uses **+3 on upright APAs ⇒ the pre-flip convention** ⇒
**run-027409 was decoded with the pre-June-2025 channel map.** Since imaging closes
against the wire file, the **wire file is effectively on the same old (pre-flip)
convention as this data** — mutually consistent, so no wire mismatch bites and sim
(same file) is consistent too.

**Actionable caveat:** the official map was **sign-flipped 2025-06-30**. If
production data is re-decoded with the **post-flip** map, (i) the `coh_group_shift`
sign must **flip**, and (ii) the wire file must be re-confirmed — an old-convention
wire file against post-flip data *would* produce a real 3-wire U/V mismatch.

**Simulation caveat:** PDHD `wct-sim` injects **only incoherent
`EmpiricalNoiseModel` noise — no coherent FEMB noise** (no `CoherentAddNoise` /
`GroupNoiseModel`; sim NF pipes commented out). So the ±3 is **benign on sim**, and
the **FEMB grouping / saturation logic cannot be validated on this sim** — that
would require injecting coherent noise on the true (electronics) FEMB grouping.

**Residual (needs LArSoft runtime):** the *absolute* offline→wire convention of the
wire file in isolation needs a `PD2HDChannelMapService` dump; the imaging-closure +
sign evidence pins the practical answer above.

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
