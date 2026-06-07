# Audit: toolkit wire geometry & channel→wire mapping vs official ProtoDUNE (HD/VD)

**Status:** investigation only — *no code/config/wire-file changes made*.
**Date:** 2026-06-07.
**Author:** geometry audit (read-only).

## 1. Scope & method

This compares the WireCell **toolkit** wire-geometry files for ProtoDUNE-HD
(PDHD) and ProtoDUNE-VD (PDVD) against the **official** DUNE channel maps /
geometry, and examines the channel→wire mapping for each, with focus on:

- the **PDHD "3-wire shift"** (a known issue, reportedly fixed upstream), and
- a suspected **stale/outdated PDVD** wire file or channel mapping.

**Files compared**

| Side | Path |
|---|---|
| Toolkit PDHD wires | `wire-cell-data/protodunehd-wires-larsoft-v1.json.bz2` |
| Toolkit PDVD wires | `wire-cell-data/protodunevd-wires-larsoft-v3.json.bz2` (also v1, v1-drift-y present, unused) |
| Official channel maps | `protodunecode/hd/ChannelMap/`, `protodunecode/vd/ChannelMap/` (symlink → `duneprototypes/.../Protodune`) |
| Wire-file generator | `wire-cell-python/wirecell/util/wires/gdml.py` (`wirecell-util gdml-to-wires`) |

`protodunecode` resolves to `duneprototypes` tag **v10_20_08d00**.

**Method** was read-only: decompressing the wire JSON and tallying
anodes/planes/channels/positions, reading the official channel-map text files and
their `mapmakers/`, and tracing git provenance of each file. **No GDML files are
present on disk** (`dunecore`/`dunereco` are not cloned here), so the one check
that would *prove* a PDHD channel↔wire shift — regenerating from the current GDML
and diffing — **cannot be run locally yet**. Conclusions below are therefore
labelled **proven** vs **strongly-evidenced (needs regenerate-and-diff)**.

---

## 2. Disambiguation: there are TWO different "±3 shifts"

These are easy to conflate; they are **different things** on different axes.

### (A) Official "visible-wire" ±3 channel shift — channel ↔ wire-position
A *cyclic permutation of offline channels within each U/V plane* so the offline
channel numbering follows where wires **emerge from under the head boards**
(the visible/active portion, matching the offline geometry's wire endpoints)
rather than their **soldered anchor points**.

- Introduced in offline channel map **v3** — `duneprototypes` commit `6bbccf25`
  (2022-06-29): *"V3 PDHD channel map with +-3 channel shifts in U and V to
  account for anchor-point vs active area shifts."*
- Carried into the WIBEth electronics era as
  `PD2HDChannelMap_WIBEth_visiblewires_v1.txt` (commit `76aadccd`, 2024-08-05).
- **Sign-flipped 2025-06-30** — commit `c8f43809`, PR #89
  (`trj_flip_threewiresign_pdhd_viswirefix_jun30_2025`):
  *"flip the sign of the +/- 3 wire shift — looks like a global offset"*; the
  mapmaker header records *"flip the sign of the +/- 3 wire correction, to fix the
  cathode-crossing inefficiency and a sign mistake."*
- Sign convention (mapmaker `MakePD2HDChannelMap_WIBEth_v1_visiblewires.C`):
  `+3` on inverted (North/Lower) APAs, `-3` on upright (South/Upper) APAs, with a
  cyclic wrap inside each plane's 800-channel block. **U and V only; W (collection)
  is unaffected** (collection wires do not wrap / have no anchor-vs-emergence
  ambiguity).

### (B) Toolkit `coh_group_shift` — channel ↔ FEMB electronics grouping
A noise-filtering grouping patch in
`cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet` (commit `9e151a19`, 2026-04-28):

```jsonnet
// coh_group_shift: cyclic offset (in offline channels) applied to U and V
// group boundaries ...  Default magnitude 3 corrects the FEMB-edge
// misassignment identified in the 027409-evt0-apa0 coherent-noise audit
// (U/V only; W is unchanged).  Sign flips by anode: +shift on anodes 0 & 2,
// -shift on anodes 1 & 3.
local shift = if n == 0 || n == 2 then coh_group_shift
              else if n == 1 || n == 3 then -coh_group_shift else 0;
local u_group(u) = std.map(function(j) n*2560 + std.mod(40*u + shift + j, 800), std.range(0,39));
local v_group(v) = std.map(function(j) n*2560 + 800 + std.mod(40*v + shift + j, 800), std.range(0,39));
```

This is **pure offline-channel arithmetic** that decides which 40 offline channels
share a coherent-noise group (a property of FEMB/WIB *electronics*, not wire
position). It was found empirically from a run-027409 correlation audit.

**Why they are different but parallel:** (A) is channel↔**wire-position**;
(B) is channel↔**FEMB grouping**. Yet both are magnitude **3**, both **U/V only**
(W exempt), and both **sign-flip by APA orientation**. That parallel motivates the
hypothesis in §5.

---

## 3. PDHD — the toolkit wire file is STALE relative to the upstream fix

**Config reference (only one):**
`cfg/pgrapher/experiment/pdhd/params.jsonnet:167`
`wires: "protodunehd-wires-larsoft-v1.json.bz2"`.

**Provenance — the heart of the problem:**

| Event | Date | Ref |
|---|---|---|
| Toolkit PDHD wire file added (single commit, never updated) | **2022-06-12** | `eb40668` "add wire geo for protodune hd" |
| Official visible-wire ±3 shift introduced (channel map v3) | 2022-06-29 | `6bbccf25` |
| WIBEth `visiblewires_v1` map | 2024-08-05 | `76aadccd` |
| **Sign-flip fix** (cathode-crossing inefficiency) | **2025-06-30** | `c8f43809` / PR #89 |

The toolkit's `protodunehd-wires-larsoft-v1.json.bz2` is from **2022-06-12 and was
never regenerated**. It therefore predates the *entire* visible-wire correction
chain — most importantly the **2025-06-30 sign flip**.

Upstream confirms this file is *meant* to be the GDML conversion output:
`wire-cell-python/.../test_gdml_integration_hd.py` uses
`dunecore/.../Geometry/gdml/protodunehd_v8_refactored.gdml` as input and
`dunereco/.../pdhd/protodunehd-wires-larsoft-v1.json.bz2` as the **reference** —
i.e. the same filename the toolkit ships. If the current `v8` GDML embeds the
corrected (post-2025-06-30) channel map, the 2022 file no longer matches it.

**Structure (decoded from the file):** 4 APAs / 8 faces / 24 planes / 22 208 wire
segments / 44 416 points; per APA U=V=800 channels, W=960 (10 240 channels total).
Within a plane, channel-by-position is **perfectly sequential**.

> **Important caveat (why this isn't yet *proven*):** internal sequentiality does
> **not** by itself reveal shift (A). The ±3 is the *offset between* the
> anchor-point ordering and the wire-emergence ordering; *each* ordering is
> internally sequential. Discriminating the two requires comparing
> channel→**physical position** against position truth (the GDML wire-emergence
> endpoints). The collection (W) plane and the per-plane channel *set* are
> identical under either convention, so neither can reveal it — **only U/V
> channel↔position can.**

**Conclusion (strongly-evidenced; not locally proven):** the PDHD wire file's U/V
channel→wire association reflects the **2022 convention** and almost certainly does
**not** carry the corrected post-2025-06-30 visible-wire mapping. The likely
observable symptom of a residual ±3 mismatch is exactly what PR #89 was fixing —
**cathode-crossing inefficiency / track-matching offsets in U/V**. Definitive
confirmation = regenerate from the current `protodunehd_v8` GDML and diff U/V
channel↔position (see §6–§7).

---

## 4. PDVD — channel-count-correct, but built from v4 geometry while v5 exists

**Config reference (only one):**
`cfg/pgrapher/experiment/protodunevd/params.jsonnet:180`
`wires: "protodunevd-wires-larsoft-v3.json.bz2"` (the `pdvd/`, `pdvd_sim/`
symlink trees carry no wire reference of their own).

**v1 → v3 reconciliation (decoded):**

| | v1 (2022-10-18) | v3 (2026-04-21, in use) |
|---|---|---|
| anodes | 16 (idents 110–117, 120–127) | **8 (idents 0–7)** |
| total channels | 12 304 | **12 288 (0–12 287)** |
| U plane (ident 0) | 3 816 ch | 3 808 ch |
| V plane (ident 1) | 3 816 ch | 3 808 ch |
| W plane (ident 2) | 4 672 ch | 4 672 ch (identical) |

v3's **12 288** channels (0–12 287) **exactly match the official scheme**
(`vd/ChannelMap/mapmakers/addcrpnum.cxx`: 0–3071 → CRP5, 3072–6143 → CRP4,
6144–9215 → CRP2, 9216–12287 → CRP3). The 16→8 anode renumber and the 16 fewer
U/V channels (8 per plane) are a **convention/numbering correction, not missing
wires** — W is byte-for-byte identical in count, and the channel range collapses
from v1's spurious 12 304 onto the official 12 288.

> **So v3 is the cleaner, current-*convention* file — it is *not* obviously
> broken.** This refines the user's "may have a problem / not latest" suspicion:
> the issue is not that v3 is malformed, but a **geometry-version** question (below).

**The version-axis confusion (the likely "not the latest version" issue):**
the wire-FILE name "v3" is a *wire-cell-data file version*, **distinct from the
GDML *geometry* version**. Upstream VD integration tests show:

- `test_gdml_integration_v4.py`: GDML `protodunevd_v4_refactored.gdml` →
  **reference** `protodunevd-wires-larsoft-v3.json.bz2` (the toolkit's file).
- `test_gdml_integration_v5.py`: GDML `protodunevd_v5_ggd.gdml` → **no reference
  wire file** yet. A `protodunevd_v5` generator config already exists in
  `gdml.py` (its U/V wire LVs carry a mandatory `_N` suffix the v4 patterns can't
  match — i.e. a genuinely different geometry).

⇒ **The toolkit PDVD wires derive from `v4` geometry; a `v5`-geometry wire file
has not been produced/adopted.** If `protodunevd_v5_ggd.gdml` is the intended
production geometry, the toolkit is one geometry generation behind.

**Channel-map currency:** official VD channel-map activity is mostly *before*
v3's 2026-04-21 commit (top-CRP `PD2VD…ChannelMap_v2`, 2025-07-10/16; TDE DAQ map
`PD2VDTopTDEChannelMap_v2`, 2025-11-28). Because v3 (2026-04) post-dates the v2
map, it **likely** incorporates it — but this is **unverified** against
`PD2VDTPCChannelMap_v2` directly, and the Nov-2025 TDE map is a separate DAQ
readout path.

**Conclusion:** v3 is channel-count-correct and current-*convention*, but is
pinned to **v4 geometry**. Recommend (a) confirming whether **v5** is the intended
production geometry and, if so, generating/adopting a v5 wire file, and
(b) verifying v3's channel↔wire against `PD2VDTPCChannelMap_v2`.

---

## 5. Cross-cutting hypothesis (well-supported — flagged as hypothesis, not asserted)

The 2026-04 `coh_group_shift` NF patch **(B)** and the PDHD wire-file staleness
**(A)** may share **one root cause**: the reconstructed data being processed
(run 027409) carries the **post-2025-06-30 visible-wire offline-channel
convention**, while the toolkit's **2022 wire file** *and* the original
coherent-noise grouping both assumed the **pre-fix** convention. A relabel of
offline channels by ±3 shifts *both* the channel↔wire-position map (A) *and* the
channel↔FEMB grouping (B) by the same amount.

Converging evidence:
- identical magnitude (**3**),
- identical scope (**U/V only; W exempt**),
- identical **APA-orientation sign-flip** (inverted vs upright APAs).

**If true**, regenerating the wire file from the current GDML may make
`coh_group_shift` redundant, or require re-tuning it (the two corrections should
not be applied twice). This is the single most valuable thing to confirm, and it
is exactly what the regenerate-and-diff in §6–§7 settles. **Stated as a hypothesis
to test — not a proven claim.** (Caveat: the empirical FEMB-edge audit that
motivated (B) is independent evidence the grouping needed a shift *for that data*;
that does not by itself pin the cause to the wire-file convention.)

---

## 6. Recommended fixes (no change made here)

1. **Clone `dunecore` + `dunereco`** (matching the duneprototypes
   `v10_20_08d00` tag) so the official GDMLs and reference wires are on disk:
   - `dunecore/dunecore/Geometry/gdml/protodunehd_v8_refactored.gdml`
   - `dunecore/dunecore/Geometry/gdml/protodunevd_v5_ggd.gdml` (and `…_v4_…`)
   - `dunereco/dunereco/DUNEWireCell/{pdhd,protodunevd}/…-wires-larsoft-*.json.bz2`
2. **PDHD — regenerate and diff** the U/V channel↔position mapping:
   ```bash
   wirecell-util gdml-to-wires -d protodunehd_v8 \
       -o protodunehd-wires-larsoft-regen.json.bz2 \
       <dunecore>/Geometry/gdml/protodunehd_v8_refactored.gdml
   ```
   Compare against the shipped `…-v1.json.bz2` (see §7). If U/V channel↔position
   differs by exactly the ±3 APA-orientation signature, **replace the wire file**
   and **re-evaluate `coh_group_shift`** (likely redundant or needs re-tuning).
3. **PDVD — decide v4 vs v5.** Confirm with the DUNE geometry owners whether
   `protodunevd_v5_ggd.gdml` is the production geometry. If yes, generate and adopt
   a v5 wire file:
   ```bash
   wirecell-util gdml-to-wires -d protodunevd_v5 \
       -o protodunevd-wires-larsoft-v5.json.bz2 \
       <dunecore>/Geometry/gdml/protodunevd_v5_ggd.gdml
   ```
   and point `protodunevd/params.jsonnet:180` at it. Independently, verify the
   current v3 channel↔wire against `PD2VDTPCChannelMap_v2`.
4. Any wire-file swap is a **geometry change** — gate behind the usual
   bit-identical-by-default discipline and rebaseline imaging/clustering snapshots.

---

## 7. Verification plan (the discriminating test)

The ±3 shift (A) is a **cyclic permutation within each U/V plane only**; the W
plane and the per-plane channel *set* are invariant. So:

- **Signal:** for each U/V plane, build `channel → rank-by-pitch-position` from
  both the shipped file and the regenerated file; the difference should be a
  constant **±3** (cyclic, mod 800), with **+3 on inverted APAs / −3 on upright
  APAs** — the visible-wire signature.
- **Null controls:** the **W plane** must show **0** difference; the per-plane
  **channel set** must be identical (only the channel↔position *pairing* moves).
- A clean ±3-only-in-U/V result **proves** the staleness; a zero result would mean
  the 2022 file already matched (it would also then fail to explain (B)).

Position truth comes from the GDML wire endpoints (head/tail in the regenerated
JSON `points`), so this requires step 6.1 first.

---

## 8. Reference appendix

**Toolkit**
- `cfg/pgrapher/experiment/pdhd/params.jsonnet:167` — PDHD `wires:`.
- `cfg/pgrapher/experiment/protodunevd/params.jsonnet:180` — PDVD `wires:`.
- `cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet` — `coh_group_shift` (B).
- `wire-cell-data/protodunehd-wires-larsoft-v1.json.bz2` — git `eb40668`
  (2022-06-12).
- `wire-cell-data/protodunevd-wires-larsoft-v3.json.bz2` — git `baeb475`
  (2026-04-21); v1 `08fcfca` (2022-10-18); v1-drift-y `f7c7e98` (2023-03-28).
- `wire-cell-python/wirecell/util/wires/gdml.py` — `gdml-to-wires`, configs
  `protodunehd_v8` / `protodunevd_v4` / `protodunevd_v5`.
- `wire-cell-python/wirecell/util/test/test_gdml_integration_{hd,v4,v5}.py` —
  GDML→wires references.

**Official (`protodunecode` → duneprototypes v10_20_08d00)**
- `hd/ChannelMap/README.txt` — v0–v5 + WIBEth map history, APA numbering.
- `hd/ChannelMap/PD2HDChannelMap_WIBEth_visiblewires_v1.txt` (+ mapmaker
  `MakePD2HDChannelMap_WIBEth_v1_visiblewires.C`) — the ±3 visible-wire logic.
- `vd/ChannelMap/PD2VDTPCChannelMap_v2.txt`, `PD2VDTopTPCChannelMap_v2.txt`,
  `PD2VDTopTDEChannelMap_v2.txt`; `mapmakers/addcrpnum.cxx` — CRP↔channel ranges.

**Key commits (duneprototypes)**
- `6bbccf25` 2022-06-29 — introduce ±3 (U/V) channel-map v3.
- `76aadccd` 2024-08-05 — WIBEth `visiblewires_v1`.
- `c8f43809` 2025-06-30 — **flip sign of ±3** (PR #89); merge `a8a2dcdd` 2025-07-07.
- `2f71e7df` 2025-07-10 / `4e92f76c` 2025-07-16 — VD top-CRP map / v2 maps.
- `a67411a4` 2025-11-28 — VD TDE map.
