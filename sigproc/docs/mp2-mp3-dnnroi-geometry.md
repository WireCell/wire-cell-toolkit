# MP2 / MP3 multi-plane ROIs for DNNROI: geometry, cushions, and the sim↔reco map

This note explains the **MP2** and **MP3** "multi-plane protection" ROIs that Signal
Processing produces as input planes for **DNNROI**, answering three questions:

1. Which **wire-geometry** information the algorithms use.
2. How much **cushion** they allow in the **time** dimension (the rebin-10 vs rebin-4
   question) and in the **wire** dimension.
3. Whether the **same wire geometry is used in simulation and in reconstruction**, and
   how the PDHD "3-channel" channel-grouping fix interacts (or not) with the
   simulation's channel→wire mapping.

All line references are to the source as of this writing.

---

## 0. What MP2/MP3 are and where they live

The multi-plane protection ROIs are produced inside `OmnibusSigProc` during ROI
refinement, gated by the `use_multi_plane_protection` flag
(`sigproc/src/OmnibusSigProc.cxx:1968-1984`):

```cpp
if (m_use_multi_plane_protection) {
    for (const auto& f : m_anode->faces()) {
        // mp3: 3-plane protection based on cleanup ROI
        roi_refine.MP3ROI(iplane, m_anode, f, m_roi_ch_ch_ident, roi_form,
                          m_mp_th1, m_mp_th2, m_mp_tick_resolution, 2, 2,
                          m_plane2layer, m_MP_feature_val_method);
        // mp2: 2-plane protection based on cleanup ROI
        roi_refine.MP2ROI(iplane, m_anode, f, m_roi_ch_ch_ident, roi_form,
                          m_mp_th1, m_mp_th2, m_mp_tick_resolution, 2, 2,
                          m_plane2layer, m_MP_feature_val_method);
    }
    save_mproi(*itraces, mp3_roi_traces, iplane, roi_refine.get_mp3_rois());
    save_mproi(*itraces, mp2_roi_traces, iplane, roi_refine.get_mp2_rois());
}
```

The algorithms themselves:

- `ROI_refinement::MP3ROI` — `sigproc/src/ROI_refinement.cxx:2598-2752`
- `ROI_refinement::MP2ROI` — `sigproc/src/ROI_refinement.cxx:2755-2896`
- signatures in `sigproc/src/ROI_refinement.h:35-43`.

Both run **after** the loose/tight deconvolution and the `CleanUpROIs` /
`generate_merge_ROIs` step, i.e. they operate on the *cleanup* ROIs of all three
planes. They are run **per anode face** and **per target plane**.

**MP3 (3-plane coincidence) — protects.** For a given target plane it looks for
real cleanup-ROI activity on **both** other planes at the **same time bin**, projects
the two reference-plane wires to a predicted target-plane wire, and—if a target-plane
ROI already exists near that wire—marks it protected. It does **not** invent new ROIs;
it only shields existing ones from later erosion. Output: `proteced_rois`
(`get_mp3_rois()`), a multimap keyed by `{chid, start_bin}`.

**MP2 (2-plane prediction) — creates.** Same geometric projection, but it stamps a
protection ROI on the predicted target-plane wire(s) **even if no target-plane ROI was
found there**. This recovers signal on a plane that locally saw nothing (e.g. a dead
or inefficient region) but where the other two planes agree a track crossed. Output:
`mp_rois` (`get_mp2_rois()`), keyed by `{roichid, 0}`.

Both are saved as trace tags `mp3_roi%d` / `mp2_roi%d`
(`OmnibusSigProc.cxx:2136-2137`) and consumed by DNNROI as input channels alongside
the loose-deconvolution waveform. The DNNROI input tag list is
`['loose_lf', 'mp2_roi', 'mp3_roi']` (with optional `tight_lf, decon_charge, gauss`
for the 6-channel models), see `cfg/layers/low/dnnroi.jsonnet:25-27` and the
per-experiment `dnnroi_pp.jsonnet`.

Note the target plane: in both algorithms `if (plane == 2) return;` — only the two
**induction** planes (U, V) get MP2/MP3 protection. The collection plane (W) is left
to its own deconvolution (it passes through DNNROI as the `gauss` passthrough).

---

## 1. Which wire-geometry information is used

The single source of truth is the **`AnodePlane` built from a `WireSchema` geometry
file**. For PDHD that file is `protodunehd-wires-larsoft-v1.json.bz2`
(`cfg/pgrapher/experiment/pdhd/params.jsonnet:167`), loaded as a `WireSchemaFile`
component and attached to each `AnodePlane` (`cfg/layers/low/anodes.jsonnet`,
`cfg/pgrapher/common/tools.jsonnet:82-85`).

Two distinct pieces of geometry are used:

### (a) Channel → wire-index lookup

For each ROI, the algorithm maps the ROI's internal channel id back to the WCT
channel ident, asks the anode for that channel's wires, keeps only the wires on the
current face, and reads each wire's **index** (the wire-attachment-number / pitch
index within its plane):

```cpp
auto chid = map_ch.at(roi->get_chid());
auto ch   = anode->channel(chid);
for (auto wire : ch->wires()) {
    if (face->which() != wire->planeid().face()) continue;
    auto pit_id = wire->index();        // <-- wire/pitch index from WireSchema
    ...
}
```
(MP3: `ROI_refinement.cxx:2619-2626`; MP2 builds the inverse `map_wireid_roichid`
the same way at `:2775-2786`, and re-reads it for the per-ROI loop at `:2799-2806`.)

### (b) RayGrid pitch geometry (the cross-plane projection)

The actual "deduce one plane's ROI from the other two planes" step is a **RayGrid**
calculation. Each contributing wire becomes a RayGrid coordinate
`{grid = wire index, layer = iplane2layer[iplane] + nbounds_layers}`, and the crossing
of the two reference-plane rays is projected onto the target plane:

```cpp
auto pitch = face->raygrid().pitch_location(c1, c2, layer);  // physical pitch of the U×V crossing
auto index = face->raygrid().pitch_index(pitch, layer);      // -> target-plane wire index
```
(MP3: `ROI_refinement.cxx:2685-2686`; MP2: `:2859-2860`.)

`pitch_location`/`pitch_index` are pure RayGrid geometry — they encode the **wire
pitch and the three wire angles** of the detector. So the geometry information actually
used is: the per-plane wire pitch, the wire orientations (via RayGrid), and the
channel↔wire-index correspondence (via WireSchema). No drift/time geometry enters;
the cross-plane match is purely transverse (within a common time bin).

### RayGrid layering and the `plane2layer` remap

- `nbounds_layers = 2` (hardcoded at the call site). This is the RayGrid convention:
  layers 0 and 1 are the two **bounding rays** of the sensitive volume; the U/V/W wire
  planes occupy layers 2, 3, 4. Hence `coord.layer = iplane2layer[iplane] + 2`.
- `plane2layer` maps a plane index to its RayGrid layer. Default is `[0,1,2]`
  (identity). **PDHD APA0 uses `[0,2,1]`** (`cfg/pgrapher/experiment/pdhd/sp.jsonnet:112`,
  GitHub issue #322 / the APA0 V-plane anomaly), which swaps the V and W layers. Because
  MP2/MP3 build `coord.layer` from this very vector, the swap is honoured **inside the
  multi-plane projection** as well as in the internal channel-map build
  (`OmnibusSigProc.cxx:254`). This is a *plane/layer bookkeeping* remap; it does not
  change which physical wire a channel reads (see §3).

---

## 2. Cushion in the time and wire dimensions

### Time cushion — this is the "rebin 10 vs rebin 4" question

The relevant knob is `tick_resolution`, fed from the config parameter
`mp_tick_resolution`. The ROI contents are walked in blocks of `tick_resolution` ticks
and each block is reduced to a single feature value:

```cpp
for (int tick = roi->get_start_bin()/tick_resolution;
     tick <= roi->get_end_bin()/tick_resolution; ++tick) {
    int content_id = tick*tick_resolution - roi->get_start_bin();
    auto feat_val = feature_val(roi->get_contents(),
                                content_id, content_id + tick_resolution, method);
    ...
}
```
(MP3: `ROI_refinement.cxx:2630-2633`; MP2: `:2809-2812`.)

A protected region is therefore quantized into **`tick_resolution`-wide blocks**
(`sta = tick*tick_resolution`, `end = sta + tick_resolution`). MP3 additionally
**merges** consecutive protected blocks on the same channel when they are separated by
no more than one `tick_resolution` (`:2711-2713`):

```cpp
if (sta - last_end <= tick_resolution) last_iter->second.second = end;  // extend
else                                   proteced_rois.insert({key,{sta,end}});  // new block
```

So the effective **time cushion is one `tick_resolution` block** — the protection is
granular to that many ticks, and gaps up to that size are bridged.

**Both rebin values are live in the repository**, depending on the pipeline:

| Value | Where |
| --- | --- |
| `mp_tick_resolution: 4` | PDVD (`protodunevd/wct-sim-*`, `protodunevd/sp.jsonnet`), PDHD sim-fans / deposplat (`pdhd/wct-sim-fans.jsonnet`, `pdhd/wcls-sim-drift-deposplat.jsonnet` active variant) |
| `mp_tick_resolution: 10` | PDHD LArSoft reco / DNNSP path (`pdhd/wcls-rawdigit-dnnsp.jsonnet`, `pdhd/wcls-nf-dnnsp-img.jsonnet`, `pdhd/wcls-sim-drift-deposplat.jsonnet` default) |

The C++ defaults are inconsistent on purpose-of-history: the config member default is
`m_mp_tick_resolution{4}` (`OmnibusSigProc.h:261`), while the *function-signature*
default in `ROI_refinement.h` is `10`. The value actually used is always the one passed
from the jsonnet config, so **check your pipeline's `mp_tick_resolution`** — PDHD
production reco runs at **10**, PDVD and the PDHD sim/active path at **4**.

> Do not confuse `mp_tick_resolution` with DNNROI's own time downsampling. DNNROI
> separately rebins **all** of its input planes by `tick_per_slice` before the network
> (`pytorch/src/DNNROIFinding.cxx:371`, also 4 or 10 depending on detector). That is a
> different parameter applied at a different stage.

### Wire cushion — effectively ±1 wire

The knob is `wire_resolution`, hardcoded to `2` at both call sites
(`OmnibusSigProc.cxx:1973-1975`). Despite the value "2", the **effective transverse
cushion is ±1 wire** on both code paths:

- **MP3** protects a target ROI when its pitch is within `wire_resolution` of the
  predicted index: `if (abs(pitch - index) < wire_resolution)` → `abs(...) < 2` →
  **±1 wire** (`ROI_refinement.cxx:2691`).
- **MP2** stamps protection across `index - wire_resolution + 1 ... index + wire_resolution - 1`,
  i.e. `index-1 .. index+1`, a 3-wide window → **±1 wire** (`:2863-2864`).

So: **time cushion ≈ one `mp_tick_resolution` block (4 or 10 ticks); wire cushion ≈
±1 wire.**

### Amplitude gating (for completeness)

Two thresholds shape what counts as activity, via a per-block feature value
(`feature_val`, `:2572-2594`; either a 3-point max or a global max, set by
`MP_feature_val_method`):

- `mp_th2` (default 500) — minimum to register a plane's activity at a tick.
- `mp_th1` (default 1000) — stricter bar to qualify a plane as a **reference** for the
  cross-plane crossing (`:2640-2656`, `:2814-2830`).

---

## 3. Same wire geometry in sim and reco? And the PDHD "3-channel" interaction

### The channel→wire geometry is identical in simulation and reconstruction — yes

Both ends of the chain read the **same `AnodePlane` / `WireSchema`** built from the
same wires file:

- **Simulation (Depo→ADC):** `DepoTransform` is constructed from `tools.anodes[n]`
  (`cfg/pgrapher/experiment/pdhd/sim.jsonnet`), i.e. the same `AnodePlane` objects.
- **Reconstruction (ADC→NF→SP→imaging):** NF, SP and imaging all take the same
  `anode` (`pdhd/nf.jsonnet`, `pdhd/sp.jsonnet`, `pdhd/img.jsonnet`), and MP2/MP3 read
  channel→wire from that anode (§1).

These `AnodePlane`s come from one `WireSchemaFile` (`anodes.jsonnet`,
`tools.jsonnet:82-85`). There is no second, alternate channel↔wire map anywhere in the
chain. So the geometry MP2/MP3 and imaging rely on is exactly the geometry simulation
projected charge onto.

### The PDHD "3-channel" fix is a noise-grouping shift, **not** a channel→wire change

The PDHD "3-channel confusion" is the `coh_group_shift` correction in
`cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet:39-56`. It applies a cyclic offset of
**3 offline channels** to the U- and V-plane **coherent-noise group boundaries**
(`+3` on APA 0 & 2, `-3` on APA 1 & 3; W unchanged):

```jsonnet
local shift = if n == 0 || n == 2 then  coh_group_shift   // default 3
              else if n == 1 || n == 3 then -coh_group_shift else 0;
local u_group(u) = std.map(function(j)
    n*2560 + std.mod(40*u + shift + j, 800), std.range(0,39));
```

This realigns the 40-channel coherent-subtraction groups with the true **FEMB edges**
(identified in the 027409-evt0-apa0 coherent-noise audit). Crucially, it changes only
**which channels are averaged together** when estimating and subtracting coherent
noise. It does **not** alter the WireSchema channel↔wire correspondence at all — the
geometry map of §1 is untouched.

### Therefore there is no interaction with the simulation's channel→wire mapping

- The shifted groups live in `femb-negpulse-groups-shifted_v2.jsonnet` and are consumed
  **only in noise filtering** (`pdhd/nf.jsonnet:82`, via the channel-noise database
  built by `chndb-base.jsonnet`). A repo-wide search finds **no simulation config**
  that references them.
- PDHD simulation noise is **incoherent-only** (per the PDHD noise-RMS validation
  study: post-NF data is compared against raw incoherent-only sim). There is no grouped
  coherent noise *generated* in sim that the shift could be inconsistent with.

In short, `coh_group_shift` is a **reconstruction-side, real-data coherent-subtraction
correction** that operates on channel *grouping*. It is orthogonal to the channel→wire
geometry that simulation and reconstruction share. Simulation is unaffected, and the
MP2/MP3 / imaging geometry is unaffected; the shift only changes how NF cleans real
data before that shared geometry is ever consulted.

### Aside: the separate APA0 `plane2layer [0,2,1]` quirk

This is a *different* PDHD peculiarity from the 3-channel grouping shift. The APA0
`plane2layer: [0,2,1]` (`pdhd/sp.jsonnet:112`, issue #322) reorders the **RayGrid layer
bookkeeping** in SP (and thus in MP2/MP3 via `coord.layer`, §1). It does **not** change
which physical wire a channel reads. Simulation's `DepoTransform` projects charge
straight through the WireSchema geometry and never consults `plane2layer`, so it is
consistent with this remap by construction — the remap only affects how SP labels the
planes internally, not the underlying channel→wire truth.

---

## Quick reference

| Question | Answer |
| --- | --- |
| Geometry used | WireSchema channel→wire index + RayGrid pitch (wire pitch & the 3 wire angles); `pitch_location`/`pitch_index` for the U×V→target projection |
| Target planes | U and V only (`plane==2` returns); W passes through |
| MP3 vs MP2 | MP3 *protects* existing target ROIs at a 3-plane coincidence; MP2 *creates* target ROIs predicted from the other 2 planes |
| Time cushion | one `mp_tick_resolution` block; **4** (PDVD, PDHD sim/active) or **10** (PDHD LArSoft reco). Distinct from DNNROI's `tick_per_slice` |
| Wire cushion | effectively **±1 wire** (`wire_resolution=2` → `abs<2` / a 3-wide window) |
| Same geometry sim vs reco? | **Yes** — one `AnodePlane`/`WireSchema` from one wires file, used by DepoTransform and by NF/SP/imaging alike |
| PDHD 3-channel fix | `coh_group_shift=3` shifts **coherent-noise group membership** in NF only; does **not** touch the channel→wire map and has **no** simulation counterpart (sim noise is incoherent-only) |

### Source map

| Item | Location |
| --- | --- |
| MP3 / MP2 algorithms | `sigproc/src/ROI_refinement.cxx:2598-2752` / `:2755-2896` |
| Signatures | `sigproc/src/ROI_refinement.h:35-43` |
| Call site + flags | `sigproc/src/OmnibusSigProc.cxx:1968-1984`, tags at `:2136-2137` |
| Config defaults | `sigproc/inc/WireCellSigProc/OmnibusSigProc.h:255-261` |
| Channel-map / plane2layer build | `sigproc/src/OmnibusSigProc.cxx:240-289` |
| DNNROI input tags / downsample | `cfg/layers/low/dnnroi.jsonnet:25-27`, `pytorch/src/DNNROIFinding.cxx:371` |
| PDHD wires file | `cfg/pgrapher/experiment/pdhd/params.jsonnet:167` |
| PDHD APA0 plane2layer | `cfg/pgrapher/experiment/pdhd/sp.jsonnet:112` |
| PDHD coh_group_shift | `cfg/pgrapher/experiment/pdhd/chndb-base.jsonnet:39-56`, `nf.jsonnet:82` |
