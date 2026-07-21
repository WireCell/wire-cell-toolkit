# SBND reco1 art-file sources (bare-ROOT, no LArSoft)

Components: `SBNDReco1FrameSource` (IFrameSource), `SBNDReco1OpFlashSource`
(ITensorSetSource), both in `root/`.  They read SBND "reco1" art/LArSoft ROOT
files **directly with the toolkit's own ROOT** — no LArSoft, no gallery, no art —
and emit the same interchange payloads the SBND standalone chain already
consumes (`frames-dnn.tar.bz2`-style tagged frames; `opflash_apa<N>.tar.gz`-style
opflash matrix tensor sets).

## Repro block

```sh
# toolkit-dev direnv env (wire-cell, wcsonnet on PATH), then:
cd sbnd_xin
./run_reco1_dump.sh -t <tag> [reco1.root]      # all events -> input_files_reco1/extracted-<tag>/
export SBND_INPUT_DIR=$PWD/input_files_reco1/extracted-<tag>
./run_img_evt.sh data <idx>                     # unchanged downstream chain
./run_ql_evt.sh  data <idx>
```

Development files: `sbnd_xin/input_files_reco1/data_filtered_decoded_reco1-
fe6033f3-07a0-4971-cea5-16ce59269fba_eventidfiltered.root` (48 real-data
neutrino-candidate events across 12 runs 18253..18409 — mostly 18255 —
SBND 2025 fall production, BNBLight stream, ~254 MB) and its
`*_eventidfiltered_frameshift.root` companion: the same 48 events with the
sbndcode `FrameShift` producer run over them (process `FRAMESHIFT`), adding
the `sbnd::timing::FrameShiftInfo` + `TimingInfo` products (~34 KB).

## How reading works (and its two hard-won gotchas)

Art files embed complete ROOT `StreamerInfo` for every stored class.
`root/inc/WireCellRoot/SBNDReco1Products.h` carries hand-trimmed **mirror
classes** (derived from `TFile::MakeProject` output on the file above) whose
names/members match that streamer info; `root/dict/LinkDef.h` generates their
dictionary into `libWireCellRoot` (the `root/` package already had the
rootcling hook).  With the dictionary loaded, plain typed reads work:

```cpp
art::Wrapper<std::vector<recob::Wire>>* w = ...;
branch->SetAddress(&w); branch->GetEntry(entry);   // w->obj is the vector
```

Gotchas (encoded in `root/src/SBNDReco1Reader.h` — do not "simplify" away):

1. **Never `tree->GetEntry(i)`.**  The Events tree has ~2500 branches including
   DAQ fragment products with no dictionary; deserializing everything
   segfaults.  Always per-branch `branch->GetEntry(i)`.
2. **Read the top-level wrapper branch** (`"<product>."`) with an
   `art::Wrapper<...>` object.  Setting a typed `std::vector` address on the
   `".obj"` member sub-branch silently reads nothing.

Products read (data reco1 branch names, all configurable):

| product | branch | used for |
|---|---|---|
| post-SP DNN wires | `recob::Wires_sptpc2d_dnnsp_Reco1.` | frame traces (tag `dnnsp`) |
| bad-channel masks | `ints_sptpc2d_badmasks_Reco1.` | cmm `"bad"` ([ch, first, last] triplets) |
| Wiener summary | `doubles_sptpc2d_wienersummary_Reco1.` | tagged-trace summary |
| OpFlash per TPC | `recob::OpFlashs_opflashtpc{0,1}__Reco1.` | opflash matrix |
| PTB decode | `raw::ptb::sbndptbs_ptbdecoder__DECODE.` | FrameShift port |
| TDC decode | `sbnd::timing::DAQTimestamps_tdcdecoder__DECODE.` | FrameShift port |
| DAQ header | `artdaq::detail::RawEventHeader_daq_RawEventHeader_*` (per-EVB scan) | FrameShift port |
| FrameShift product | `sbnd::timing::FrameShiftInfo_frameshift__FRAMESHIFT.` (only in `*_frameshift.root`) | `caf_offset_mode="product"` |
| run/subrun/event | `EventAuxiliary`, Runs-tree `RunAuxiliary` | idents, frame time |

## Output conventions (compatibility contract)

Pinned against yuhw's LArSoft dumps (`sbnd_xin/input_files/wcls-frame-dump.*`,
`wcls-flash-dump.*`) and larwirecell `CookedFrameSource`/`OpFlashSource`:

- Frames: dense float traces, `charge = frame_scale(50) x wire value`, nticks
  3427, tick 500 ns, tbin 0, identity channel numbers, frame ident = event
  number.  Frame time replicates the CookedFrameSource convention (event vs
  run-begin wall-clock difference in bare seconds — an upstream quirk kept for
  drop-in compatibility; nothing consumes it as physical time).
- Opflash tensor: f8 `[nflash, 1+312]`, col 0 = `OpFlash::Time()` (us) ×
  `units::microsecond` = ns (trigger-relative), cols 1..312 PE per OpDet
  (padded from the file's 305 entries with zeros).  Tensor-set ident = event
  number; md `run/subrun/event/nchan`.

## FrameShift port and `frame_apply_at_caf`

Newer SBND dumps carry a per-event scalar `frame_apply_at_caf` (ns) in the
opflash tensor-set metadata; `Aux::FlashTensorToOpticalPCs` (knob
`correct_flash_time`, default on) adds it to every flash/light time.  Upstream
this comes from sbndcode's `FrameShift` producer (SBN DocDB 46654) +
`FrameShiftInfo::FrameApplyAtCaf()` — the reco1 files predate that producer, so
`SBNDReco1OpFlashSource` ports the derivation
(`root/src/SBNDReco1FrameShift.h`, prototype
`sbndcode/Timing/FrameShift/FrameShift_module.cc`, constants from
`frameshift_sbnd_data.fcl`): DAQ-header timestamp − 367 us NTB correction; PTB
HLT word unmask (timestamps × 20 ns); TDC ch4 = ETRIG, ch0 = CRT T1, ch2 = RWM;
global ETRIG reference TDC > PTB > raw; closest-timestamp selection; stream
from ETRIG HLT code (beam [1,2,16] / offbeam [3,4,17] / xmuon [5,14,15]);
frames gate/etrig/crtt1 with the fcl TDC-vs-PTB shifts.

Knob `caf_offset_mode` (C++ default `"none"` = key omitted, byte-identical to
the older no-key dumps):

- `"product"` — read `FrameShiftInfo::fFrameApplyAtCaf` from the file's
  `frameshift_product` branch (throws if absent — requires a
  `*_frameshift.root` input).  **Authoritative.**
- `"auto"` — write `frame_etrig - frame_default` (ns) from the ported
  derivation.  Approximation for pre-FrameShift files.
- `"override"` — write `caf_offset_override` verbatim.

**The reduction question, answered (2026-07-21).**  A FRAMESHIFT-processed
companion of the dev file resolved what `FrameApplyAtCaf()` actually returns:
on all 48 events (`fTimingType = 0`, SPEC-TDC-ETRIG-decoded frames)
`fFrameApplyAtCaf == fFrameTdcRwm` **exactly** — the shift from the decoded
frame to the SPEC-TDC RWM signal (per-spill beam arrival; the PTB
`fFrameHltBeamGate` tracks it within ±20 ns).  The `auto` candidate
`frame_etrig − frame_default` (256-ns-quantized, 0–2560 ns) underestimates it
by 43–482 ns (mean 262): right scale — hence the passed beam-window
validation — wrong formula.  The exact value is *not derivable* from the
DECODE-stage products `auto` reads: the module's decoded-frame reference sits
a few hundred ns below the nearest SPEC-TDC ch3 (FTRIG) entry and comes from
PMT-decoder products the extraction does not consume.  So `auto` stays a
±0.5 us fallback for pre-FrameShift files; use `product` whenever the
producer has run.

## Verification — product mode (2026-07-21, `*_frameshift.root`)

- Mirror-class reads (`FrameShiftInfo`, `TimingInfo`, both split-level-1
  single-object wrappers): values identical to a dictionary-free
  `TTree::Scan` on all 48 events; the port's TDC RWM/ETRIG picks equal the
  module's `TimingInfo` timestamps u64-exactly on 48/48.
- Extraction `-caf product` (tag `extracted-2025fall-48evt-fsprod`):
  `frame_apply_at_caf` metadata equals `fFrameApplyAtCaf` on 48/48 events in
  both APA tarballs.
- Physics: raw in-time flash median −0.715 us; corrected median +1.28 us with
  **45/48 in the +0.3..1.9 us beam window** (the `auto` approximation gave
  42/48).  All 48 events are neutrino candidates.  Outliers: 214469/111412
  marginal (1.94/2.04 us), 131357 (−0.40 us, min-|t| flash likely a cosmic).

## Verification (2026-07-21)

- Format identity: extracted archives reproduce the reference member scheme
  exactly (names, shapes `(11276, 3427)` f4 / `(nflash, 313)` f8, tickinfo
  `[t, 500, 0]`, chanmask triplets — first bad channel 546 matches the
  existing data sample).
- Determinism: two independent full 48-event extractions member-hash identical
  (`hash_archive.py`): `frames-dnn.tar.bz2` rollup `46ff819f…` (240 members),
  opflash rollups `554c6924…`/`307de3f6…`.
- End-to-end: evt 256587 through unchanged `run_img_evt.sh` +
  `run_ql_evt.sh` (imaging 46 s, QL 32 s, `mabc-all-apa.zip` produced; QL log
  shows `frame_apply_at_caf offset=2048 ns` applied, 15+11 flashes matched).
- Existing-chain A/B gate PASS: data evt 686 joint QL chain (the only standard
  SBND job that loads `WireCellRoot`), baseline lib vs new lib under
  `setarch x86_64 -R`: `mabc-all-apa.zip` member-hash identical
  (`4b006453…`, 5 members; snapshots `ab_A1`/`ab_B1`/`ab_B2`, B1=B2 confirms
  determinism).  The dictionary embedded in `libWireCellRoot` does not alter
  any existing output.  Re-gated after the FrameShiftInfo/TimingInfo mirror +
  `product`-mode additions: same evt 686 rerun (setarch, fresh work root)
  reproduces `4b006453…` exactly.
- `root/` has no doctest binary (`wcdoctest-root` does not exist; package
  carries only legacy `test_*` tools) — nothing additional to run.

## Known limitations

- Real-data reco1 only (the MC counterpart uses `simtpc2d` labels — pass
  `wire_product`/`badmask_product`/`summary_product` accordingly).
- No `raw::RawDigits` survive the event filter, so NF/SP cannot be re-run from
  these files; the chain necessarily starts at imaging.
- Mirror classes are pinned to the streamer layouts listed in
  `SBNDReco1Products.h` (Wire v19, OpFlash v19, OpHit v15...).  A future
  LArSoft schema change will need a member-level review against the new
  StreamerInfo (ROOT evolves by member name; additive changes are usually
  transparent).
- Loading `libWireCellRoot` inside a LArSoft job now registers duplicate
  dictionaries for `recob::*` etc.  ROOT keeps the first-loaded (LArSoft's)
  and warns.  The standalone chain is unaffected (gate above); if the wcls
  path ever complains, move the dictionary to a separate lazily-loaded lib.
