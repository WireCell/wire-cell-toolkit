# Stage 1 — ROOT → WCT conversion of PDHD light data

Status: complete, validated on one event of each of the four example runs
(27305/150, 27980/8, 28084/74408, 29107/983).

## Components (package `root/`)

### `PDHDOpFlashSource` (`ITensorSetSource`)

Reads one event's LArSoft-reconstructed optical products and emits the opflash
tensor set of `design.md` §3.4:

- flash enumeration + y/z summaries from `opflashana/PerFlashTree`;
- per-OpDet PE from `flashopdet/flash_opdet` (OpChannel == OpDet, 0..159);
- OpHits from `opflashana/PerOpHitTree` (tensor 2; `flash_id = -1`: the dump's
  `FlashHitMatchTree.HitID` is not a usable join key);
- trigger anchor from `trigoff/trigger_offset`.

Config: `filename`, `run`, `event` (required), `subrun` (-1 = any), `nchan`
(160), `nominal_offset_us` (250.0), `include_summary`/`include_ophits` (true).

### `PDHDOpWaveformSource` (`IFrameSource`)

Reads the `decoana` TH1D snippets into one `IFrame`: one trace per snippet
(channel = OpChannel, 1024 ticks), tagged trace sets `"raw"` (ADC, faithful)
and `"deconv"` (LArSoft SPE-normalized).  `tick = 16 ns`.

Config: `filename`, `run`, `event`, `subrun`, `nominal_offset_us` (250.0),
`frame_tag` ("pdhd_light"), `min_peak` (2.0), `t0_search_window_us` (5000).

## Time anchoring — what was learned from the data

- `PerFlashTree.FlashTime`, `PerOpHitTree.PeakTimeAbs`, `flash_opdet.flash_time`
  are **absolute DTS timestamps in 16 ns ticks stored as doubles**.  At
  ~1.07e17 ticks the double ULP is 16 ticks, so all these times are quantized
  to **±8 ticks (±128 ns)** in the dump itself.  Converted times inherit this.
- Conversion to trigger-relative WCT ns: `t = (T_dts − tc_time_candidate)·16 ns`
  with the trigger candidate chosen per event as the one whose `offset_us` is
  closest to `nominal_offset_us` (250 µs); choice recorded in the tensor-set
  metadata (`tc_type`, `tc_time_dts`, `offset_us`).
- **decoana histogram x-axis**: the bin *low edge* numerically encodes the
  snippet start in 16 ns ticks **relative to the event's earliest saved
  snippet** (`t_first`); the bin *width* is written as 0.016 (µs/tick) — only
  the offset is meaningful.  `t_first` itself is not recorded anywhere, and is
  event-specific (observed from −195k to +2k ticks relative to
  `rd_timestamp`, i.e. snippets begin up to ~3 ms before the TPC readout
  start).
- `t_first` recovery: every clear deconv peak (max ≥ `min_peak`) and every
  same-channel OpHit vote `t_first = PeakTimeAbs − snippet_start − peak_bin`;
  votes within `t0_search_window_us` of `rd_timestamp` are histogrammed, the
  ±2-tick-smoothed mode wins.  The true peak is sharp (50–160 votes vs ~30
  background).  Frame `time() = (t_first − tc)·16 ns`, trace `tbin` = snippet
  start in ticks from `t_first`.

## Config + job

- `cfg/pgrapher/experiment/pdhd/flash.jsonnet` — builders `opflash_source`,
  `opwaveform_source`, `opflash_sink` (TensorFileSink, prefix `opflash_`),
  `waveform_sink` (FrameFileSink).
- `cfg/pgrapher/experiment/pdhd/pdhd-opdet-geom.json` — 160 OpDet positions
  (mm), extracted by `flash/test/extract_pdhd_light_maps.py`.
- `cfg/pgrapher/experiment/pdhd/pdhd-opch-map.json` — hardware DAPHNE channel →
  OpDet (informational; the reconstruction basis is offline OpChannel == OpDet).
- Job (wcp-porting-img repo): `pdhd/wct-light-convert.jsonnet`, run via
  `pdhd/run_light_evt.sh <run> <evt>`, writing to `pdhd/work/<RUN6>_<EVT>/`:
  - `opflash_pdhd.tar.gz` — opflash tensor set (QLMatching-ready archive);
  - `light-frames.tar.bz2` — `frame_{raw,deconv}_<evt>.npy` (+ channels,
    tickinfo) dense per-tag arrays.

## Validation — `flash/test/check_light_convert.py`

`python3 check_light_convert.py <rootfile> <run> <event> <workdir>` checks:

1. tensor-set metadata (run/event, chosen trigger);
2. opflash matrix shape `[nflash, 161]`, flash times == `(FlashTime − tc)·16 ns`,
   PE row sums == `flash_total_pe` == `PerFlashTree.TotalPE`;
3. summary y/z centers/widths (cm→mm);
4. ophit count and trigger-relative peak times;
5. frame `tick == 16 ns` and the time anchoring: ≥80% of clear saved deconv
   peaks coincide with an OpHit (within 16 ticks, or within the hit's width
   for wide merged hits).

Results: all checks pass on the four runs.  Notes:

- the dump caps `decoana` at **400 snippets/event**, so only a subset of hit
  waveforms is present (the check runs snippet→hit, not hit→snippet);
- busy events show a small (<20%) population of secondary deconv peaks
  clustered ~±5 µs from the nearest hit — explained in stage 3: the
  production's OpHit PeakTimeAbs collapses every sub-pulse of a snippet onto
  the snippet head (DTS-tick TimeStamp + µs offset units mix-up, see
  `validation.md`), so secondary pulses genuinely sit away from their stamped
  hit time — structured, not an alignment failure (the t_first vote uses the
  dominant first-pulse population, which is unaffected);
- run 29107 carries three trigger types (14/29/31); the nearest-to-250 µs
  selection picks tc_type 14 there and the anchoring checks pass.
