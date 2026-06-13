# Stage 3 — Validation of the WCT light reconstruction vs LArSoft

Tool: `flash/test/check_light_reco.py <workdir>...` — compares the stage-2
products against the stage-1-converted LArSoft reference in the same
`pdhd/work/<RUN6>_<EVT>/`, writes plots to `<workdir>/light_validation/`
(`decon_compare.png`, `hit_compare.png`, `flash_compare.png`).

Run on one event of each example run (2026-06-12, with the production
run28368 v1 SPE templates + run27950 noise spectra, see
`stage2-reconstruction.md`):

| run/evt | decon corr (med) | decon peak ratio | peak shift | hit exact-pulse match (of ref) | hit surplus (of wct) | flash matches (wct/ref) | flash dt med | flash PE ratio* |
|---|---|---|---|---|---|---|---|---|
| 27305/150 | 1.0000 | 1.000 | 0 | 417/419 | 87/504 | 16 (46/28) | +4.1 µs | 1.96 |
| 27980/8 | 1.0000 | 1.000 | 0 | 1081/3828 | 2/1083 | 60 (83/265) | +3.2 µs | 1.00 |
| 28084/74408 | 1.0000 | 1.000 | 0 | 1755/9685 | 4/1759 | 80 (107/444) | +2.9 µs | 0.92 |
| 29107/983 | 1.0000 | 1.000 | 0 | 868/8089 | 4/872 | 69 (104/396) | +3.2 µs | 0.95 |

\* PE ratio on channels where both flashes have light, median over matched
pairs (matching: nearest in time within 10 µs, largest-first).

## Reading of the numbers

- **Deconvolution is exact**: with the production v1 per-channel SPE
  templates and the run27950 noise spectra, every snippet matches the
  in-file reference at correlation 1.0000, peak ratio 1.000, shift 0
  (residuals at the 1e-6 float level).
- **Hits are exact pulse-by-pulse**.  Hits are matched on the exact pulse
  parameters (channel, width, amplitude, area), because the reference
  PeakTimeAbs is unusable: the production stamps
  `TimeStamp + TickPeriod·t_max` where the DAPHNE decoder `TimeStamp()` is
  in DTS *ticks* while the offset is in *µs*, so every sub-pulse of a
  snippet collapses onto the snippet head (verified: the collapsed stamps
  reproduce exactly under the 16-tick double quantization).  Our hit times
  are the per-pulse times and are *more* correct than the reference.
  - On busy runs 99.5–99.8% of our hits have an exactly-equal reference
    pulse (2–4 hits/event differ); the unmatched *reference* hits are the
    waveforms the 400-snippet/event dump cap never gave us.
  - On 27305/150 the surplus is concentrated on channels 41/43/44/45/47:
    LArSoft drops these with "unrecognized channel number" (geometry
    `IsValidOpChannel` fails); we keep them — they are real OpDets in
    `opdet_geo` and contribute light to matching.
- **Flashes**: common-channel flash PE agrees at the 5–10% level on the
  busy runs and y/z centroids track.  Counts/times are not directly
  comparable: the reference flashes were assembled from *all* readout
  channels (no dump cap) using the collapsed hit times above, and the
  +3–4 µs PE-weighted offset follows from that time structure.

## Q-L matching contract (smoke check)

Both `opflash_pdhd.tar.gz` (converted LArSoft flashes) and
`opflash_pdhd-wct.tar.gz` (our reconstruction):

- round-trip through `TensorFileSource` → `TensorFileSink` in a wire-cell
  job (the same source step that feeds `FlashTensorToOpticalPCs` for SBND);
- satisfy the `FlashTensorToOpticalPCs` input contract: tensor 0 is 2-D f8
  with `ncol == nchan+1 == 161` (checked after the WCT round-trip).

With `FlashTensorToOpticalPCs {nchan: 160}` + `QLMatching {nchan: 160}` (and
a future `semi-analytical-pdhd.json` built from `pdhd-opdet-geom.json`), the
PDHD charge-light matching can consume either archive unmodified.

## Follow-ups recorded

1. ~~Obtain the run28368 per-channel v1 SPE templates~~ — DONE: retrieved
   from the DUNE StashCache cvmfs (via dunegpvm), together with the
   run27950 noise spectra; deconvolution now exact.
2. The trigger-type heuristic (nearest 250 µs) and the 400-snippet dump cap
   are properties of the temporary exchange format; both disappear with the
   future light-data format.
3. PDHD Q-L matching itself (SemiAnalyticalModel parameters or photon-library
   route, `pdhd/docs/photon-detector-chain.md` §7) is the next subsystem.
