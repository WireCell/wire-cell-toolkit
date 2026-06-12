# Stage 3 — Validation of the WCT light reconstruction vs LArSoft

Tool: `flash/test/check_light_reco.py <workdir>...` — compares the stage-2
products against the stage-1-converted LArSoft reference in the same
`pdhd/work/<RUN6>_<EVT>/`, writes plots to `<workdir>/light_validation/`
(`decon_compare.png`, `hit_compare.png`, `flash_compare.png`).

Run on one event of each example run (2026-06-12, toolkit @ stage-2 commit):

| run/evt | decon corr (med) | decon peak ratio | peak shift | hit eff | hit surplus | flash matches (wct/ref) | flash dt med | flash PE ratio* |
|---|---|---|---|---|---|---|---|---|
| 27305/150 | 0.969 | 0.984 | +2 ticks | 401/419 | 736/827 | 16 (55/28) | +4.1 µs | 2.19 |
| 27980/8 | 0.957 | 0.932 | +1 | 835/1113 | 1624/1907 | 59 (89/265) | +3.3 µs | 1.06 |
| 28084/74408 | 0.978 | 1.015 | +1 | 1408/1769 | 2782/3205 | 90 (132/444) | +2.7 µs | 0.91 |
| 29107/983 | 0.971 | 1.031 | +1 | 750/919 | 1539/1744 | 70 (114/396) | +3.2 µs | 0.96 |

\* PE ratio on channels where both flashes have light, median over matched
pairs (matching: nearest in time within 10 µs, largest-first).

## Reading of the numbers

- **Deconvolution**: shift-aligned correlation 0.96–0.98 and peak ratios
  within a few % everywhere.  The +1–2 tick shift and the residual shape
  difference are the documented SPE-template provenance issue
  (`stage2-reconstruction.md`): the reference was deconvolved with the
  per-channel run28368 v1 templates, which are not retrievable from the
  accessible cvmfs products; we use the official 2024 NP04 FBK/HPK
  templates.  The C++ implementation itself reproduces an independent numpy
  replication of the LArSoft module to ~1e-6.
- **Hits**: 75–96% of reference hits are found at the same time.  The large
  surplus (small hits, median ~1 PE) and the inflated matched-hit PE (long
  pulse trains kept alive by the template-mismatch baseline absorb area that
  LArSoft splits into several hits) are both downstream symptoms of the same
  template difference; on the *reference* deconv input the ported sliding
  window finds 504 pulses vs the 419 PerOpHitTree entries of evt 27305/150.
- **Flashes**: the matching-relevant quantities are reasonable: on the busy
  runs the common-channel flash PE agrees at the 5–10% level and y/z
  centroids track (see plots).  The +3–4 µs PE-weighted time offset follows
  from the extra late-light hits.  Absolute flash counts are not comparable:
  the reference flashes use *all* readout channels while the waveform dump is
  capped at 400 snippets/event, so our chain only sees a subset of the light.

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

1. Obtain the run28368 per-channel v1 SPE templates (dune_pardata release or
   collaborator) and regenerate `pdhd-spe-templates.json` — expected to
   collapse the decon shift/baseline and the hit-level surplus/PE inflation.
2. The trigger-type heuristic (nearest 250 µs) and the 400-snippet dump cap
   are properties of the temporary exchange format; both disappear with the
   future light-data format.
3. PDHD Q-L matching itself (SemiAnalyticalModel parameters or photon-library
   route, `pdhd/docs/photon-detector-chain.md` §7) is the next subsystem.
