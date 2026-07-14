# data evt 59415 — a cathode crosser missed because its flash is split across APAs

**Status: NOT an algorithm bug — a light-data problem.** Diagnosis of the
cross-TPC case the lan-reco2 data hand-scan flagged in event 59415 (file 1 of the
3-file lan-reco2-Run1 BNB set). The hand-scan note was:

> Two clusters did not get xTPC: (x, y, z) = (31.8, −57.0, 451.1) vs.
> (x, y, z) = (−2.5, −81.2, 450.1)

Conclusion: the two halves are a genuine cathode crosser, but the **flash that
should tie them together is reconstructed as two separate flash objects 4.3 µs
apart with a 58× PE asymmetry**. Cross-TPC pairing requires the two halves' flash
times to agree to within 80 ns, so they are never co-evaluated. This is a flash
(light) reconstruction defect, not a Q/L matching shortcoming — the matching
behaves correctly given the flashes it is handed.

All numbers from `sbnd_xin/work/ql_evt59415/{calib-evt59415.json,
wct_ql_evt59415.log}` (drift speed 0.1563 cm/µs).

## The crosser

| half | cluster | raw x (cm) | npts |
|------|---------|------------|------|
| TPC0 | APA0 ident 9 (uid 9)        | [−234.1, −197.8] (beyond −200 anode) | 862 |
| TPC1 | APA1 ident 7 (uid 1000007)  | [199.1, 232.5]   (beyond +200 anode) | 311 |

Both sit outside their TPC box in raw coordinates (the signature of a cathode
crosser before T0 correction), at the same z ≈ 451, with Hough axes parallel to
**2.0°**. Geometrically an easy cross-TPC pair.

## What the matcher did

| half | selected flash | t0 (µs) | strength | pred light / meas PE |
|------|----------------|---------|----------|----------------------|
| APA0 clus 9        | flash 12      | −1257.9 | 0.163 | 3369 / 213 |
| APA1 clus 1000007  | flash 1000009 | −1120   | 0.973 |  522 / 405 |

The two halves landed on flashes **138 µs apart** → never seen as the same event
→ no xTPC flag, no pairing.

## Root cause — the flash is split across the two APAs

At the matched t0 (≈ −1260 µs) the two cathode ends meet:

```
t = −1258 µs (off 196.6 cm):  TPC0 end @ −1.2 cm   TPC1 end @ +2.5 cm   gap 3.6 cm
t = −1262 µs (off 197.2 cm):  TPC0 end @ −0.5 cm   TPC1 end @ +1.8 cm   gap 2.4 cm
```

That is a trivial scenario-1 match (closest approach 2.4–3.6 cm ≪ dmax = 5 cm).
But the crosser's flash is reconstructed as **two separate objects**:

- APA0: flash **12** at **−1257.9 µs**, **213 PE**
- APA1: flash **1000003** at **−1262.2 µs**, **12395 PE**

→ 4.3 µs apart, 58× PE asymmetry, for what physics says is a single flash (one
track, one t0; a01 = 2°, y/z aligned).

`cull_cross_tpc` only pairs two candidate bundles when their flash times agree to
within `flash_group_window = 80 ns` (`QLMatching.cxx:2095`). 4.3 µs is ~54× that
window, so cluster 9 (on flash 12) and cluster 1000007 (on flash 1000003) — both
valid `window_truncated` candidates, both contained at that t0 — are **never
co-evaluated**. In the log cluster 9 is paired only at t = −365, −1120, −507 µs
(flashes that *do* have a sub-µs cross-APA partner), and never at −1258/−1262.

### Downstream symptom

The "TPC1 matched the wrong flash" seen in the scan is the consequence, not the
cause:

1. flash time split (4.3 µs) → no coincident pairing →
2. xtpc flag never sets → `cull_inconsistent` never forces cluster 1000007 onto
   flash 1000003 →
3. the truncated TPC1 stub free-matches the faint flash 1000009 (−1120, 405 PE),
   whose light (522 predicted) it fits well.

### Validation

At t = −1120 the log reports closest approach `d = 46.78 cm`, which equals the
pure drift-x gap at that (wrong) t0 — 47 cm. The two cathode ends therefore
already coincide in y,z; at the correct t0 the 3-D distance collapses to the ~3 cm
x-gap. So this is unambiguously a scenario-1 crosser that the pairing step simply
never got to look at.

### Red herring — the −1351 µs flash pair

The faint flashes gid 9 (APA0, 408 PE) / gid 1000007 (APA1, 475 PE) at −1351 µs
*are* a tight 80 ns cross-APA pair, but correcting to that t0 over-shoots both
halves ~13 cm past the cathode and neither cluster has a candidate bundle there
(`require_containment` drops the over-the-cathode protrusion). Not the operative
flash.

## Why this is a data problem, not an algorithm problem

A cathode-crossing track has a single t0. Both PMT arrays should report the same
flash time; the per-TPC PE asymmetry is expected (each array sees mostly its own
TPC's light, and most of this track's length is in TPC1), but a **4.3 µs time
disagreement is not physical** — it points to a corrupted / mis-timed light
waveform or flash-finder output for this event's faint TPC0 flash. The Q/L matcher
is doing the right thing with the flashes it is given:

- both crosser-half bundles are built, contained, and retained;
- they are correctly *not* paired, because the inputs say their flashes are 4.3 µs
  apart;
- widening the 80 ns coincidence window to mask this would be wrong — the window is
  reused elsewhere (`:1821`, `:1862`) and a 4 µs tolerance would manufacture false
  cross-TPC coincidences across genuinely different flashes.

No code change is warranted. The fix, if any, belongs upstream in flash
reconstruction / light data quality. Recorded here so the next hand-scan does not
re-litigate it as a matcher bug.

## Contrast with evt12 (the other xTPC case)

MC evt12 (`qlmatching-evt12-xtpc-diagnosis.md`) was a genuine *algorithm* ordering
bug — the crosser bundle was built but culled before pairing, and was fixed by
reordering flag-before-cull + xtpc-priority steering. evt59415 is the opposite: the
algorithm is correct; the **light input is split across APAs**, so the pair is
never pairable in the first place.
