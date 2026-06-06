# SBND two-track upstream-boundary tag (`clustering_separate`)

## Problem

Data event 139220, TPC1 (APA1): two beam-inclined cosmics that imaging merges into a
single cluster (`cluster_id=4`, 12127 pts, 505 cm). `clustering_separate` is the only
pass that can cut them apart (it runs per-APA only; the all-APA pipeline has no
`separate`), and it declined to:

- `JudgeSeparateDec_2` → **false**
- `JudgeSeparateDec_1` → true, but the secondary PCA/angle gates fail
  (`angle_beam=49.3°`, `pca_ratio1=0.17`, `flag_top=0`) → `flag_proceed=0`.

The merged cluster's PCA looks like one beam-inclined cosmic, so relaxing the PCA/angle
gates would over-split genuine single cosmics. The discriminating feature is **topology,
not PCA**: the cluster has three charge-dense endpoints sitting on detector boundaries,
**two of them on the upstream (z-min) wall**. A single track touches a planar face at
most once, so two well-separated dense endpoints on one face imply two tracks.

### Why `JudgeSeparateDec_2` misses it

`JudgeSeparateDec_2` counts a hull extreme as "outside" only past the **inner-FV minus
margin** threshold. SBND's FV is inset 1 cm plus a 2.5–3 cm margin, so the upstream
"outside" test is `z < FV_zmin − zmin_margin = 0.85 − 3 = −2.15 cm`. The two real cosmic
endpoints land right at the wire edge (measured: z ≈ 0.6 and z ≈ 0.9–2.85 cm), i.e.
*inside* that threshold. Only one of them dips below `FV_zmin`, so Dec_2's per-axis
expansion captures a single upstream extreme — never the two it needs for its
`independent_surfaces.size() > 1` criterion.

## The tag — `JudgeSeparateDec_SBND_boundary`

SBND-only, default **OFF** (`sbnd_boundary_tag`), purely additive: it can only flip
`flag_proceed` false→true, never alters an existing accept. On a trigger it seeds the
detected tips into `independent_points` so `Separate_1`'s endpoint selection picks the
real prongs; the existing `flag_dec1` recursive `Separate_1`/`Separate_2` path does the
actual cut.

Trigger (all required), measured/pinned against evt 139220:

| step | cut | 139220 |
|------|-----|--------|
| long cosmic | `cluster_length > 250 cm` | 505 cm |
| upstream tips | ≥2 dense (`nnearby>75` in 15 cm) hull pts in band `z < FV_zmin + 10 cm`, grouped at 25 cm | 2 groups |
| separation | the two tips ≥ 40 cm apart | 64.6 cm |
| gap guard | midpoint of the two tips is ≥ 10 cm from cluster charge (a single grazing track keeps it ~0–3 cm) | 14.3 cm |
| 3rd exit | ≥1 dense hull pt within 5 cm of a non-upstream face | bottom y≈−199 |

The scope (upstream-only double-exit + a required 3rd exit elsewhere) is the
lowest-false-split form; it was chosen over a general "any face" / "two-tips-alone" rule.

## Verification

- evt 139220: the 12127-pt overcluster splits into two long tracks (7001 pt / 365 cm and
  5060 pt / 438 cm); the two hand-label points (81.9,10.6,72.3) and (102.8,−23.2,69.0)
  land in separate clusters. All other clusters in the event are unchanged.
- **OFF byte-identical**: with `sbnd_boundary_tag=false` the apa0/apa1/all-APA clustering
  is byte-identical to the pre-change build.
- **Regression**: a 32-event data+MC spread (incl. the cathode-connect accept events),
  tag ON vs OFF, is per-APA byte-identical on all but evt 139220. One all-APA-only
  discrepancy on evt 138670 was traced to the pre-existing intermittent `clus_all_apa`
  non-determinism (the tag does not fire on 138670; a fresh tag-ON run matches tag-OFF),
  i.e. unrelated to this change.

## Files

- `clus/src/clustering_separate.cxx` — `JudgeSeparateDec_SBND_boundary` (static helper) +
  call in the secondary gate + `sbnd_boundary_tag` config plumbing.
- `cfg/pgrapher/common/clus.jsonnet` — `separate()` factory gains `sbnd_boundary_tag`
  (conditional key; omitted when false → existing configs bit-identical).
- `cfg/pgrapher/experiment/sbnd/clus.jsonnet` — per-APA `cm.separate(...,
  sbnd_boundary_tag=true)`.
