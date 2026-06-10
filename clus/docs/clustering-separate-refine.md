# Post-separation refinement: collinear tip recovery + two-band re-carve

Two knob-gated refinement steps in `ClusteringSeparate`
(`clus/src/clustering_separate.cxx`), applied to the "family" of clusters
produced by separating one original cluster.  Both default **OFF** (existing
configs are bit-identical); PDVD and PDHD enable both.

```jsonnet
cm.separate(use_ctpc=true, max_hull_points=100000,
            collinear_recover=true, band_recarve=true)
```

## Motivation (PDVD run 39324 evt 0, drift group 4567)

After the hull-cap fix (`max_hull_points=100000`, see
`clustering-separate-hull-cap.md`) the 29k-pt over-merged cluster — one
drift-aligned "vertical" track plus two wide crossing isochronous bands —
finally separates, but the carve was suboptimal in two structural ways:

1. **Stranded collinear tip.**  The vertical track's tip (190 points lying on
   its PCA axis, 3.2° collinear, <4 cm perpendicular) sits across a real
   6.07 cm imaging gap.  `Separate_1` carves its "path" cluster by 2D proximity
   to a graph shortest path, which stops at the gap; the final `Separate_2` is
   pure connectivity with a 5 cm relink, so the tip stayed chained to the
   leftover bin instead of the track.
2. **Arbitrary band fragmentation.**  The two bands cross at only ~13° in the
   y–z plane and overlap over a wide region.  `Separate_1` lumps everything not
   "saved" into ONE leftover cluster per recursion level, and only that bin is
   re-separated; the result was one 20k-pt mash (27 disconnected components at
   3 cm, containing chunks of BOTH bands) plus arbitrary slices of one band.
   Local directions inside a wide band vary by 40–60° (hatched isochronous
   strips), so no local angle cut can assign fragments per-band — only a
   global 2-line model can.

## Step A — `collinear_recover` (recover_collinear_tips)

For each family member `T` that is a long, thin, non-isochronous track
(`length > 50 cm`, `pca eval1/eval0 < 0.15`, axis at least 15° away from
perpendicular-to-drift — wide bands are mutually near-collinear and would steal
each other's far blobs), scan the other members' blobs for ones that continue
`T`'s PCA axis beyond its `get_main_axis_points()` endpoints:

- blob center within **4 cm** perpendicular of the axis,
- blob local direction (`vhough_transform`, 15 cm radius, sign-free) within
  **15°** of the axis,
- contiguous axial run: walking outward from the endpoint, consecutive axial
  gaps ≤ **15 cm**.

Qualifying blobs are moved into `T` (per-donor `Grouping::separate` with
`remove=false` + `take_children` + `destroy_child`, then the scope triple
`set_default_scope`/`set_scope_filter`/`set_scope_transform` re-applied to `T`
— `take_children` does not invalidate the facade PCA/length cache; the
`set_default_scope` round-trip performs the full `clear_cache()`).

### Interior fragment reclaim (`collinear_interior`, default OFF)

A separate knob (only effective when `collinear_recover` is also on; PDHD and
PDVD enable it).  At a track crossing, `Separate_1`'s carve can shed small
mid-track fragments of `T` — interior bites, not tips — which the later
proximity merge passes then attach to the WRONG track (PDHD 27409 evt 40900
with DNN-SP imaging: a ~24 cm chunk of one crossing track, holding the
user-flagged point, rode away with the other track's cluster).  The tip walk
cannot reach them (interior to `T`'s span) and `track_recarve` cannot either
(arms must be ≥20 %/60 cm).

For each qualifying track `T` (same gates as tip recovery), absorb whole
SHORT sibling fragments lying along `T`'s axis inside its span:

- fragment `length < 50 cm` (claimers are ≥ `TRACK_LEN_MIN` = 50 cm, so
  claimer/claimee sets are disjoint and mutual claiming is impossible),
- EVERY fragment blob center within **8 cm** perpendicular of `T`'s axis and
  axially interior to `T`'s span,
- fragment main PCA axis within **30°** of `T`'s axis when the fragment is
  ≥ 6 cm long (the evt-40900 chunk reads ~19°); below 6 cm the pca direction
  is noise and the geometry gates alone decide.

Fragment-level tests only — per-blob `vhough_transform` directions on a 20 cm
fragment are noise, and a per-blob variant either leaks crumbs or (with a
donor-spine gate) never fires because a carve fragment IS its own spine.
Donors already claimed by the tip walk are skipped.  The OFF-check (knob
absent) is byte-identical by construction and was verified empirically
(old-vs-new binary content-identical mabc zips on 27409 evt 40908).

## Step B — `band_recarve` (recarve_two_bands)

Gate: ≥2 family members that are band-like (PCA axis within **10°** of
perpendicular-to-drift, `length > 60 cm`), with a seed pair that **touches**
(closest points < **2 cm**) and has distinct y–z projected axes (angle in
**[15°, 45°]**; widest qualifying pair seeds the fit).  The pool grows from the
seed pair over band-like members by the same 2 cm touch (deterministic BFS in
family order).

Seed gates (revised with the FV-inset work, see `clustering-separate-fv.md`):
both seed members must be ribbon-like (rms width across the main axis ≥ 6 cm)
and the seed angle window is [15°, 45°] — thin crossing line-tracks that
survive the width gate pair up at ~11° and must never seed the x-blind
re-carve (PDHD 27409 evt 40900), while genuine two-band seeds sit at 25–39°.
A `track_recarve` step (same doc) handles crossing thin tracks instead.

Fit: k=2 line fit in the (y,z) plane — isochronous bands are extended
transverse to drift, so the carve is an x-independent partition.  **10 fixed
iterations** of assign-each-blob-to-nearer-line (blob `center_pos()`,
perpendicular distance to the infinite line, tie → line 0) / refit each line as
the npoints-weighted 2D PCA of its assigned blob centers (closed-form 2×2:
θ = ½·atan2(2c_yz, c_yy − c_zz)).

Aborts keeping the existing carve when: no valid seed pair; one side empties
during iteration; or the final npoints-weighted minority side is below **10%**
of the pool (degeneracy guard, checked *before* any mutation).  On success the
pooled members are merged into a fresh cluster and re-separated into exactly
two clusters by nearer-line blob assignment.

## Ordering and determinism

Step A runs before Step B (it pulls the vertical-track tip out of the leftover
bin before the re-carve pools the band-like members).  Family tracking captures
every cluster created while separating one original cluster — the four
`Separate_1` call sites plus the two final `Grouping::separate`
materializations (whose return values were previously discarded) — with
consumed intermediates value-erased.  All iteration is in vector (creation)
order; sorts carry (family index, blob index) tie-breaks; iteration counts are
fixed.  No pointer-keyed iteration.

## Result on the pinned event (39324 evt 0 group4567)

- Vertical track: whole (1451 → 1565 pts; the 114-pt tip recovered across the
  6.07 cm gap), pca eval ratio 0.012, separate from the bands.
- Band complex: exactly **2** clusters, one per physical band; the previous
  arbitrary {mash, slice, slice} carve (old clusters 51/52/57) is gone.  The
  user-flagged point pairs (−127.3,−178.1,153.5)/(−118.6,−178.8,148.8),
  (−154.4,−307.7,132.0)/(−153.8,−253.8,135.1) and
  (−155.1,−277.6,235.1)/(−153.8,−281.1,210.6) each land in one cluster.
- Drift group 0123 of the same event: one clean tip recovery (a long
  drift-aligned track adopts 42 collinear continuation points from a leftover
  380-pt fragment).  Nothing else changes.

## Regression

With both knobs OFF the outputs are byte-identical (verified on 39324 evt 0:
all three clustering JSONs md5-equal to the pre-change reference), so SBND and
any other existing config are unaffected.

Knobs-ON sweep (PDVD runs 39252/39253 evts 0–4 + 39324 evts 0–10; PDHD runs
027409 evts 0–8 + 027380 evts 0–7): the majority of events are identical;
every change is confined to a previously-separated family and is one of the
two intended flavors — a thin track reassembled from carve fragments
(tip recovery: e.g. 39324 evt 2, a 290 cm drift track rebuilt from 3 pieces,
final pca eval ratio 0.044) or a band complex consolidated to one cluster per
band (recarve: e.g. 39324 evt 4, 5 arbitrary pieces → 2 bands).  Determinism:
two knobs-ON runs of 39324 evt 8 (both steps firing on one family) are
byte-identical.  Full per-event detail in
`pdvd/docs/clustering-iso-overcluster-39324.md` (wcp-porting-img repo).
