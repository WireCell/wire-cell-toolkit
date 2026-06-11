# Post-separation refinement: collinear recovery/rejoin, repartition, band merge-back, re-carve, slab split

Knob-gated refinement steps in `ClusteringSeparate`
(`clus/src/clustering_separate.cxx`), applied to the "family" of clusters
produced by separating one original cluster.  All default **OFF** (existing
configs are bit-identical); PDVD and PDHD enable all of them.

```jsonnet
cm.separate(use_ctpc=true, max_hull_points=1000000,
            collinear_recover=true, collinear_interior=true,
            collinear_member_merge=true,
            track_repartition=true, band_merge_back=true, band_recarve=true,
            iso_slab_split=true)
```

Run order within one family: `collinear_recover` (+`collinear_interior`) →
`collinear_member_merge` → `track_repartition` → `band_merge_back` →
`band_recarve` → `track_recarve` (in `clustering-separate-fv.md`) →
`iso_slab_split`.

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

## Step A1b — `collinear_member_merge` (merge_collinear_members)

The carve can cut ONE straight cosmic into two (or more) long thin pieces
that nothing rejoins: `collinear_recover` claims *blobs* (a ≥50 cm sibling is
never an interior-reclaim donor, and the 4 cm tip-walk gate loses a
slightly-bent cosmic within a couple of metres), and the downstream
`connect1` does not always reach across (PDHD 27409 evt 40924: two pieces of
one cosmic — axes 6.5° apart, *touching* at 0.3 cm — stayed two clusters).
Worse, the leftover fragments around such breaks can mislead `connect1` into
merging two *different* cosmics (evt 40920's over-merge disappeared once the
fragments were rejoined into their proper tracks first).

Merge a pair of family members when ALL hold (iterated until no pair merges —
a track cut into three needs two rounds):

- both long (≥100 cm) and thin (pca eval1/eval0 ≤ 0.05), touching (≤5 cm),
- main-axis angle ≤ 10°, each centroid ≤ 15 cm perpendicular from the OTHER's
  pca main line,
- the union (npoints-weighted blob centers) is still one thin straight track:
  perp rms about its own 3D principal line ≤ **7 cm** and eigenvalue ratio
  ≤ 0.010.

Measured operating point: across the full PDHD+PDVD regression suite every
genuine one-cosmic rejoin reads union rms ≤ 5.7 cm (12 firings, angles up to
8.9°); the one pair that must NOT rejoin — two fork-adjacent (kinked) pieces
in 27409 evt 40904, whose rejoin caused the downstream `connect1` to fuse the
two full fork prongs — reads 8.1 cm.  The 7 cm gate splits the two
populations with margin on both sides.

## Step A2 — `track_repartition` (repartition_crossing_tracks)

At a crossing, the carve + `Separate_2`'s 5 cm relink can fuse a mid-track
segment of track B INTO track A's cluster (PDVD 39324 evt 339990 group0123: a
48 cm chunk of one crosser, 26 cm off the host's axis, embedded in the other's
cluster; the genuine crossing is at 49.6°, closest approach 3.7 cm).
`collinear_interior` cannot reach it — it moves whole sibling fragments, and
this chunk is fused into the host.

For each PAIR of family members that are both thin (`eval1/eval0 < 0.15`) and
long (>100 cm) and whose pca main lines genuinely cross (closest approach
<10 cm, interior [0.15, 0.85] to ≥1 member, within [−0.1, 1.1] of both): pool
the pair's blobs, seed `side` by current membership, run the same fixed-count
k=2 3D line assign/refit and validation gates as `track_recarve` (arm ≥20 %
npoints and ≥60 cm, rms residual ≤10 cm, crossing interior), then reassign
each blob to the nearer fitted line.

**Surgical no-op guard**: the would-be moves are computed first and the pair
is mutated ONLY if the moved blobs include a run with ≥15 cm axial extent
sitting >15 cm off its host's line.  Crossing-region dribble (close to both
lines by construction) never qualifies — on the five verified PDHD 27409
separation events the guard skips 5 of 8 candidate pairs and the three that
fire move genuine misplaced runs (73–110 cm extents); all five user-point
checks still pass.

## Step A3 — `band_merge_back` (merge_back_bands)

With separation firing on band-topology clusters (the FV-inset round),
`Separate_1`'s path carve can hatch ONE wide isochronous band into 2+ clusters
of interleaved alternating chunks (PDVD 39324 evts 339870/339930/339990/340010
— six user-flagged over-separations; verified visually and by the common
x-slab).  The pieces are mutually near-collinear so no recarve seed exists,
and nothing re-assembles them.

Pool the touching (<2 cm, transitive) band-like family members, in two tiers:

- **wide members** (axis within 10° of ⊥-drift, len >60 cm, rms width ≥6 cm —
  the recarve band test + its seed width gate) can ANCHOR a band;
- **thin debris** (width 2–6 cm, len 60–100 cm, same axis gate) may join a
  pool and merge into a nearby anchor but can never anchor or merge on its
  own.  The width gate is what keeps thin crossing-track pieces from forming
  bands, and the 100 cm debris length cap is what keeps long thin iso-aligned
  TRACKS out (PDHD 27409 evt 40908: 318/503 cm tracks at 9.5° off-⊥, width
  2–3.5 cm, were briefly eaten without the cap).

For each pool, fit ONE line in (y,z) to the npoints-weighted blob centers and
characterize each member by its signed transverse offset mean μ and rms σ:

- **crossing guard**: if the union residual GROWS toward the axial ends
  (outer/inner rms > 1.75 over the t-quartiles) the pool is skipped — two
  genuine bands crossing as an X read ~2.0, the worst hatched single band
  1.48.
- **anchor grouping**: anchors single-link iff |Δμ| < 85 cm AND
  |Δμ| < 1.7·(σᵢ+σⱼ) AND |Δx̄| < 20 cm.  Hatched-piece anchor gaps measure
  ≤67 cm; the two genuine side-by-side parallel bands of evt 339850 (which
  touch at 0.8 cm and must stay TWO clusters) sit 130 cm / 2.27·(σᵢ+σⱼ)
  apart.  Plain single-linkage over all pieces is the wrong shape:
  overlap-region debris chains the two parallel bands together (this
  re-merged 339850 in the first iteration).  The **x-slab gate** (|Δx̄|,
  npoints-weighted blob-center x means) exists because the y-z fit is
  x-blind: pieces of ONE hatched band share its ~12 cm drift-time slab, while
  two bands whose y-z projections overlap can sit 50 cm apart in x (39324
  evt 339890) and must never link.
- **debris assignment**: each non-anchor joins the nearest *same-slab* anchor
  by |Δμ| (absolute gate only — debris σ is unreliable); too-far debris stays
  its own cluster.
- **union-rms cap**: before mutating, each multi-member group's own y-z line
  is refit and its npoints-weighted transverse rms must stay band-sized
  (≤ 50 cm).  A genuine re-assembled band reads 17–43 cm across the suite;
  two distinct same-slab complexes that the μ-gates wrongly grouped read
  57–64 cm (27409 evt 40900 — a regression caught in round 5 — and 39324
  evt 339890).  Rejected groups print `rejected group of N (union rms ...)`.

Each surviving multi-member group merges into its lowest-family-index member
(`take_children` + the scope round-trip cache clear).

**Band interior steal** (same knob): a chunk of a band can stay FUSED inside a
non-band sibling — e.g. a crossing drift-angled track (evt 339930: a ~70 cm
band chunk rode in a 248 cm track member whose overall PCA is 48° off-⊥, so
the member-level merge can never take it without dragging the track along).
For each band-like member B and each long (>100 cm) non-band sibling H: take
H's blobs sitting >15 cm off H's own main line, seed with those touching B
(<5 cm), grow within the off-line set (15 cm radius), and move the grown run
into B when its spatial extent is ≥15 cm.  Blobs on H's line — including the
genuine crossing region — never move.  **x-gate** (round 5): the selected run
must live in B's x-slab — run x-extent ≤ 20 cm AND run x-mean within 20 cm of
B's x-mean.  A genuinely fused band chunk is isochronous (x-narrow); a
drift-track run is not (39324 evt 339990: a 193 cm run spanning 66 cm in x
was a real drift track, wrongly stolen into the band — both of that event's
steals are now rejected and the track stays whole).

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

## Step C — `iso_slab_split` (split_iso_slabs)

An isochronous band lives at ONE drift time — a narrow, dense x-slab — while
a drift-direction track spans x and can touch several bands.  The final
`Separate_2` (pure 5 cm connectivity) then chains band–track–band into one
cluster FOREVER: every y-z-projected mechanism above is x-blind, and a drift
track nearly vanishes in the y-z projection.  PDVD 39324 evt 339890 group4567
(the round-4 "known residual"): a 44.7k-pt member holding TWO bands at
x-slabs 50 cm apart (24k pts at x∈[59,67], 14k at x∈[115,123]) plus three
drift tracks; evt 339990: one band plus a drift track fused by the carve
itself (not by the steal — verified by rejecting the steal and watching the
fusion persist).

Partition each family member (≥100 cm, ≥20 blobs, ≥8 x-bins of extent) by
x-slab occupancy:

1. **slabs**: weighted blob-center x histogram (4 cm bins); a slab is a
   contiguous run of bins each holding > 5× the uniform density, total width
   ≤ 25 cm, ≥ 20% of member npoints, spanning ≥ 60 cm in (y,z) (band-sized).
2. **tracks**: out-of-slab blobs form connected components (10 cm);
   components with 3D extent ≥ 30 cm AND x-extent ≥ 20 cm (the drift
   signature — band shoulders just outside the slab core read ≤ 12 cm in x)
   seed tracks; collinear seeds (angle ≤ 15°, mutual line offset ≤ 10 cm)
   join — a track pierced by a slab reappears beyond it.
3. **assignment**: track-component blobs keep their track; any other blob
   within 6 cm of a track line joins that track (the track's continuation
   THROUGH a band); remaining in-slab blobs form one band cluster per slab;
   leftover debris joins the nearest slab in x.

Fires only when ≥1 valid slab AND ≥1 track exist, so a pure band, a pure
track, an isochronous *track* (a slab with no crossers) and a band plus
dribble are all structural no-ops — verified across the full regression
suite, where it fires exactly five times: 39324 evt 339890 (4 tracks +
2 bands — all five user-flagged structures land in distinct clusters), evt
339990 (1+1), evt 340010 (1+1), 27409 evt 40900 (1+1, peeling a drift track
off the band complex) and 27380 evt 1 (1+1).

One systematic side effect, by design of the downstream chain: a structure
freed from a big track cluster is now judged on its own by the group-stage
deghost/isolated passes, where previously it hid inside the big cluster.
Display point counts can therefore change beyond the immediate split
(largest observed: 27380 evt 1 drops a 5.5k-pt single-drift-time sheet
spanning the full y-z face — the classic ghost signature; 27409 evt 40920's
event gains back 1.2k points that used to be dropped).  On the user-scanned
runs the net deltas are ≤0.2%.

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
