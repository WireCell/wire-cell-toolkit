# Group-stage connect1 + deghost (PDHD/PDVD stage 3)

2026-06-10.  In the MicroBooNE single-TPC chain `separate` is followed by
`connect1` (reconnect "dashed-line" track fragments using the induction-plane
angles) and `deghost` (remove ghosts by 2D-projection overlap).  When the
PDHD/PDVD chains moved `separate` to the per-drift-group stage (stage 3), no
equivalent of those two follow-up passes could run there: `connect1` raised on
any multi-wpid grouping and `deghost` raised on any multi-APA grouping.  Both
are now generalized to the drift-group scope and inserted right after
`separate` in both detectors:

```jsonnet
// pdhd stage 3                      // pdvd stage 3
cm.separate(...),                    cm.separate(...),
cm.connect1(),                       cm.connect1(allow_mixed_faces=true),
cm.deghost(empty_view_unique=true),  cm.deghost(allow_mixed_faces=true, empty_view_unique=true),
cm.examine_x_boundary(),             cm.examine_x_boundary(allow_mixed_faces=true),
```

Motivating case: PDVD run 39324 event 339850 — a drift-direction track split
into two collinear fragments ~96 cm apart in x at the same (y,z)
((152.1, −262.2, 75.9) and (248.1, −260.8, 69.1) cm).  Stage-1 per-face
connect1 ran before the fragments were assembled by the group-stage merges;
the group-stage connect1's prolonged branch (3× = 150 cm skeleton extension)
now merges them (verified: both points in one cluster of `mabc-group4567.zip`).

## Validation of the grouping

Both passes accept a multi-wpid grouping only when the LIVE wpids form one
drift volume, via the shared `validate_drift_group()`
(`clus/src/ClusteringFuncs.cxx`), copying the `examine_x_boundary` semantics:
identical `FV_xmin/FV_xmax(+margins)` metadata across wpids, same face unless
`allow_mixed_faces` (PDVD: an anode's two faces are the y-halves of one CRP
and share one drift volume).  GOTCHA: validate the **live** `wpids()`, not
`dv_wpids()` — a face-restricted drift-group DetectorVolumes (PDHD
`detector_volumes(anodes, face)`) still reports both faces of its anodes, and
the opposite faces legitimately carry different FV_x (this raised on PDHD
group13 in the first implementation round).

## connect1 design decisions

- **Per-volume geometry**: `wpid_params` (drift dir + U/V/W angles per wpid)
  and `af_dead_{u,v,w}_index[apa][face]` dead maps are built for every wpid,
  exactly the deghost pattern.
- **Prolonged classification is per cluster**: the direction's transverse
  alignment with U/V/W wire directions (7.5°) is evaluated against each volume
  the cluster's blobs occupy (`wpids_blob_set()`, cached) and OR-ed.  A
  direction prolonged in any occupied view suffers the same 2D ambiguity the
  3× skeleton extension exists to handle.  Single-wpid groupings degenerate to
  the legacy computation bit-for-bit (`eval_prol` keeps the original U→V→W
  early-exit operation order).  The parallel test uses only the drift axis and
  is unchanged.
- **Per-point routing**: the overlap-counting loop fetches
  `cluster->wire_plane_id(j)` (cached array read) and routes both the
  dead-wire lookup and the 2D skeleton query through that point's
  (face, apa).
- **Extrapolation bucketing**: `make_points_linear_extrapolation` gains a
  trailing `seed_wpid` (the wpid of the extreme point being extrapolated
  from; computed only for multi-wpid groupings, 2 kd-knn per cluster).  With
  a multi-volume `wpid_params`, each synthetic point is re-bucketed into the
  volume containing it (`dv->contained_by`, grid-accelerated box test) so a
  ray crossing an APA (PDHD z-neighbors) or face (PDVD y-halves) boundary
  lands in the 2D KD trees that the target cluster's points actually query;
  out-of-volume points keep the seed.  The single-volume path is the legacy
  code verbatim.
- `vhough_transform` still pools all volumes in 3D — accepted, direction
  estimation is face-insensitive in practice (same reliance as
  `extend`/`regular` at this scope).

## deghost design decisions

deghost was already per-(apa,face) internally (per-wpid angle map, per-point
`test_wpid` routing); the guard relaxation above is the only structural
change.  One semantic addition was REQUIRED at group scope, gated by the
default-OFF knob `empty_view_unique`:

- An empty per-(plane,face,apa) 2D index returns distance −1 from
  `get_closest_2d_point_info`, which satisfies the `<= dis_cut/3` (and
  `<= dis_cut*2`) comparisons and tallies a bogus `nullptr` "overlap".  At
  group scope the clouds are seeded volume-by-volume in cluster-length order,
  so the first (longest) cluster of every not-yet-seeded volume would count
  zero unique points and be **wrongly destroyed**.  (Within one APA this is
  masked by the wrapped-wire cross-face fills; the rare existing case is the
  `find_max_cluster` nullptr-skip.)
- With `empty_view_unique=true` such a point counts as unique evidence
  instead: the longest cluster per volume self-seeds the live/skeleton clouds
  and the volume thereafter behaves exactly like the single-APA case.

Ghost adjudication itself is per-volume by construction (each point tests in
its own (plane,face,apa) KD tree), and the merge fallback uses 3D
closest-points — both volume-correct for clusters spanning APAs.

## Knobs

| knob | C++ default | PDHD stage 3 | PDVD stage 3 | meaning |
|---|---|---|---|---|
| `allow_mixed_faces` (connect1, deghost) | false | — | true | waive same-face in drift-group validation |
| `empty_view_unique` (deghost) | false | true | true | empty per-volume view counts as unique (required at group scope) |

Both knobs default OFF == byte-identical; existing stage-1 connect1 and
stage-2 deghost instances are untouched.

## Verification (2026-06-10)

- **Byte-identity** (new C++, old configs; A/B against a baseline binary built
  from HEAD with the changes stashed): SBND mc 10 evts + data 10 evts, PDHD
  027409 evt 0, PDVD 39324 evts 0+4 — every mabc zip content-identical.
- **Debug case**: PDVD 39324/339850 fragments merge (above).
- **Deghost sanity** (new config vs baseline): PDVD 39324 evt 0 group4567
  48→48 clusters (−474 ghost points), evt 4 group4567 51→50 (−725 points),
  group0123 unchanged/-1pt; PDHD 027409 evt 0 group02 59→59 (−2 points),
  group13 unchanged.  No whole-volume cluster loss anywhere.
- **Timing** (MABC perf, busiest event 39324 evt 4, group4567):
  separate 3.9 s, connect1 1.25 s, deghost 1.98 s.
