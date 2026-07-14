# MC evt12 — a real cathode crosser whose TPC0 half matches the WRONG flash (Q/L)

**Status: FIXED** (see "Fix (shipped)" below). Diagnosis of the SBND MC event-12
cross-TPC case the hand-scan flagged. The two
halves are **TPC0 cluster ident 2** (user pos (-219.7, 23.9, 436.8), raw) and
**TPC1 cluster ident 3** (user pos (234.1, 80.2, 416.7), raw) — Bee/display ids
"25" and "2". (NB cluster idents are per-APA; an earlier draft of this doc
analysed the wrong pair — group-29 cluster3↔cluster5 — see the appendix.)

Question: *flagging failure, or something else?* Answer: **a cull/pairing
ORDERING problem** — `cull_inconsistent` collapses TPC0 cluster 2 onto a spurious
flash *before* the cross-TPC pairing runs, so the real crosser bundle (which would
pass xtpc scenario-1 at d=1.31 cm) is gone before it can be flagged, and TPC0
matches the wrong flash.

All numbers from a `-calib -cathode-diag` rerun with added debug logs; see
`sbnd_xin/work/ql_evt12/wct_ql_evt12.log`.

## The crosser

One physical track crosses the cathode at t0 ≈ **-1350.66 µs** (flash **group 1**):
- **TPC1 cluster 3** → flash gid 1000005 (t=-1350.66): ks=0.039, chi2=2097,
  **strength 0.996, auto-matched, high-consistent**. Correct. ✓
- **TPC0 cluster 2** → should match flash gid 5 (t=-1350.66) but instead matched
  flash gid 14 (t=**-1190**, group 4). ✗

Both raw clusters sit *outside* their TPC box (TPC0 raw x∈[-234,-212] beyond the
−200 anode; TPC1 raw x∈[212,234] beyond the +200 anode). After the **correct** T0
correction (off=±211 cm at t=-1350.66) they move to the cathode and meet:

```
TPC0 clus2 @ -1350.66:  corrected x∈[-23.0, -0.5]  y∈[-4.4,36.8]  z∈[429.5,444.9]
TPC1 clus3 @ -1350.66:  corrected x∈[  0.5, 23.0]  y∈[34.9,80.2]  z∈[416.4,431.7]
closest approach d = 1.31 cm   →  xtpc scenario-1 (dmax=5cm) PASS
```

So geometrically this is an easy, unambiguous cross-TPC match. The user's mental
model — "outside the box before T0, a clean crosser after T0" — is exactly right.

## Why TPC0 cluster 2 matches the wrong flash

TPC0 cluster 2 builds a bundle on every candidate flash. Two matter:

| flash | t0 (µs) | grp | corrected x | ks | chi2/ndf | high-consist | outcome |
|-------|---------|-----|-------------|-----|---------|--------------|---------|
| gid 5  | -1350.66 | 1 | [-23.0,-0.5] (cathode) | 0.113 | 3027/49 | **no** | str 0, not selected |
| gid 14 | -1190    | 4 | [-48,-26]              | 0.064 |  532/49 | **yes** | str 0.47, auto-matched |

The **true crosser bundle (gid 5) is window-truncated and under-predicts → not
high-consistent**; a coincidental bundle on flash 14 *is* high-consistent. Then:

```
QLCULLINC apa0 cluster 2 drop bundle flash 5 (cluster kept high-consistent on flash 14)
```

`cull_inconsistent` (run in prefit) keeps cluster 2's single high-consistent
bundle (flash 14) and **drops every other bundle of that cluster — including the
gid-5 crosser bundle.** `cull_cross_tpc` runs *after* that, so when it gathers
candidates, TPC0 cluster 2 only exists at t=-1190; its true partner TPC1 cluster 3
is at t=-1350.66; the two are **not flash-coincident → never paired → never
flagged xtpc** (confirmed: no `QLXTPC coincident 0/2 1/3` line; the only cluster-2
coincidence is the wrong `0/2 1/1` at t=-1190). The crosser is broken and cluster 2
keeps the spurious -1190 match.

## Why the xtpc flag could not save it even if pairing ran (the user's
"shouldn't flag_xtpc make it high-consistent?")

Two reasons, both real:
1. **Ordering:** the gid-5 bundle is culled by `cull_inconsistent` *before*
   `cull_cross_tpc`, so there is nothing left to flag.
2. **The flag is too weak even if set:** `flag_xtpc_consistent` only culls a
   cluster's rivals that are *neither* high- nor xtpc-consistent. Cluster 2's
   rival on flash 14 *is* high-consistent, so it would survive, and the flag does
   NOT promote the gid-5 bundle in the ladder/LASSO/strength. So flagging alone
   would not flip cluster 2 from flash 14 to flash 5. The user's expectation —
   that an xtpc-consistent pair should be *promoted* like a high-consistent match
   — is **not** how the code currently works (`flag_xtpc_consistent` is read only
   inside `cull_cross_tpc` + the dump; grep-confirmed).

## Fix (shipped) — flag before culling + xtpc-priority steering

Both root causes are addressed. The fix is gated by the existing `m_xtpc_flag`
(SBND-on, default OFF), so pdhd/pdvd/uboone and any `xtpc_flag:false` config stay
bit-identical — no new config toggle.

1. **Ordering — flag xtpc BEFORE culling.** `cull_inconsistent` was moved out of
   `run_one_apa_prefit`; `operator()` now runs *build all APAs → `cull_cross_tpc`
   (flag only) → `cull_inconsistent` all APAs → fit all APAs*. The cross-TPC pairing
   therefore sees the full pre-cull bundle set, so cluster 2's gid-5 crosser bundle is
   present and pairs with cluster 3 at **d = 1.31 cm → scenario 1**. To bound the
   pairing cost on the larger set, candidates are restricted to bundles with
   `flag_at_x_boundary || flag_window_truncated` (both crosser halves carry both).
2. **Steering — xtpc-priority cull, scenario-1-scoped.** `cull_cross_tpc` now only
   *flags* (it sets `flag_xtpc_consistent`, and `flag_xtpc_scenario1` for the tight
   case); its old rival-cull block was folded into a unified `cull_inconsistent`:
   *a cluster with a scenario-1 xtpc bundle keeps that bundle and drops ALL others
   (incl. its high-consistent bundle on the wrong flash)*; otherwise the historical
   keep-(high|xtpc) rule applies (so scenario-2 keeps today's weaker behaviour and the
   non-xtpc path is bit-identical). Cluster 2 then keeps only gid 5; the LASSO
   charge-conservation DOF gives it strength 0.989 and it matches the crosser flash.

**Validation (10 data + 10 mc hand-scan, auto_selected vs `.scan_state` selected,
keyed (flash_gid, main_cluster); A/B vs the pre-fix binary):** the fix is fully
surgical — across all 20 events the **only** change is mc evt12, where cluster 2
moves from flash 14 (wrong) to flash 5 (the crosser): `added (5,2) / removed (14,2)`.
Agreement **MC 100→101 (+1, missing 12→11), DATA 95/95 unchanged**; every other event
byte-identical in `auto_selected`. evt12 now flags both halves xtpc scenario-1
(cluster 2 strength 0.989, cluster 3 strength 0.996).

Edits: `match/src/QLMatching.cxx` (operator() reorder; `run_one_apa_prefit`/`run_one_apa`;
unified `cull_inconsistent`; `cull_cross_tpc` flag-only + candidate restrict;
`xtpc_pair_consistent` returns scenario code; calib dump adds `xtpc_scenario1`),
`match/inc/WireCellMatch/TimingTPCBundle.h` (`flag_xtpc_scenario1`),
`match/inc/WireCellMatch/QLMatching.h` (signature/comment; removed the now-unused
`consistent_bundles`). The cleaner long-term answer remains joint matching (combine
TPC0+TPC1 charge for the shared flash; `project_group_aware_ql_matching`).

## Instrumentation added (debug-level, output-neutral)

`match/src/QLMatching.cxx`: per-pair + coincident-pair + candidate logging in
`xtpc_pair_consistent`/`cull_cross_tpc`; removal logging in `cull_inconsistent`;
per-candidate/decision logging in `rescue_empty_flashes`. All `log->debug`, so
production output is unchanged.

## Appendix — the group-29 red herring (different track)

Flash group 29 (t≈1224 µs) contains a *separate*, genuinely doubly-truncated
crosser (TPC0 cluster 3 ↔ TPC1 cluster 5). Its T0-corrected closest approach is
406 cm (both near-cathode ends truncated away → halves hug opposite anodes),
collinear (3.6°) but past the scenario-2 300 cm purity cap. That is a real but
distinct efficiency limit of the xtpc distance cut; it is **not** the case the
hand-scan flagged here.
