# QLMatching on the bright PDHD outlier (run 29107, evt 1015) — where the 12–15 GB / ~21 min goes, and how to cut it

**Assumption stated first** (per CLAUDE.md §1): the user's "run 27409, evt 1015" is read as **run 29107, DAQ event 1015** — the single bright, very-high-light PDHD outlier documented in
[`pdhd/docs/pdhd-pipeline-resource-profile.md`](../../pdhd/docs/pdhd-pipeline-resource-profile.md) §5 and
[`pdhd/docs/pdhd-pd-activity-per-event.md`](../../pdhd/docs/pdhd-pd-activity-per-event.md) §6 (454 flashes, 1.33 M PE). No run 27409 exists in this study. If a different event is meant, the per-phase method below still applies; the numbers are 29107/1015.

**Scope.** This is an *investigation + ideas* doc for review — **no code is changed**. It supersedes, for this event, the optimization lead in `pdhd-pipeline-resource-profile.md` §6, which guessed the cost was the `build_bundles` "flashes × clusters × points" loop. The logs say otherwise.

---

## 1. The measurement (from the production run's own `QLtiming` debug log)

The per-phase `QLtiming` markers already in `QLMatching.cxx` were emitted in the all-30-event reprocess. Event 1015 is matched as two per-volume runs (`group13` = +x APAs 1/3, `group02` = −x APAs 0/2), each a full single-APA pipeline. Verbatim from `pdhd/work/.batch_clus_029107_4.log`:

```
group13  build_bundles : nflash 454  ngroups 55  nbundles 24970  vis_loop  95758.6 ms over 28 146 638 candidate-points
group13  fit_round1    : nbundle 9017   nflash 450  nopdet 79  matrix_build 259681.8  lasso_solve 172904.4 ms  (X 9467 x 9467)
group13  fit_round2    : nbundle 154    nflash 94   nopdet 79  matrix_build     17.4  lasso_solve      0.4 ms  (X 154 x 154)

group02  build_bundles : nflash 454  ngroups 54  nbundles 24516  vis_loop  73596.8 ms over 52 523 260 candidate-points
group02  fit_round1    : nbundle 12534  nflash 428  nopdet 34  matrix_build 196841.2  lasso_solve 449630.8 ms  (X 12962 x 12962)
group02  fit_round2    : nbundle 213    nflash 112  nopdet 34  matrix_build     16.8  lasso_solve      1.3 ms  (X 213 x 213)
```

### Time

| phase | group13 (s) | group02 (s) | both (s) | share |
|---|--:|--:|--:|--:|
| `build_bundles` vis_loop | 95.8 | 73.6 | **169** | 14 % |
| `fit_round1` matrix_build | 259.7 | 196.8 | **457** | 37 % |
| `fit_round1` lasso_solve | 172.9 | 449.6 | **623** | 50 % |
| `fit_round2` (both phases) | <0.1 | <0.1 | ~0 | — |
| **total** | **528** | **720** | **~1248 (≈21 min)** | |

**The cost is `fit_round1`, not `build_bundles`.** matrix_build + lasso_solve = **87 %** of the wall; the SemiAnalyticalModel vis_loop the previous doc fingered is only **14 %**. `fit_round2` is free — round 1's strength prune collapses ~9–12 k bundles to ~150–210, so the second LASSO is trivial. **All the cost is the *first* LASSO on the full candidate universe.**

### Memory — accounting for the 12.4 GB peak

`fit_round1` builds dense Eigen matrices (`Ress::matrix_t = Eigen::MatrixXd`) and then the LASSO solver (`util/src/LassoModel.cxx`) makes more. Sizes from the logged dims (`nbeta = nbundle + nflash`):

| allocation | where | group13 (nbeta 9467) | group02 (nbeta 12962) |
|---|---|--:|--:|
| `P` dense `(nopdet·nflash) × nbeta` | `fit_round1` | 2.69 GB | 1.51 GB |
| `PT = P.transpose()` | `fit_round1` | 2.69 GB | 1.51 GB |
| `X = PᵀP + PFᵀPF` dense `nbeta²` | `fit_round1` | 0.72 GB | 1.34 GB |
| `X` copied **by value** into `Ress::solve(matrix, …)` | `Ress.cxx:7` | 0.72 GB | 1.34 GB |
| `Eigen::MatrixXd X = GetX()` (another copy) | `LassoModel::Fit` | 0.72 GB | 1.34 GB |
| `tripletList.reserve(nbeta²)` worst-case | `LassoModel::Fit` | 1.43 GB | 2.69 GB |
| **live peak (these overlap in `fit_round1`)** | | **≈9.0 GB** | **≈9.7 GB** |

Plus the ~0.8 GB clustering base and the event's charge/flash point clouds → the **~12.4 GB** the operator saw. The two groups run sequentially, so the event peak ≈ the larger group, not the sum.

**Two structural facts drive everything below:**
1. `P` (and `PF`) are **block-sparse**: a bundle is one (flash, cluster-group) pair, and its `pred_flash` column has nonzeros only in *its own* flash's `nopdet`-row block — 1 of `nflash` blocks. Density ≈ `1/nflash` ≈ **0.2 %**. They are stored **dense**.
2. The solver re-derives the Gram of the matrix it is handed and reserves a **fully-dense** triplet list for it, regardless of actual sparsity.

---

## 2. Why the obvious idea (split into independent sub-problems) does **not** help this event

The tempting structural fix — the LASSO separates across connected components of the bundle-coupling graph, so solve each component alone and never build the full `nbeta²` — **gives no win on 1015**, and the logged numbers prove it.

The coupling graph is bipartite: flash-nodes ↔ cluster-group-nodes, one edge per candidate bundle. From the log (group02 round 1): **428 flashes, 54 cluster-groups, 12 534 edges** out of 428×54 = 23 112 possible → **54 % density**. Average cluster degree = 12534/54 ≈ **232 flashes per cluster** (54 % of all flashes).

Pigeonhole: two clusters each adjacent to >50 % of the flashes **must** share a flash (232 + 232 − 428 = 36 > 0) → every pair of the 54 clusters is connected → **one component spanning ~all 428 flashes**. group13 (cluster degree ~36 %, ~164/cluster) gives expected ~60 shared flashes per cluster-pair → still **effectively one component**. Decomposition would shrink `nbeta` from 12962 to ~12900 — nothing.

**This is the key distinction to keep straight:** "the graph splits into disjoint pieces" (**false** here) is *not* the same as "the matrix is sparse *within* the one piece" (**true** here). Inside the single component, each bundle couples only to its ~29 same-flash bundles + ~232 same-cluster bundles ≈ **261 of 12962 → ~2 % nonzero**. The solver already *stores* its inner Gram (`XdX`) sparse — it just *builds* it densely and reserves it densely. That 2 % is the real lever, not graph partitioning.

> Keep component-decomposition as a *general-case* tool (events that genuinely fragment will benefit), but it is **demoted** for this event: a single dense component is exactly the regime it cannot help.

---

## 3. Ideas, ranked (evidence-driven)

Project rule ([[feedback_toggleable_behavior_changes]]): result-changing variants must be jsonnet-togglable, default OFF, existing configs bit-identical. The first two ideas are **bit-identical** (pure dataflow/storage); the rest change results and need a toggle + MC/data validation.

### A — Sparse linear algebra in `fit_round1`/`fit_round2` (bit-identical) — the headline

`P` and `PF` are ~0.2 % dense but stored as `Eigen::MatrixXd`. Build them as `Eigen::SparseMatrix` and form `X` / `y` with sparse products. Two payoffs, both on the dominant phase:

- **matrix_build (457 s, 37 %) collapses.** Today `X = PᵀP` multiplies two dense `(nopdet·nflash) × nbeta` matrices: ≈ `nbeta² · (nopdet·nflash)` ≈ 2.4×10¹² flops (group02), almost all `0·0`. Sparse `PᵀP` only touches same-flash column pairs: ≈ `Σ_flash nopdet · (bundles_per_flash)²` ≈ 428 · 34 · 29² ≈ **1.2×10⁷** flops — ~5 orders of magnitude fewer. matrix_build should drop to well under a second.
- **The `P`/`PT` memory spike (2.69 + 2.69 GB on group13) disappears** — sparse `P` at 0.2 % is a few MB. And `X = PᵀP + PFᵀPF` comes out **sparse (~2 %)**, so it can stay sparse downstream.

The one subtlety: `Ress::solve` and `LassoModel` take a **dense** `matrix_t`. Realising the matrix_build and `P/PT` wins needs no solver change (form `X` densely from sparse `P` at the boundary). Realising the **solve-time** win needs the solver to accept a sparse design matrix (idea B2). Recommend doing the sparse assembly first (self-contained, bit-identical, ~74 % of the *time* is matrix_build+solve and matrix_build alone is 37 %).

### B — Remove the redundant dense copies + fix the worst-case reserve (bit-identical) — the memory headline

In `util/src/LassoModel.cxx::Fit` and `util/src/Ress.cxx`, three allocations are pure waste at `nbeta ≈ 13 000`:

1. **`tripletList.reserve(size_t(nbeta) * size_t(nbeta))`** (`LassoModel.cxx`) = **2.69 GB** reserved (group02) for a matrix that ends up ~2 % full. The code comment shows this was deliberately raised from `nbeta²/2` to `nbeta²` to avoid one growth realloc "at hd-max sizes" — i.e. it trades **2.7 GB** to save one `std::vector` doubling. Replace with a **counted** reservation: the structural nonzero count of `XdX` is cheap to pre-tally (it is the number of (bundle,bundle) pairs at graph-distance ≤ 2), or reserve per-column and let Eigen pack. Caveat below.
2. **`Eigen::MatrixXd X = GetX();`** (`LassoModel.cxx::Fit`) — a full dense copy of the `nbeta²` matrix (1.34 GB). Use a `const&`.
3. **`Ress::solve(Ress::matrix_t matrix, …)`** takes the matrix **by value** (`Ress.cxx:7`) — another 1.34 GB copy at the call. Take it by `const&`.

Together ≈ **−5.4 GB** on group02 with zero change to the result. **Blast-radius flag:** `LassoModel.cxx`/`Ress.cxx` are shared with imaging (`img/src/CSGraph.cxx:190`, `img/src/BlobSolving.cxx:121`). Any "bit-identical" edit here must be re-validated against an imaging run, not just QLMatching. The reserve heuristic in particular has a known time/space tradeoff at "hd-max sizes" — prefer a *counted* reserve (exact, no realloc, no over-allocation) over simply lowering the constant.

### B2 — Let the solver consume a sparse design matrix (bit-identical, larger lift) — the solve-time lever

`lasso_solve` is 623 s (50 %). Inside `Fit`, `XdX(i,j) = X.col(i).dot(X.col(j))` over all `i ≤ j` is an **O(nbeta³)** dense Gram reconstruction (~10¹² ops at nbeta 12962 → the 450 s). If `LassoModel` accepted `X` as `Eigen::SparseMatrix`, those column dot-products would skip the ~98 % zero structure and `XdX` would assemble directly from the sparse pattern. This is the same arithmetic on the same nonzeros → bit-identical, but it edits the **shared** solver interface (overload, don't replace) and must pass the imaging regression. Bigger than A/B; do it only if A+B+C leave the solve as the bottleneck.

> A deeper observation, **explicitly out of scope** (do not implement without a separate study): QLMatching pre-forms the *normal equations* `X = PᵀP`, `y = PᵀM` and hands `X` to the LASSO, which then re-Grams it (`XdX = XᵀX`). Every other `Ress::solve` caller (`CSGraph`, `BlobSolving`) passes the **design** matrix directly. Passing `P`/`M` straight through would skip both Grams — but it changes the conditioning the L1 penalty sees and is a ported, validated MicroBooNE-ToyMatching convention. Flagging only; the matching tune depends on it.

### C — Shrink `nbeta` at the source with a pre-fit cull (toggle, changes physics) — the real handle on solve time

`fit_round1` enters with **12 534** candidate bundles, ~233× the **154–213** that survive round 1. Because cost scales like `nbeta²`–`nbeta³`, halving the candidate set is a 4–8× win — far more leverage than any constant-factor dataflow fix. Most candidates are a noise flash paired with a geometrically impossible cluster. Two existing knobs already gate exactly this and are simply **OFF for PDHD**:

- `require_containment` — drop a bundle whose cluster leaves the TPC box once the flash-T0 x-offset is applied.
- `reject_overpred` (`overpred_total_ratio`/`overpred_maxch_ratio`) — drop a bundle predicting far more light than the flash shows.

A third, cheap, PDHD-appropriate gate: a **flash-time / drift-x window** — a flash at time *t* implies a drift-x; clusters whose x is outside `[anode, cathode]` for that *t* cannot match, and the candidate need never be built (cut it in `build_bundles` before the vis_loop, which also trims the 14 % there). All three are post-`examine_bundle` physics changes → default OFF, validate on MC + the data hand-scans (the `ql_scan` tool already exists). This is the lever that attacks the **O(nbeta³)** solve at its root rather than its constant.

### D — Coarsen the vis_loop for very large blobs (toggle, ~14 %) — secondary

`build_bundles` evaluates the SemiAnalyticalModel per charge point: 28–52 M point-opdet evaluations. For large blobs the per-point light is smooth, so sampling 1-in-k points and scaling would cut this ~k×. Approximation → toggle, validate. Only 14 % of the wall and the SBND doc's "fix A" already covers the safe bit-identical micro-opts, so this is a distant secondary.

### E — Component decomposition (general case only) — demoted

Exact and free *when the graph fragments*, but §2 shows evt 1015 is a single dense component, so it is **inert here**. Worth a cheap guard (if components > 1, solve per component) for events that do fragment, but it does not touch this outlier.

---

## 4. Recommendation

1. **A + B together** (bit-identical, no toggle, no physics risk): sparse `P`/`PF` assembly kills matrix_build (37 %) and the `P/PT` spike; removing the solver copies + the `nbeta²` reserve kills ~5 GB. Expected: **~12 GB → ~3–4 GB, ~21 min → ~12 min**, output byte-for-byte unchanged. Validate against **both** QLMatching and an imaging run (shared util).
2. **C** (toggle, default OFF) for the solve-time tail: a flash-time/drift-x or containment pre-fit cull, validated on MC + `ql_scan` data. This is what bends the `nbeta³` curve and helps *every* busy event, not just 1015.
3. **B2 / D** only if 1–2 leave a bottleneck. **E** as a cheap general-case guard. The gram-of-gram note stays an observation, not a task.

The headline correction for the next reader: on this event the LASSO `fit_round1` — its dense Gram build and the solver's dense copies/reserve — is the whole story; `build_bundles` is a sideshow, and splitting the problem into independent pieces does not apply because there is only one.

---

## 5. Implemented — Option A: sparse `P`/`PF` assembly (`match/src/QLMatching.cxx`)

**Re-baseline.** Measured against the updated flash reconstruction (HEAD: ADC-saturation
veto + per-PD `min_fired_pe`, `min_fired_pds:5`/`min_total_pe:20`). Those cuts barely move evt
1015 — still **443 flashes** — so it remains the stress case. The dense baseline (this build):
**14.32 GB / 1148 s**, `fit_round1` matrix_build **420.6 s** + lasso_solve **522.0 s**.

**Change.** `fit_round1`/`fit_round2` now fill `P`/`PF` as `Eigen::Triplet` lists and realise
them either **dense** (default — `Ress::matrix_t(P_sp.toDense())`, then the *unchanged*
`X = PᵀP + PFᵀPF` GEMM) or **sparse** (`m_sparse_lasso` — sparse `PᵀP`/`PᵀM`, then `X.toDense()`
into the unchanged `Ress::solve`). New `sparse_lasso` config flag, **C++ default `false`**
(SBND / every existing config byte-identical); PDHD `cfg/.../pdhd/qlmatching.jsonnet` sets it
`true`. Sparse and dense matrix products accumulate FP sums in different orders, so this is a
toggle per [[feedback_toggleable_behavior_changes]], not an unconditional refactor.

**Validation (all 30 events of run 29107, vs the dense baseline reference):**

| path | mabc production output | calib diagnostic (`strength`) |
|---|---|---|
| dense default (`sparse_lasso:false`) | **29/29 byte-identical** | **29/29 byte-identical** |
| sparse (`sparse_lasso:true`, PDHD) | **30/30 byte-identical** | 25/30 byte-identical; **5/30 strength-only** |

The dense path is byte-identical down to the diagnostic `strength` float → the triplet refactor
changed nothing. The sparse path leaves every **production** output (cluster t0, matched-flash
index, `op` predictions) byte-identical on all 30 events; the only movement is ULP-level drift in
the diagnostic LASSO `strength` of 5 busy events (assignments, flags, KS/chi2 unchanged) — exactly
the documented reason it is gated.

**Resource win, evt 1015 (sparse vs dense, same build):**

| metric | dense baseline | + Option A (sparse) | Δ |
|---|--:|--:|--:|
| `fit_round1` matrix_build (both groups) | 420.6 s | **0.5 s** | **−99.9 %** |
| `fit_round1` lasso_solve (both groups) | 522.0 s | 524.6 s | ~0 (Option B's target) |
| peak RSS | 14.32 GB | **11.03 GB** | **−3.3 GB** |
| total QLMatching wall | 1148 s | **733 s** | **−36 %** |

Option A does exactly what §1/§3-A predicted: the dense `P`/`PT` spike (~3 GB here) and the
O(nbeta²·nopdet·nflash) Gram build vanish. What remains — the **~11 GB** and the **~520 s
lasso_solve** — is the dense `X`, its four copies, and the `nbeta²` triplet reserve inside
`LassoModel`, i.e. **Option B**.

---

## 6. Implemented — Option B: drop the redundant dense copies + right-size the reserve (`util/`)

**Change (all bit-identical — pure dataflow/storage, no arithmetic touched):**
- `Ress::solve(matrix_t matrix, …)` → `const matrix_t&` (`util/src/Ress.cxx`,
  `util/inc/WireCellUtil/Ress.h`): the GB-scale Gram is no longer copied at the call.
- `LinearModel::SetData`/`SetX` and the `ElasticNetModel::SetX` override → `const Eigen::MatrixXd&`
  (`util/inc/WireCellUtil/{LinearModel,ElasticNetModel}.h`): the one surviving copy is the `_X`
  member. `LassoModel::Fit` / `ElasticNetModel::Fit` `Eigen::MatrixXd X = GetX();` →
  `const Eigen::MatrixXd& X` (`util/src/{LassoModel,ElasticNetModel}.cxx`). Net: the `nbeta²`
  matrix is copied **once** instead of four times.
- `LassoModel::Fit` triplet/`XdX` reservations: the `tripletList.reserve(nbeta²)` (~2.7 GB) and
  `XdX.reserve(nbeta/2 per col)` (~1 GB) are right-sized by a one-pass `nnz(X)` check — a **dense**
  X (imaging) keeps the historical nbeta-scale reserve unchanged; a **sparse** X (the QL Gram)
  reserves a small seed and grows. Capacity-only → the assembled `XdX` and the fit are unchanged.

These are **shared util** (used by imaging `CSGraph`/`BlobSolving`), so default-OFF toggling does
not apply — they must be unconditionally bit-identical, and are.

**Validation:**
- **Imaging regression** (run 29107 evt idx 9, `wct-img-all`/`run_img_evt.sh -d on`): the four
  `clusters-apa-apa{0..3}-ms-active` archives are **byte-identical** pre- vs post-B (every inner
  numpy member) — the shared-util edits do not perturb imaging.
- **QL** (29 events of run 29107 vs the post-A reference): **29/29 mabc AND calib byte-identical**,
  including the raw `strength` float — B is numerically a no-op, as intended.

**Resource win, evt 1015:**

| metric | dense baseline | + A (sparse) | + A + B | A+B Δ vs baseline |
|---|--:|--:|--:|--:|
| `fit_round1` matrix_build | 420.6 s | 0.5 s | 0.4 s | **−99.9 %** |
| `fit_round1` lasso_solve | 522.0 s | 524.6 s | 541.5 s | ~0 (compute unchanged) |
| peak RSS | 14.32 GB | 11.03 GB | **8.41 GB** | **−5.9 GB (−41 %)** |

B removes ~2.6 GB on top of A (three redundant `X` copies + the multi-GB reserves). The **remaining
~8.4 GB / ~540 s** is structural and outside A+B's scope: the dense `X` still lives twice (the
`fit_round1` local + the `_X` member) feeding the **unchanged** `Ress::solve`, and the O(nbeta³)
coordinate-descent Gram-of-Gram is the solve time. Cutting those is **B2** (teach `LassoModel` to
consume a sparse `X` — removes both the second dense copy and the dense Gram reconstruction) and
**C** (a pre-fit cull to shrink `nbeta` at the source); both are follow-ups, not part of A+B.

---

## 7. Remaining headroom after A+B (what's left, and is it result-preserving?)

After A+B evt 1015 is **8.41 GB / ~733 s**, of which **~540 s is `lasso_solve`** and the memory
floor is the dense `X` (still materialised for the solver) plus the charge point clouds. Measuring
`lasso_solve` vs `nbeta` across run-29107 events (from the logs, no new runs):

| nbeta | 1429 | 2225 | 3082 | 4119 | 9289 | 12688 |
|---|--:|--:|--:|--:|--:|--:|
| solve (s) | 0.4 | 2.1 | 6.1 | 16.9 | 169.8 | 371.7 |

This is **≈ O(nbeta³)** (fit exponent ~2.8–3.1), which fingers the **dense Gram-of-Gram build**
in `LassoModel::Fit`: `XdX(i,j) = X.col(i).dot(X.col(j))` over all `i ≤ j` runs `nbeta²/2` dot
products, each over **dense** columns of a matrix `X` that is only ~2 % non-zero. The coordinate-
descent iterations (over the already-sparse `XdX`, via `InnerIterator`) scale only ~nbeta² and are
the small term. So essentially all of the remaining solve time is dense arithmetic on ~98 % zeros.

**B2 — feed the solver a sparse `X` (result-preserving, under the existing `sparse_lasso` toggle).**
The lever: `fit_round1` already has the sparse Gram `Xs` (it currently calls `Xs.toDense()` only to
satisfy `Ress::solve`'s dense `matrix_t`). If `LassoModel` accepted a sparse `X`, it would (i) never
materialise the dense `X` (removes the last ~2.7 GB dense matrix and its `_X` copy), and (ii) build
`XdX` from the sparse pattern — skipping the ~98 % zero dot-products and only visiting column pairs
that actually overlap. Same non-zero arithmetic → byte-identical modulo the same ULP drift A already
carries, so it lives under `sparse_lasso` (default OFF, PDHD on) and is validated the same way
(all-30 mabc byte-identical, calib strength-only). Expected: the ~540 s solve → tens of seconds and
peak RSS → ~4–5 GB. **This is the single largest remaining result-preserving win.** Cost: it edits
the shared `LassoModel`/`Ress` interface (add a sparse overload, keep the dense path for imaging),
so it needs the imaging regression again.

**C — pre-fit cull to shrink `nbeta` (NOT result-preserving → toggle + physics validation).** Orthogonal
and compounding: `fit_round1` enters with ~12 k candidate bundles vs ~200 survivors. A flash-time /
drift-x or containment cull (the existing `require_containment` / `reject_overpred` knobs, off for
PDHD) cuts `nbeta` at the root, and because the solve is O(nbeta³) even a 2× cut is ~8×. But it drops
candidate matches → it changes physics, so it is a default-OFF toggle validated on MC + `ql_scan`
hand-scans, *not* a byte-identical change.

**Smaller / not worth it.** Within the current math everything else is minor: `ydX`/`norm` are
O(nbeta²); the XdX symmetry is already exploited (compute `i≤j`, mirror). The deeper structural note
(QLMatching hands the solver the pre-formed normal-equations `X = PᵀP` and the solver re-Grams it to
`XᵀX`) is a genuine redundancy, but undoing it (pass `P`/`M` directly) changes the conditioning the L1
penalty sees — it is the ported, tuned MicroBooNE convention, so it is **not** result-preserving and
stays an observation.

**Bottom line:** the one big *result-preserving* lever left is **B2** (sparse `X` through the solver),
gated by the same `sparse_lasso` flag; **C** is a larger lever still but only as a validated physics
toggle. A+B already removed the matrix-build cost and ~6 GB; B2 would take the solve from minutes to
seconds and the memory to a few GB. **B2 is now implemented — see §8.**

---

## 8. Implemented — Option B2: feed the solver a sparse `X` (`util/` + `match/`)

Realises §7. Under `sparse_lasso`, `fit_round1`/`fit_round2` now pass the sparse Gram `Xs` **directly**
to a new `Ress::solve(const Eigen::SparseMatrix<double>&, …)` overload (no `toDense()`), and
`LassoModel::Fit` gained a sparse path that builds `norm` / `ydX = Xᵀy` / `XdX = XᵀX` by **sparse
products** instead of the dense `X.col(i).dot(X.col(j))` loop. The dense path (imaging, SBND, any
caller of the dense `solve`) is wrapped verbatim in the `else` branch — untouched.

- `util/inc/WireCellUtil/{Ress.h,LassoModel.h}`, `util/src/{Ress.cxx,LassoModel.cxx}`: sparse `solve`
  overload → `LassoModel::SetXsparse` (stores `_Xsp`, sets `_use_sparse`); `Fit` branches on it.
- `match/src/QLMatching.cxx`: the `sparse_lasso` branch builds `Xs` and calls the sparse `solve`; the
  dense branch is byte-for-byte the old code.

**Validation:**
- **Imaging** (idx 9, 4 anodes) **byte-identical** pre/post — the dense `Fit` path survived the refactor.
- **QL** (29 events vs the dense baseline): **29/29 mabc production byte-identical**; the diagnostic
  `strength` float is now strength-only different on all 29 (the sparse `XᵀX` accumulates in a
  different FP order than the dense dot loop — more events move than under A alone, but still *only*
  `strength`; every assignment, flag, KS and chi2 is identical). Result-preserving, as intended.

**Resource win, evt 1015 — the full A → B → B2 progression:**

| metric | dense baseline | + A | + A+B | **+ A+B+B2** |
|---|--:|--:|--:|--:|
| `fit_round1` matrix_build | 420.6 s | 0.5 s | 0.4 s | **0.14 s** |
| `fit_round1` lasso_solve | 522 s | 525 s | 541 s | **9.2 s** |
| peak RSS | 14.32 GB | 11.03 GB | 8.41 GB | **3.94 GB** |
| total QLMatching wall | 1148 s | 733 s | 748 s | **220 s** |

B2 takes the solve from **~540 s to ~9 s (≈59×)** — group13 170 s→2.7 s, group02 372 s→6.5 s — and
drops another **4.5 GB** (the dense `X` is gone; `XdX` is built from the sparse pattern). **Combined
A+B+B2: 1148 s → 220 s (5.2×) and 14.32 GB → 3.94 GB (3.6×) on evt 1015, production byte-identical.**

What's left is now genuinely structural for the *solver*: the ~9 s solve is the coordinate-descent
iterations over the sparse `XdX` (~O(nbeta²·sweeps)). **But the solver is no longer the bottleneck —
see §9: a fresh re-profile shows `build_bundles` (the per-point visibility loop) is now ~90 % of
QLMatching.** The next result-preserving lever is there, not in the LASSO.

---

## 9. Post-B2 re-profile (run 29107, evt 1015) — the bottleneck is now `build_bundles`, not the solver

Re-measured on the current production build (A+B+B2 all on), evt 1015 = charge idx 4, joint matcher
(two per-drift-side runs). Numbers straight from the run's own `QLtiming` debug log
(`pdhd/work/029107_4/wct_clus_029107_4.log`):

| phase | run A | run B | total | share of QL |
|---|--:|--:|--:|--:|
| **`build_bundles` vis_loop** | 106.2 s | 66.3 s | **172.5 s** | **~90 %** |
| `cull_cross_tpc` (pairing loop) | — | — | 19.4 s | ~10 % |
| `fit_round1`+`fit_round2` (LASSO) | 0.52 s | 0.03 s | 0.56 s | ~0.3 % |
| output (tensor build) | — | — | 0.38 s | ~0.2 % |
| **operator total** (prefit+cull+fit+out) | | | **~197.5 s** | |

(process wall 236 s / peak RSS 1.38 GB; the gap to 197 s is frame I/O + clustering + finalize.) **The
A+B+B2 work made the LASSO free, which is exactly why the picture flipped:** the `build_bundles`
vis_loop — the ported `SemiAnalyticalModel` per-point/per-opdet kernel — was 169 s ≈ 14 % of the *old*
1248 s; it is unchanged at ~172 s but is now ~90 % of the *new* ~197 s. The §3 "Option D (~14 %)"
ranking was the pre-optimization fraction; it is now the headline. `cull_cross_tpc` is the only other
nontrivial cost (19.4 s over **1.86×10¹⁰** point-pairs ≈ 1.1 ns/pair — already at the scope-hoisted
floor, §`xtpc_pair_consistent`; the math is irreducible without an AABB/k-d pre-filter, which is a
separate physics-touching change).

**Where the vis_loop time is wasted (the result-preserving lever).** `build_bundles` builds a bundle
for **every** (flash × cluster-group) pair — run A is 451 flashes × 52 groups = 23 452 bundles, run B
451 × 54 = 24 354 — and computes the full per-point visibility **before** the candidate filters run.
`require_containment` is already hoisted *above* the vis_loop (it `continue`s before the loop), but the
opaque-cathode **cross-side mismatch drop** (`cross_side_mismatch_drop`) runs *after* the loop, even
though its decision (flash lit-side vs the run's fixed cluster side, plus the pre-loop `at_x_boundary`
flag) **never touches the predicted light**. PDHD feeds the full global flash list to *each* per-side
run, so a large share of each run's bundles are cross-side and are dropped only after paying full
visibility cost (survivors `pre_bundles`: run A 4831 of 23 452, run B 1433 of 24 354). Hoisting that
check above the vis_loop — the same pattern `require_containment` already uses — skips the visibility
of bundles that are discarded anyway → **byte-identical matching, lower wall time**. The two in-loop
side effects are safe: `total_charge_blob`/`total_charge_point` are debug-log only, and the calib-dump
bundle table already applies the same `cross_side_mismatch_drop` before reading any pred. The realized
saving (= the vis_loop time spent on cross-side-but-contained bundles) is measured directly by the
`vis_loop` A/B in §10.

## Reproduce

```
# the per-phase numbers above come straight from the production debug log:
grep -E 'QLtiming (fit_round|build_bundles)|preselected' \
    pdhd/work/.batch_clus_029107_4.log

# matrix/solver sources cited:
#   match/src/QLMatching.cxx   fit_round1 (~1286), fit_round2 (~1388)
#   util/src/Ress.cxx:7        solve(matrix BY VALUE, …)
#   util/src/LassoModel.cxx    Fit(): GetX() copy + tripletList.reserve(nbeta^2)
# graph-density / single-component arithmetic: §2, from the logged nflash/ngroups/nbundle.
```
