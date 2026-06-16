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
