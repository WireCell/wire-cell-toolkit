# TrackFitting `spd_solver`: CG + IncompleteCholesky for the four normal-equation solves

Status: **knob off = byte-identical** (gates below); **knob on = NOT bit-identical** (every event's fit moves in its last bits;
the physics read on the pr_ml tuning-event fit lane is neutral). Default off everywhere; no detector config sets it.

## Repro

```bash
# the knob: the lane's TrackFitting parameter file with one more key (TaggerCheckNeutrino trackfitting_config_file)
#   { ..., "spd_solver": 1 }            # optional: "spd_solver_tol": 1e-10 (the default)
# the unit test (new symbol: link against this build's own libraries, see "Verification")
./build/clus/wcdoctest-clus -tc="*spd_solver*"
# the census of the solves (log only, pre-existing): WCT_DQDX_SOLVE_CENSUS=1 with -L debug prints one DQDX_SOLVE line per
# dQ/dx solve; with the knob on it ends in " solver=cgic" (or " solver=cgic_failed_bicgstab" when the legacy solve ran)
```

## Symptom

On the pr_ml doc 87b tuning-event fit lane (external graphs, `pr_fit_only = "external"`), one FD-HD event
(`nue_b11099_e2`) took 1,270-2,032 s on the 1e5 pair head's G-ki graphs, against 11 s for the median event. Every lane arm
waited on it. Its dQ/dx fit also left 7,617 fit points with no charge.

## Root cause

`dQ_dx_multi_fit` solves the normal equations A x = b, A = R^T M R + F^T F, which are symmetric positive (semi-)definite.
It used `Eigen::BiCGSTAB` with Eigen's defaults: tolerance at machine epsilon, cap 2n iterations, inherited from the
prototype (`PR3DCluster_multi_dQ_dx_fit.h:803`). The other three solves (`fit_point`, `trajectory_fit`, `dQ_dx_fit`) use the
same call.
- On the slow event the graph puts 1,704 object segments (trunks on shower and low-energy pieces) over one large shower. The system is n = 25,024 with 18.0 M
  nonzeros (median 787 couplings per point, against 25 on a typical graph) and condition number about 2.5e8.
- BiCGSTAB ran to its cap (50,048 iterations, 1,148 s) and stopped at a relative residual of 8.8e-4.
- The answer was off by up to 680 times the largest exact value: 7,587 negative points, where the exact solution
  (SimplicialLDLT) has 4,624.

## Why it hid

- On ordinary systems BiCGSTAB converges, the census reports nothing unusual, and the solve takes milliseconds.
- A capped solve still returns a finite vector, and `solver.error()` is not NaN, so the NaN fallback never fires.
- Only an event with a very dense, ill-conditioned dQ/dx system shows it, and then as time, not as an error.

## Fix

Two TrackFitting parameters, both read through `set_parameter` / `get_parameter`:
- `spd_solver` (default 0 = off);
- `spd_solver_tol` (default 1e-10, read only when on).

The public static `TrackFitting::spd_cg_ic_solve` runs `Eigen::ConjugateGradient<SparseMatrix, Lower|Upper,
IncompleteCholesky<double, Lower, AMDOrdering<int>>>` from the same guess, at that tolerance (Eigen's cap of 2n).
- Before it, every diagonal entry is stored explicitly, because IncompleteCholesky reads the first stored entry of each lower
  column as its diagonal.
- It returns false, leaving x at the guess, when the preconditioner fails or the answer or its error is not finite.

**With the knob on**, each of the four solves calls the helper first. If the helper returns false, the legacy BiCGSTAB statements
run unchanged.

**With the knob off**, the legacy statements run as before.
- The census and dump lines print the same counters in the same format.
- `trajectory_fit`'s leave-one-plane-out debug path reads `solver.error()` only when BiCGSTAB ran.

Eigen only; no new dependency.

## Verification

- **Unit test:** `clus/test/doctest_trackfitting_spd_solver.cxx`, 4 cases, 40 assertions.
  - The defaults and the round trip.
  - A 24 × 24-grid smear-plus-smoothness system (n = 576). BiCGSTAB at the legacy settings runs to its cap and misses the
    direct answer by more than 10 % of its largest value; the helper converges and matches it to 1e-4.
  - The guards: an empty column, a zero right-hand side, a NaN in A, and an empty or mis-sized system.
  - Singular 3 × 3 position systems (rank 2 and rank 1, with empty columns): the helper gives the legacy solve's own answer.
- **The clus suite:** 550 of 550 cases (1 skipped), 1,112,416 assertions.
  - The test calls a new symbol, so `wcdoctest-clus` must link against this build's libraries ahead of an older
    `local/lib`.
- **Knob absent, byte-identical** (`/ab-verify`). Side A was pin3 (toolkit `0ed0754d`) and side B this change, at `0ed0754d`.
  pr_ml doc 87b has the files.

| gate | labels | result |
|---|---|---|
| uBooNE qlport two-gate | `d87cgA_ub` / `d87cgB_ub` | zips 35 / 35 content-identical; tagger 34 / 35, 6805 bistable on both sides (3 repeats each, doc 79's two values on each side) |
| PD-HD / PD-VD `abtest` clus-only, events.txt | `d87cgA_pd` / `d87cgB_pd` | OVERALL: PASS |
| SBND PR chain, 60 events | `work-d87cg{A,B}-mcp60` | 120 / 120 archives byte-identical |
| pr_ml dev lane null, 8 events | `n2_pin3off` / `n2_pin4off` | 8 / 8 equal by value |

- **Knob on, smoke** (the slow event, census on):
  - the heavy solve takes 993 CG iterations, error 9.7e-11, 63.7 s (BiCGSTAB: 50,048 iterations, error 8.8e-4, 1,148 s);
  - the no-charge fit points fall from 7,617 to 4,660;
  - the job's CPU time falls from 1,794 s to 517 s.
- **Knob on, physics read** (pr_ml doc 87b read CG1: the 1e5 pair head's two seeds on 159 tuning events, the paired fit read):
  - no cell for or against on either seed; HOLDS, NEUTRAL;
  - the slow event takes 511 / 456 s instead of 1,753 / 1,270 s.
- **What is not shown:** the knob on in any production chain (PD-HD, PD-VD, SBND, uBooNE). Every event's fit moves in its last
  bits, so turning it on anywhere needs that chain's own revalidation.
