# Drifting & Transport -- Code Examination

## Summary of Findings

| File | Bugs | Efficiency | Notes |
|------|------|-----------|-------|
| Drifter | 2 confirmed (comparator), 1 minor | 1 | Core drift physics component |
| DepoSetDrifter | 1 (division-by-zero on empty input) | 1 | Thin adapter |
| TrackDepos | 1 (charge sign logic) | 1 | Test/config source |
| BlipSource | 1 (memory management) | 0 | Blip depo source |
| RecombinationModels | 2 (division-by-zero) | 0 | Physics models |
| TransportedDepo | 1 (division-by-zero) | 0 | Header-only transport |

Total: 8 actionable findings.

---

## Drifter.h / Drifter.cxx

### Algorithm Overview

The Drifter takes individual depositions (`IDepo`) ordered by absolute time, assigns each to one of several configured X-regions (drift volumes), and drifts it to the response plane. Each X-region has three planes: cathode, response, and anode. Depos between cathode and response ("bulk") are drifted forward with full physics (absorption, diffusion). Depos between response and anode ("near") are anti-drifted backward in time with no physics, as a best-effort correction.

The Drifter buffers drifted depos per X-region in a time-sorted `std::set`, flushing "ripe" depos (those with drift time < current input time) to output as new inputs arrive. On EOS, all remaining depos are flushed.

### Bug 1 (confirmed): Comparator compares `lhs` to itself

**File:** `gen/src/Drifter.cxx`, lines 23-32

```cpp
bool Gen::Drifter::DepoTimeCompare::operator()(const IDepo::pointer& lhs, const IDepo::pointer& rhs) const
{
    if (lhs->time() == rhs->time()) {
        if (lhs->pos().x() == lhs->pos().x()) {   // BUG: compares lhs to lhs
            return lhs.get() < rhs.get();
        }
        return lhs->pos().x() < lhs->pos().x();   // BUG: compares lhs to lhs
    }
    return lhs->time() < rhs->time();
}
```

**Line 26:** `lhs->pos().x() == lhs->pos().x()` -- compares `lhs` X position to itself. Should be `rhs->pos().x()`. This condition is always true (barring NaN), so the tie-break by X position on line 29 is dead code.

**Line 29:** `lhs->pos().x() < lhs->pos().x()` -- also compares `lhs` to itself. This is always false and never executes (due to the bug on line 26), but if the line-26 bug were fixed, this line would also need fixing to `rhs->pos().x()`.

**Impact:** The `std::set` ordering falls through to pointer comparison when times are equal, which still provides a strict weak ordering. The set will not crash or corrupt, but the ordering is not by X position as clearly intended. This means depos with identical times are ordered by memory address rather than spatial position. In practice the physics output is likely unaffected since the ordering of same-time depos through `flush_ripe` does not change the physics, but it makes the output non-deterministic across runs with different allocation patterns.

### Bug 2 (minor): `flush_ripe` flushes depos with time < now, not <= now

**File:** `gen/src/Drifter.cxx`, lines 212-217

```cpp
while (depoit != depo_end) {
    if ((*depoit)->time() < now) {
        ++depoit;
        continue;
    }
    break;
}
```

The set is ordered by time, and this loop walks forward while `time < now`. Depos with `time == now` are NOT flushed. This is likely intentional (conservative -- only flush depos strictly older than the current input time) but worth noting in case the intent was `<=`.

### Efficiency 1: Re-sorting after merge from already-sorted sets

**File:** `gen/src/Drifter.cxx`, lines 191-198 (`flush`) and line 229 (`flush_ripe`)

Both `flush` and `flush_ripe` insert depos from multiple X-regions into `outq` and then call `std::sort`. Since each X-region's depo set is already time-sorted, an N-way merge would be O(N log K) instead of O(N log N). For small K (number of X-regions, typically 1-2), the difference is negligible.

### Key Algorithmic Details

- **Absorption:** `absorbprob = 1 - exp(-dt / m_lifetime)` (line 163). Charge is reduced by binomial sampling of absorbed electrons when `m_fluctuate` is true.
- **Diffusion:** Longitudinal and transverse extents grow as `sqrt(2*D*dt + sigma^2)` (lines 178-179), correctly adding in quadrature with pre-existing extent.
- **Anti-drift:** Depos in the "near" region (between anode and response) are moved backward in time with `direction = -1.0` (line 140), with no diffusion or absorption applied.
- **Scale factor:** A configurable `charge_scale` is applied to all drifted charge (line 183).

---

## DepoSetDrifter.h / DepoSetDrifter.cxx

### Algorithm Overview

A thin adapter that wraps an `IDrifter` (per-depo drifter) to implement the `IDepoSetFilter` interface. It takes a `DepoSet`, iterates over its depos (plus an appended EOS marker to flush the inner drifter), collects all output, removes the trailing EOS nullptr, and wraps the result in a new `SimpleDepoSet`.

### Bug 1: Division by zero when `charge_in == 0`

**File:** `gen/src/DepoSetDrifter.cxx`, line 67

```cpp
log->debug("call={} drifted ndepos={} Qout={} ({}%)", m_count, all_depos.size(), charge_out, 100.0*charge_out/charge_in);
```

If all input depos have zero charge (or the input set is empty), `charge_in` is 0.0, producing a floating-point division by zero (NaN or Inf in the log message). This is a minor issue since it only affects logging, but it could produce confusing debug output.

### Efficiency 1: Redundant sorting

As noted in the file's own header comment (lines 4-5), the per-depo drifter maintains time ordering internally, which is unnecessary work when the DepoSetDrifter could do a single sort at the end. The comment acknowledges this.

### Edge Case: EOS propagation

On receiving EOS (line 38), the DepoSetDrifter returns `out = nullptr` and `true` but does NOT forward the EOS to the inner drifter. If the inner drifter had been fed depos from a previous call and the framework sends an EOS, those buffered depos would be lost. However, in practice each `DepoSet` call includes its own EOS flush (line 46), so this is safe as long as the DepoSetDrifter is only used with complete DepoSets.

---

## TrackDepos.h / TrackDepos.cxx

### Algorithm Overview

A configurable source of depositions that generates depos along straight-line tracks. Each track is defined by a time, a ray (start/end points), and a charge parameter. Depos are placed at regular step intervals along the track, with time advancing at a fraction of the speed of light. An optional `group_time` parameter chunks the output into time-delimited groups separated by EOS markers.

### Bug 1: Ambiguous charge sign handling

**File:** `gen/src/TrackDepos.cxx`, lines 85-91

```cpp
double charge_per_depo = units::eplus;  // charge of one positron
if (charge > 0) {
    charge_per_depo = -charge / (length / m_stepsize);
}
else if (charge <= 0) {
    charge_per_depo = charge;
}
```

When `charge > 0`, it is interpreted as total charge distributed across all steps, and negated (electrons are negative). When `charge <= 0`, the raw value is used per depo. But the default initialization `units::eplus` (line 85) is immediately overwritten by one of the two branches (since `charge <= 0` covers the `charge == 0` case). The initial value of `units::eplus` is dead code.

When `charge == 0`, this produces depos with zero charge, which the Drifter will silently drop (Drifter line 121-124). This is arguably correct but undocumented.

When `charge < 0` (the default `-1.0`), the absolute value is used as electrons per depo. The header comment says "if < 0 then it gives the (negative of) absolute amount of charge per deposition" -- so `-1.0` means 1 electron per depo. This is consistent but the negation convention is confusing.

### Efficiency 1: Re-sorting after each track addition

**File:** `gen/src/TrackDepos.cxx`, line 102

```cpp
std::sort(m_depos.begin(), m_depos.end(), ascending_time);
```

Called after every `add_track`. If multiple tracks are added in `configure()`, the sort runs after each one. A single sort after all tracks are added would suffice. Since track count is typically small, this is low impact.

### Key Algorithmic Details

- Track depos are evenly spaced at `m_stepsize` intervals (default 1mm).
- Time of each depo: `time + step / (clight_fraction * c)`.
- The `group_time` feature inserts EOS markers to delimit time groups, useful for chunked processing.

---

## BlipSource.h / BlipSource.cxx

### Algorithm Overview

Generates random point-like depositions ("blips") in a configurable volume, with configurable charge distribution and time distribution. Supports mono-energetic charge, PDF-sampled charge, and exponentially distributed inter-arrival times (modeling radioactive decay, e.g., Ar39).

### Bug 1: Raw `new`/`delete` for polymorphic members

**File:** `gen/src/BlipSource.cxx`, lines 152, 155, 166, 177 (new) and lines 26-31 (delete)

The `m_ene`, `m_tim`, and `m_pos` members are raw pointers allocated with `new` and freed in the destructor with `delete`. If `configure()` is called twice, the old allocations are leaked because there is no `delete` before the second `new`.

```cpp
// In configure():
m_ene = new ReturnValue(ene["value"].asDouble());  // potential leak if m_ene already set
```

Should either use `std::unique_ptr` or explicitly delete before reassigning.

### Edge Case: `m_stop` not initialized in constructor

**File:** `gen/inc/WireCellGen/BlipSource.h`, line 42

```cpp
double m_time, m_stop;
```

`m_stop` is not initialized in the constructor (only `m_time` is set to 0.0). It gets set in `configure()` from `tim["stop"]`. If `operator()` is called before `configure()`, `m_stop` is uninitialized, leading to undefined behavior. In practice, the framework always calls `configure()` first, so this is low risk.

### Edge Case: Pdf division by zero

**File:** `gen/src/BlipSource.cxx`, line 117

```cpp
double rel = (cp - cumulative[ind - 1]) / (cumulative[ind] - cumulative[ind - 1]);
```

If two consecutive cumulative values are equal (i.e., a PDF bin has zero probability), the denominator is zero. This could happen if the user supplies a PDF with zero-valued bins. In practice `rng->uniform` is unlikely to land exactly on such a boundary, but it is theoretically possible.

### Key Algorithmic Details

- Inter-arrival times follow exponential distribution with rate = activity (Bq).
- The PDF sampler uses inverse CDF with linear interpolation within bins.
- Blip positions are uniform within a box.

---

## RecombinationModels.h / RecombinationModels.cxx

### Algorithm Overview

Three recombination models converting energy deposition (dE) to ionization charge (dQ):

1. **MipRecombination:** `dQ = (Rmip / Wi) * dE` -- constant recombination factor, independent of dE/dX.
2. **BirksRecombination:** `R = A3t / (1 + (dE/dX) * k3t / (Efield * rho))`, then `dQ = R * dE / Wi`.
3. **BoxRecombination (Modified Box):** `R = ln(A + tmp) / tmp` where `tmp = (dE/dX) * B / (Efield * rho)`, then `dQ = R * dE / Wi`.

Each model also provides an inverse `dE(dQ, dX)` method.

### Bug 1: Division by zero when dX == 0 in Birks model

**File:** `gen/src/RecombinationModels.cxx`, line 54

```cpp
const double R = m_a3t / (1 + (dE * units::cm / dX) * m_k3t / (m_efield * m_rho));
```

When `dX == 0`, this produces division by zero. The function signature has `dX = 0.0` as default (inherited from the interface). The MipRecombination model ignores dX, but Birks and Box models require non-zero dX.

### Bug 2: Division by zero when dX == 0 in Box model

**File:** `gen/src/RecombinationModels.cxx`, line 97

```cpp
const double tmp = (dE /units::MeV*units::cm/ dX) * m_b / (m_efield * m_rho);
const double R = std::log(m_a + tmp) / tmp;
```

Same issue. Also, when `tmp == 0` (which happens when `dE == 0`), the `log(A)/0` division produces infinity. The limit of `ln(A+x)/x` as `x->0` is `1/A` for small x, but the code does not handle this.

### Edge Case: Box model `dE()` inverse can produce negative dE

**File:** `gen/src/RecombinationModels.cxx`, line 104

```cpp
const double a_exp = std::exp(dQ/dX*units::cm * coeff * m_wi);
const double numerator = (a_exp - m_a) * units::MeV/units::cm *dX;
```

If the exponent is small enough that `a_exp < m_a`, the numerator is negative, producing a negative dE. This can happen for very small dQ/dX values. Whether this is physically meaningful depends on the use case.

### Key Algorithmic Details

- Wi = 23.6 eV per electron (work function of argon).
- Default Rmip = 0.7 (70% recombination for MIP).
- Birks parameters: A3t = 0.8, k3t = 0.0486 g/(MeV cm^2) kV/cm.
- Box parameters: A = 0.930, B = 0.212 g/(MeV cm^2) kV/cm.
- Default E-field: 500 V/cm, LAr density: 1.396 g/cm^3.
- All three models use the BNL LAr properties reference: http://lar.bnl.gov/properties/pass.html#recombination

---

## TransportedDepo.h (header only)

### Algorithm Overview

A simple depo wrapper that transports a depo from its current X position to a target `location` at a given `velocity`. The new time is calculated as `t_new = t_old + (x_old - location) / velocity`. Charge, id, pdg, and energy are forwarded from the original depo. The original depo is stored as `prior()`.

### Bug 1: Division by zero when velocity == 0

**File:** `gen/inc/WireCellGen/TransportedDepo.h`, line 23

```cpp
m_time = from->time() + dx / velocity;
```

If `velocity == 0`, this produces floating-point division by zero. There is no guard. In practice, callers should always provide a non-zero drift velocity, but the constructor does not validate this.

### Design Note: No diffusion extents

`TransportedDepo` does not override `extent_long()` or `extent_tran()`, which default to 0.0 from `IDepo`. This is intentional -- this class is used for simple transport without diffusion physics, as documented in the Drifter's "near region" anti-drift path.

---

## Cross-cutting Observations

1. **Comparator bug in Drifter is the most impactful finding.** The `DepoTimeCompare` comparing `lhs` to itself on lines 26 and 29 is clearly a copy-paste error. While the `std::set` still functions (the pointer-based fallback provides strict weak ordering), the spatial ordering intent is defeated.

2. **Division-by-zero in recombination models** is a design issue: the interface declares `dX = 0.0` as a default parameter, but two of three implementations require `dX != 0`. This is a latent bug that would manifest if a caller uses the default.

3. **BlipSource memory management** uses raw pointers where `std::unique_ptr` would be safer and would prevent leaks on re-configuration.
