# SBND optical (PMT) detector simulation chain

How `sbndcode/OpDetSim/` turns propagated optical photons into digitized
`raw::OpDetWaveform`s, and exactly how the `PMTCoated*Eff` / `PMTUncoatedEff`
efficiencies in `digi_pmt_sbnd.fcl` are applied.

> Source tree referenced here: `sbndcode/sbndcode/OpDetSim/` (a linked checkout at
> `./sbndcode`). This is the **LArSoft/sbndcode** simulation — a different codebase from
> the WCT `QLMatching` model documented in [`sbnd-pmt-efficiencies.md`](sbnd-pmt-efficiencies.md).
> The two use the same *numbers* (0.0392 / 0.026 / 0.0357); for the uncoated PMT both now
> drop the direct/VUV term (see the consistency note under Q1).

## The efficiency parameters

```fcl
PMTCoatedVUVEff_tpc0: 0.0392   # coated PMT, direct (VUV 128 nm) light
PMTCoatedVISEff_tpc0: 0.026    # coated PMT, reflected (VIS ~430 nm) light
PMTUncoatedEff_tpc0:  0.0357   # uncoated PMT (single efficiency)
# ...identical _tpc1 block
```

PMTs are classified by `sbndPDMapAlg::pdType()` (from `sbnd_pds_mapping.json`) into
`pmt_coated` / `pmt_uncoated` (and `xarapuca_vuv` / `xarapuca_vis`). TPC is selected by
channel parity: `ch % 2 == 0` → tpc0, else tpc1 (`DigiPMTSBNDAlg.cc:244,302,319,376,437,458`).

## The user's three questions

### Q1 — Is `PMTUncoatedEff` only for the VIS (reflected) part?

**Yes.** An uncoated-PMT waveform is built **only** for reflected photons. The worker's
dispatch has a branch for coated PMTs (both Direct and Reflected) but the uncoated branch
is gated on `Reflected`:

```cpp
// opDetDigitizerWorker.cc:220
else if( (Reflected) && (pdtype == "pmt_uncoated") ) { // Uncoated PMT channels
  pmtDigitizer->ConstructWaveformLiteUncoatedPMT(ch, litesimphotons, ...);
}
```

There is **no `(Direct && pmt_uncoated)` branch** — direct (VUV) photons that land on an
uncoated channel are simply never digitized. So the single `PMTUncoatedEff` (0.0357)
multiplies only the reflected/VIS stream; the uncoated PMT has no VUV response at all in
simulation. This matches the physics: an uncoated PMT has no TPB wavelength-shifter, so it
cannot detect the 128 nm direct light, only the already-shifted reflected visible light.

> **Consistency with the WCT `QLMatching` model.** In
> [`sbnd-pmt-efficiencies.md`](sbnd-pmt-efficiencies.md) the uncoated PMT now carries
> `VUVEff = 0`, `VISEff = 0.0357` — its VUV term is zeroed, matching the simulation, which
> never digitizes direct-on-uncoated. (The original larwirecell port used
> `VUVEff = VISEff = 0.0357`, a nonzero VUV term; that has been corrected so predicted and
> simulated light agree on uncoated channels.)

### Q2 — Do VUV and VIS share one optical-photon collection, or separate ones (before efficiency)?

**Separate collections.** Direct and reflected light arrive as **two distinct art
products**, distinguished by the product instance name, not by a per-photon flag:

```cpp
// opDetDigitizerWorker.cc:154
const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
```

They are sorted into two separate maps and routed by PD type
(`opDetDigitizerWorker.cc:204-216`):

```cpp
if( pdtype == "pmt_coated" ){
  if(Reflected) ReflectedPhotonsMap.insert({ch, litesimphotons});
  else          DirectPhotonsMap.insert({ch, litesimphotons});
}
```

For a coated PMT both maps are handed to one constructor
(`opDetDigitizerWorker.cc:269`, `ConstructWaveformLiteCoatedPMT(ch, …, DirectPhotonsMap,
ReflectedPhotonsMap, …)`), and the **two efficiencies multiply the two streams
independently** inside `DigiPMTSBNDAlg.cc`:

```cpp
// Direct stream → VUV efficiency
mean_photons = directPhotons.second * _PMTCoatedVUVEff;     // :442
// Reflected stream → VIS efficiency
mean_photons = reflectedPhotons.second * _PMTCoatedVISEff;  // :463
```

So it is **not** one shared photon collection split per-photon into VUV/VIS — the split
already happened upstream (in the photon-library / `PDFastSim` stage), and each stream
keeps its own efficiency right through detection.

### Q3 — The optical-photon simulation chain, end to end

Module `opDetDigitizerSBND_module.cc` → per-channel `opDetDigitizerWorker.cc` →
`DigiPMTSBNDAlg` (PMTs) / `DigiArapucaSBNDAlg` (X-ARAPUCAs).

1. **Input.** `sim::SimPhotons` (full per-photon list) or `sim::SimPhotonsLite` (photons
   binned per 1 ns time bin); `UseSimPhotonsLite: true` is the default
   (`opdetdigitizer_sbnd.fcl`). Producer `PDFastSim`, two instances: default (`Direct`) and
   `Reflected`. Loaded in `produce()` via `getMany<…>()`.
2. **Routing.** Each channel's `pdType()` selects coated / uncoated / xarapuca_vuv /
   xarapuca_vis; Direct vs Reflected from the product instance name (Q2).
3. **Detection (the efficiency draw).** Two equivalent mechanisms depending on input type:
   - **SimPhotons (full):** per-photon Bernoulli trial —
     `if (fFlatGen.fire(1.0) < eff) { … }` (`DigiPMTSBNDAlg.cc:246, 306, 323`).
   - **SimPhotonsLite (default):** Poisson on the scaled count —
     `mean_photons = count * eff; accepted_photons = fPoissonQGen.fire(mean_photons)`
     (`DigiPMTSBNDAlg.cc:382, 443, 464`).
   - **Key nuance — `ScintPreScale`.** The fcl efficiencies are *not* used raw: at
     construction each is divided by `larProp->ScintPreScale()`
     (`DigiPMTSBNDAlg.cc:13-18`), which undoes the pre-scaling already applied when the
     photon library was generated, so the net detection probability is correct. A guard
     throws if any result exceeds 1.0 (`:34-46`).
4. **Photon timing.** For each accepted photon: transit-time spread (`Transittimespread`,
   a Gaussian of width `TTS`), TPB re-emission time `fTimeTPB->fire()` (sampled from a
   histogram), plus fixed `CableTime` / `TransitTime` offsets.
5. **Single-PE response (and PE nonlinearity).** `AddSPE()` adds the per-channel SPE
   template scaled by the photo-electron count; template chosen by `PMTSinglePEmodel`
   (`data` / `ideal` / `testbench`, loaded from `digi_pmt_sbnd_v2int0.root`), with optional
   gain fluctuations. **The PE count fed to `AddSPE` is first passed through the PMT
   nonlinearity (saturation) model** when enabled — see
   [§ PMT nonlinearity](#pmt-nonlinearity-saturation-npe_true--npe_observed) below.
6. **Noise.** Gaussian baseline/line noise (`AddLineNoise`, width `PMTBaselineRMS`),
   optional data-driven noise templates (`AddDataNoise`), and dark counts
   (`AddDarkNoise`, exponential inter-arrival at `PMTDarkNoiseRate`, each an SPE).
7. **Saturation + digitization.** Clamp to `PMTBaseline ± PMTADCDynamicRange`
   (`CreateSaturation`); convert charge→ADC via `PMTChargeToADC` (−25.97). Output is a
   16-bit ADC sample vector.
8. **Trigger / windowing.** `opDetSBNDTriggerAlg` finds threshold crossings and slices
   ~10 µs windows around them (`FindTriggerLocations` → `sliceWaveforms`).
9. **Output.** `produces<std::vector<raw::OpDetWaveform>>()` — the full waveforms plus a
   `"slicedWaveforms"` instance (the triggered windows), with trigger bookkeeping products.

> X-ARAPUCAs follow the analogous SiPM path in `DigiArapucaSBNDAlg` (`xarapuca_vuv` uses
> the direct map, `xarapuca_vis` the reflected map, each with its own efficiency).

## PMT nonlinearity (saturation): NPE_true → NPE_observed

This models the fact that a PMT's response stops being proportional to the light when many
photo-electrons pile up in a short time (space-charge / dynode saturation). **It is
distinct from ADC saturation** (step 7 / `CreateSaturation`, which merely clamps the final
waveform to the digitizer's `PMTBaseline ± PMTADCDynamicRange`). The nonlinearity acts
earlier, on the *PE count* before it becomes a pulse.

### Where and how it is applied

`DigiPMTSBNDAlg` accumulates accepted PE into a 1 ns-resolution vector `nPE_v` (size = 2×
the 2 ns waveform, `DigiPMTSBNDAlg.cc:243`). Then, per time bin, instead of injecting the
raw PE it injects a *rescaled* count:

```cpp
// DigiPMTSBNDAlg.cc:260-265 (same pattern at :337, :397, :480 — all four PMT branches)
if(fParams.SimulateNonLinearity){
  AddSPE(t, wave, ch, fPMTNonLinearityPtr->NObservedPE(ch, t, nPE_v));
} else {
  AddSPE(t, wave, ch, nPE_v[t]);
}
```

The tool is built only when the `NonLinearityParams` block is present
(`SimulateNonLinearity = nonLinearityParams.get_if_present(...)`, `DigiPMTSBNDAlg.cc:755`;
construction + `ConfigureNonLinearity()` at `:146-148`). **It is ON by default** —
`digi_pmt_sbnd.fcl:50-52` sets `NonLinearityParams: @local::PMTNonLinearityTF1ChannelByChannel`.
Comment that line out to disable.

### The functional form

The response is a ROOT `TF1`. The active form (`PMTAlg/pmtnonlinearity_config.fcl`) is

```
NObserved(x) = x / sqrt( 1 + (x/p0)^p1 )
```

where **`x` = the true PE accumulated over the preceding `AttenuationPreTime` window**
(`bin-4 .. bin` inclusive = `PreTime+1` = **5 one-ns samples**, spanning ~4 ns),
`p0` is a PE-saturation scale and `p1` a steepness exponent. For small `x` it is
≈ `x` (linear); for large `x` it bends over and grows only slowly. (`pmtnonlinearity_config.fcl`
also keeps a commented-out alternative form; the one above is what runs.)

At configure time a per-PE lookup of the *attenuation factor* is precomputed
(`PMTNonLinearityTF1ChannelByChannel_tool.cc:100-113`):

```cpp
fPEAttenuation_V[opch][pe] = TF1.Eval(pe) / pe;        // factor < 1, for pe in [range0,range1)
// pe < NonLinearRange[0]=100  -> factor stays 1   (linear)
// pe >= NonLinearRange[1]     -> fully saturated constant = round(TF1.Eval(range1))
```

and applied per bin (`:116-123`):

```cpp
npe_acc = sum(pe_vector[bin-PreTime .. bin]);          // 5-sample (PreTime+1) inclusive sum
if (npe_acc < NonLinearRange[1]) return pe_vector[bin] * fPEAttenuation_V[opch][npe_acc];
else                            return fPESaturationValue_V[opch];
```

So each bin's PE is scaled by the attenuation factor evaluated at the *local accumulated
load*. Note the consequence (important for the curve below): the saturation is driven by an
**instantaneous rate** — the running sum over a 5-sample window — **not by the pulse's grand
total**. If all PE of a pulse land within that window, the observed total equals
`TF1.Eval(total)`; if the same PE are spread over a longer time, each window saturates
independently and the net attenuation is *much milder*. This is exactly why the two curves
below (envelope vs realized) differ so dramatically.

### Per-PMT curves (your Q3 — yes)

Two tools implement the interface:

| tool | params | range | per-channel? |
|------|--------|-------|--------------|
| `PMTNonLinearityTF1` | fixed `[p0,p1]=[269, 1.84]` in fcl | `[100, 701]` | no — one curve for all |
| `PMTNonLinearityTF1ChannelByChannel` *(default)* | `p0,p1` per channel from PMT calibration DB | `[100, 7000]` | **yes** |

The default channel-by-channel tool sets, for each `opch`
(`PMTNonLinearityTF1ChannelByChannel_tool.cc:96-97`):

```cpp
SetParameter(0, fPMTCalibrationDatabaseService->getNonLineatiryPESat(opch));  // p0
SetParameter(1, fPMTCalibrationDatabaseService->getNonLineatiryAlpha(opch));  // p1
```

so **every PMT has its own saturation curve**, pulled from `PMTCalibrationDatabase`.

**Provenance of the per-channel `(PESat, Alpha)`.** They are not in the source tree — the
provider (`Calibration/PDSDatabaseInterface/PMTCalibrationDatabaseProvider.cxx`) queries a
remote LArSoft conditions DB via `lariov::DBFolder`: table **`pds_calibration`**, columns
**`nonlinearity_pesat`** / **`nonlinearity_alpha`**, tag **`v3r1`** (SBND Fall-2025
production, `calibration_database_PDS_TagSets_sbnd.fcl`). If a channel is absent the
provider returns `pesat = alpha = 0` (nonlinearity effectively off for it). Obtaining the
120 per-PMT values therefore requires a LArSoft session with DB access; export them to a
CSV (`opch,pesat,alpha`) to drive the standalone tool below.

### From observed PE to reconstructed PE (NPE_reco)

The nonlinearity is applied **on the waveform** (per bin, `AddSPE` injects `npe × SER`), so
there is **no direct `NPE_true → NPE_reco` function in sbndcode** — the total→total relation
only exists as an integral over the time profile. But the round-trip is an *identity*: reco
is **linear** in the recorded PE (`PE = OpHit.Area / SPEArea`, and the SER area is calibrated
to `SPEArea` per 1 PE — see [`sbnd-opdetreco-chain.md`](sbnd-opdetreco-chain.md)). Building
the waveform `Σ_t observed(t)·SER` and integrating it back therefore gives, exactly,

> **NPE_reco = Σ_t observed(t)** &nbsp; (in the unclipped regime)

i.e. summing the per-bin observed PE *is* the reconstructed PE; the deconvolution+area÷SPEArea
step adds nothing nonlinear. `pmt_nonlinearity_curve.py` (`roundtrip_reco_pe`) verifies this
numerically. The one separate effect is **ADC saturation** (`CreateSaturation`, a hard clip
at `PMTBaseline ± PMTADCDynamicRange`): with `digi_pmt_sbnd.fcl` values (baseline 14250,
dynamic range 14250) it needs a >14250-ADC excursion, but even the worst case — a single
burst capped at one SER × `round(Eval(7000))≈349` PE ≈ 8–9k ADC — never reaches it. So ADC
saturation is **inert** in the relevant range; the PE-accumulation model above is the whole
story.

### Scope: what the standalone curve reproduces vs. skips

The curve is **not** produced by running the full SBND sim+reco chain — there is no Geant4
optical sim, no `SimPhotons`, no SER-convolved digitizer waveform, and **no Wiener
deconvolution / OpHit finding**. It maps onto the three chain stages as:

| chain stage | real SBND | standalone tool |
|-------------|-----------|-----------------|
| 1. photon arrival (`nPE_v`) | Geant4/`PDFastSim` → `SimPhotons` → efficiency → per-1 ns `nPE_v` | **assumed** time profile (single-burst, or fast+slow scintillation mixture) of `N_true` detected PE |
| 2. nonlinearity (`NObservedPE`) | per-bin TF1 on the 5-sample running sum | **reproduced exactly** from the C++ |
| 3. waveform → reco PE | `AddSPE` (Σ observed·SER) → ADC → Wiener deconv → `Area/SPEArea` | **replaced by the identity** `NPE_reco = Σ_t observed(t)` (+ the minimal round-trip above) |

Skipping stage 3's machinery is exact, not an approximation: every step from per-bin observed
PE to reconstructed PE is **linear** (linear SER superposition, SER area calibrated to
`SPEArea`/PE, linear deconvolution + integration), and ADC saturation never engages — so the
full chain would reproduce this curve identically. The nonlinearity lives *entirely* in stage 2,
which is copied verbatim. The **only modelling assumption is the stage-1 time profile**
(`f_fast/tau_fast/tau_slow`): hence the *realized* curve is "assumed profile" while the
*envelope* (parameter-free TF1) is authoritative. Running the full chain would add only
the real photon timing, hit-finder noise/threshold scatter, and trigger windowing — none of
which shift the mean mapping.

### How to get the NPE_true ↔ NPE_observed curve

A standalone reproduction lives at
[`sbnd_xin/pmt_nonlinearity_curve.py`](../../sbnd_xin/pmt_nonlinearity_curve.py) — it
reproduces `NObservedPE` exactly (TF1 + 5-sample window + per-bin scaling + cap), plus the
explicit waveform round-trip above, and needs no LArSoft. Single-channel mode sweeps
`NPE_true → NPE_observed` for **two** profiles (defined next), plus the numeric inverse /
saturation correction; `--all-pmt --params-csv <file>` overlays one curve per PMT from a
per-channel `(opch,pesat,alpha)` CSV. The committed all-PMT plot
(`sbnd_xin/pics/pmt_nonlinearity_allpmt.png`) is the **realized (scintillation)** curve —
the physical NPE_true→NPE_reco for real flashes, not the worst-case envelope — generated
from the real `v3r1` per-channel values in `sbnd_xin/pmt_nonlin_params_v3r1.csv` (120 PMTs:
104 with a saturation curve, PESat 130–400 / Alpha 1.5–2.2; 16 with `pesat=alpha=0` →
nonlinearity off, treated as linear). It stays within ~1% of linear up to NPE≈1000 and
attenuates ~8–14% by NPE≈5000 (lower-PESat channels bend earliest). With no `--params-csv`
it falls back to a clearly-labelled *illustrative* synthetic spread.

**Envelope vs realized — why they look so different.** Because saturation acts on the
instantaneous 5-sample rate (not the total), the mapping depends entirely on the time profile:

- **Envelope (single burst):** all `N` PE land in one window → `npe_acc = N` →
  `observed = TF1.Eval(N)`. This is the **worst case** and is the bare TF1 shape; it bends at
  `N≈100` (range_lo) and plateaus near `round(Eval(range_hi))` (≈268 global / ≈349 per-channel).
  **It is *not* the physical NPE_total→NPE_reco map** — it is the unphysical all-light-in-one-bin
  limit.
- **Realized (scintillation):** PE follow a fast(singlet)+slow(triplet) mixture; the ~1.5 µs
  triplet spreads them over hundreds of windows, so the per-window rate stays below the 100-PE
  onset until `N` is large → the curve is **nearly linear**. This is the physical mapping in
  normal operation, and its headline is that *the nonlinearity barely engages*. Its one caveat:
  it depends on the assumed `f_fast/tau_fast/tau_slow` (the tool's defaults, not a single
  authoritative sbndcode constant), whereas the envelope (pure `Eval`) is authoritative.

What you call "NPE total vs NPE_nonlinear" is exactly this mapping. Two practical routes:

1. **The analytic per-channel curve (recommended, exact, no MC run)** = the envelope above.
   Plot `y = x / sqrt(1+(x/p0)^p1)` over your PE range, using that channel's `(p0,p1)`. Get
   `(p0,p1)` from the PMT calibration DB (`getNonLineatiryPESat` / `getNonLineatiryAlpha`),
   or use the global `[269, 1.84]` for a single representative curve.
   - **Forward** (true→observed): `y = TF1.Eval(x)`.
   - **Inverse** (observed→true, i.e. the saturation *correction*): the form is monotonic,
     so invert numerically (bisection) or build the inverse lookup from the same
     `TF1.Eval(pe)` table — it is literally the tool's `fPEAttenuation_V` table read the
     other way round, no need to re-derive it.

2. **The realized MC curve (pulse-shape-aware)** = the scintillation curve above. Because the
   simulation saturates on the 5-sample running sum, the *effective* total-PE relationship
   depends on the light's time profile. To measure it as it actually manifests:
   - Run `opdetdigitizer` **twice with the same RNG seed**, once with
     `NonLinearityParams` commented out (linear / NPE_true) and once with it on
     (NPE_observed), and scatter the two per-channel integrated PE (from OpHit/OpFlash, see
     [`sbnd-opdetreco-chain.md`](sbnd-opdetreco-chain.md)); **or**
   - Instrument `DigiPMTSBNDAlg` to dump, per channel per pulse, `sum(nPE_v)` (true) against
     `sum(NObservedPE)` (observed) — the cleanest ground-truth pair, no reconstruction in
     the loop.

   Route 2 captures the real spread around the analytic curve; route 1 is the upper-envelope
   (single-window burst) limit of it.

> Minor caveat: `NObservedPE` computes `start_bin = bin - AttenuationPreTime` with unsigned
> arithmetic, so the first `AttenuationPreTime` (4) bins of a waveform underflow the index;
> harmless in practice (no signal expected in the first 4 ns of the readout window) but
> worth knowing if you instrument it.

### Applying the non-linearity in QLMatching (predicted PE)

QLMatching predicts **true** per-PMT PE (`pred_flash[idet] += q·QtoL·(dir_vis·VUVEff +
ref_vis·VISEff)`, `match/src/QLMatching.cxx`), but compares against the **measured**
post-saturation `op_pes`. To bring the prediction into the same observed space, the realized
curve above is fitted per PMT to a **monotone log-quadratic capped power law** and applied to
each PMT's predicted total:

```
pred' = pred                                   for pred <= knee (700 PE)
      = knee · exp(beta·L + gamma·L²),  L = ln(pred/knee),   for pred > knee
```

(beta=1, gamma=0 ⇒ identity). The per-PMT (beta, gamma) are fitted by
`sbnd_xin/pmt_nonlinearity_curve.py --emit-qlmatching` to NPE_true=10⁵ (≤2% over the data
regime; see `sbnd_xin/pics/pmt_nonlinearity_fit.png`) and written to
`cfg/pgrapher/experiment/sbnd/pmt_nonlinearity_params.jsonnet`. QLMatching reads
`pmt_nonlinearity` / `pmt_nl_knee` / `pmt_nl_beta` / `pmt_nl_gamma`. **As of 2026-06-04 the SBND
default is ON**: the canonical `qlmatching.jsonnet` `matching(...,pmt_nl=true,...)` /
`matching_joint(...,pmt_nl=true,...)` bake the `nl_on` overlay by default (so SBND LArSoft
production also gets it); **`pmt_nl=false` reproduces the old OFF/identity baseline byte-identical**
(verified with `wcsonnet`: no `pmt_nl*` keys emitted). The sbnd_xin `qlmatching.jsonnet` wrapper
just forwards `pmt_nl`; `run_ql_evt.sh` and `run_clust_QL_evt.sh` default `PMT_NL=true`. Non-SBND
detectors are unaffected (SBND-only file).

**OFF-vs-ON validation** (`sbnd_xin/ql_nonlin_compare.py`,
`sbnd_xin/pics/ql_pmt_nonlin_compare.png`; 10-event MC + data; matching unbroken, 106/106 and
101/101 common matches, ON pred ≤ OFF by construction):
- **MC**: detrended pred/meas *rises* with brightness (1.27→1.31 vs the <700 baseline) — a
  genuine mild-saturation signature — and the correction **flattens it toward 1** in every
  bright bin (→1.25→1.20). The correction acts in the right direction.
- **Data**: detrended pred/meas *falls* with brightness (0.95→0.62) — the brightest PMTs
  measure **more light than the reconstructed charge predicts** (a charge/light effect, e.g.
  recombination or near-PMT visibility, **not** PMT saturation), so lowering the prediction
  cannot help and slightly steepens the slope.
- The per-channel effect is small relative to the pred/meas scatter (|median error| ~0.1–0.5).

So the correction is physically sound and verified. On these samples its effect is small (it
only engages per-PMT above the 700 PE knee, rare in practice — see the PE-error study below,
where the effective light error moves only ~1 pt with NL on vs off). It is now shipped
**default-ON for SBND** by request; `pmt_nl=false` recovers the old behavior. A real
quantitative benefit would still require re-tuning the efficiencies / `QtoL` in tandem (a
separate study), since the current tuning already absorbed the average response.

### Per-PMT light error and the data/MC normalization (PE-error study)

`sbnd_xin/ql_pe_error.py` (see `sbnd_xin/docs/17_pe-error-study.md`) uses the 10 data + 10 MC
hand-scans to fit the per-PMT light error `a` from `E[(pred−meas)²]=meas+a·pred²` (the matcher
assumes 30% ⇒ `a=0.09`). Findings: the robust per-PMT error is **~30% for data** (≈33% on clean
matches) and **tighter for MC** (~15–24%). The notable result is a **data/MC normalization
mismatch of ~30%**: the same photon model **over**-predicts data light (median meas/pred≈0.86)
but **under**-predicts MC (≈1.14).

**The data recipe (implemented, SBND default):** two SBND-config corrections fix this for data
(sim is left alone):
1. **DATA prediction scale `QtoL = 0.86`** (sim 1.0) — since `pred = q·QtoL·vis·eff`, lowering
   the data `QtoL` removes the ~16% over-prediction.
2. **PE-dependent error `σ² = meas + max(5 PE, 0.25·pred)²`** (`pe_err_floor=5`, `pe_err_frac=0.25`,
   `pe_err_knee=20`, both data and sim) — a constant floor at low PE, 25% fractional at high PE.
   `pe_err_on_pred=true` makes the **bundle χ²** use the *predicted* PE for this error
   (`TimingTPCBundle::examine_bundle`); the **LASSO** keeps the per-flash measured-based weight
   (its measured vector is shared across candidate bundles).

After the recipe (`sbnd_xin/ql_recipe_compare.py`, hand-scan kept matches): median **χ²/ndf drops
from ~5 → 1.72 (data) / 1.04 (MC)**, **data meas/pred → 1.01** (MC stays 1.17, unscaled), and
**ks is ~unchanged** (it is shape-normalized, so scaling `pred` leaves it invariant). These
recalibrated χ²/ks plus the hand-scan labels are the inputs for further QLMatching cut tuning.

## Where this lives in code

- Module / worker: `OpDetSim/opDetDigitizerSBND_module.cc`, `opDetDigitizerWorker.cc`.
- PMT digitizer: `OpDetSim/DigiPMTSBNDAlg.cc` / `.hh`.
- Config: `OpDetSim/digi_pmt_sbnd.fcl`, `opdetdigitizer_sbnd.fcl`.
- PD-type map: `OpDetSim/sbndPDMapAlg_tool.cc`, `sbnd_pds_mapping.json`.
- Nonlinearity: `OpDetSim/PMTAlg/PMTNonLinearity.hh` (interface),
  `PMTNonLinearityTF1_tool.cc`, `PMTNonLinearityTF1ChannelByChannel_tool.cc`,
  `pmtnonlinearity_config.fcl`.
- Trigger: `OpDetSim/opDetSBNDTriggerAlg.cc` / `.hh`.
