# Charge-light matching port (`WireCellMatch`)

## What this subpackage does

`wire-cell-toolkit/match/` is a wire-cell-toolkit subpackage that ports the
charge-light matching code originally living in
`larwirecell/larwirecell/qlmatch/` (`QLMatching`) and its larsoft
dependency `larsim/PhotonPropagation/SemiAnalyticalModel` into WCT with no
larsoft / `art` / `fhicl` runtime dependency.

After this port, `wire-cell` can run the SBND `QLMatching` step from a
standalone jsonnet on already-clustered inputs + an opflash tensor, with
all geometry and parameterisation data supplied as JSON loaded via
`WireCellUtil/Persist`.

## Layout

```
match/
├── inc/WireCellMatch/
│   ├── SemiAnalyticalModel.h   # WCT port of larsim phot::SemiAnalyticalModel
│   ├── Opflash.h               # moved from larwirecell, identical interface
│   ├── TimingTPCBundle.h       # moved from larwirecell, identical interface
│   ├── QLMatching.h            # ITensorSetFanin component, JSON-driven config
│   └── Util.h                  # BEE-JSON dump helpers
├── src/                        # implementations
└── wscript_build               # WireCellMatch waf target
```

`bld.smplpkg('WireCellMatch', use='WireCellClus WireCellAux WireCellIface WireCellUtil')`.

## What was ported

### SemiAnalyticalModel (minimal scope)

Ported from `larsim/PhotonPropagation/SemiAnalyticalModel.{h,cxx}` —
covers only the code paths that QLMatching exercises in SBND:

- **dome PMTs (type=1)** + **flat (X)Arapucas (type=0)** at the
  anode/cathode plane (orientation=0).
- VUV direct visibilities and VIS reflected visibilities, with
  Gaisser-Hillas correction + border corrections + Omega-dome /
  rectangle / disk solid-angle helpers.

Not ported (SBND does not use them; add back if a future detector needs them):

- lateral PD corrections (`FlatPDCorrLat`, lateral VUV/VIS tables)
- anode reflections (`IncludeAnodeReflections`)
- Xe absorption
- field-cage transparency
- vertical-border correction mode

### larsoft dependencies replaced

| larsoft surface                              | WCT replacement                                                  |
|----------------------------------------------|------------------------------------------------------------------|
| `geo::Point_t`, `geo::Vector_t`              | `WireCell::Point` (`D3Vector<double>`)                          |
| `geo::GeometryCore` (OpDet array)            | `std::vector<SemiAnalyticalModel::OpticalDetector>` from JSON   |
| `geo::GeometryCore` (active volume / cathode)| `SemiAnalyticalModel::Geometry` struct from JSON                |
| `larg4::ISTPC::extractActiveLArVolume`       | numeric `active_center_y/z`, `active_size_y/z`                  |
| `detinfo::LArPropertiesService::AbsLengthSpectrum` | scalar `vuv_absorption_length` in JSON (cm)                |
| `lar::providerFrom<>` service handles        | nothing — geometry is plain data                                 |
| `phot::OpticalPath` tool + `SBNDOpticalPath` | inlined same-TPC visibility check (same X-sign)                  |
| `fhicl::ParameterSet`                        | `WireCell::Configuration` (`Json::Value`) + `Persist::load`     |
| `cet::exception`                             | `raise<ValueError>` / `WireCell::Exceptions`                     |
| `mf::LogInfo`                                | `Aux::Logger` / spdlog                                           |
| `art::make_tool`                             | not used; `OpticalPath` is inlined                              |
| `CLHEP::pi`                                  | `M_PI` (`units::pi` is `static const`, not constexpr)            |

### QLMatching component

The component itself was already a `WIRECELL_FACTORY` (`ITensorSetFanin`,
`IConfigurable`) — only the bits coupling it to larsoft needed to be cut:

- `fhicl::ParameterSet::make(...)` + `art::make_tool` replaced with a
  single `WireCell::Persist::load("semimodel_file")` at `configure()`
  time. The JSON's top-level keys (`VUVHits`, `VISHits`, `Geometry`,
  `OpDets`) are passed straight to the new `SemiAnalyticalModel` ctor.
- All references to `phot::SemiAnalyticalModel` updated to
  `WireCell::Match::SemiAnalyticalModel`.
- Default 312-entry `VUVEfficiency` / `VISEfficiency` arrays kept in
  the header so the standalone jsonnet does not need to ship them
  (overridable via `cfg["VUVEfficiency"]` / `cfg["VISEfficiency"]`).

## SBND data file

The component reads a single JSON at configure time:

```
cfg["semimodel_file"] = "sbnd/photodet/semi-analytical-sbnd.json"
```

This file lives in
`/exp/sbnd/app/users/yuhw/wire-cell-data/sbnd/photodet/semi-analytical-sbnd.json`
and is discovered via `WIRECELL_PATH`. Schema:

```json
{
  "VUVHits": { ... copied verbatim from semimodel_sbnd.fcl ... },
  "VISHits": { ... copied verbatim from semimodel_sbnd.fcl ... },
  "Geometry": {
    "active_center_y": 0.0,
    "active_center_z": 250.5,
    "active_size_y":   407.465,
    "active_size_z":   501.0,
    "cathode_x":       0.0,
    "vuv_absorption_length": 85.0
  },
  "OpDets": [
    {"x": ..., "y": ..., "z": ..., "h": ..., "w": ...,
     "type": 0|1, "orientation": 0},
    ...   // 312 entries for SBND v02_06
  ]
}
```

### How the JSON is generated

Tools live in
`wcp-porting-img/sbnd/standalone-sample/build-semi-analytical-data/`:

1. **`SBNDOpDetDumper_module.cc`** — one-off larsoft analyzer. Lives in
   `larwirecell/larwirecell/qlmatch/`, is built by the same
   `CMakeLists.txt` as `WireCellQLMatch`, and emits CSV `OPDET:` lines +
   `GEOM:` lines for the active volume / cathode / drift constants.
   Mirrors what `SemiAnalyticalModel::opticalDetectors()` and the
   constructor do at runtime, so the WCT-side JSON matches.
2. **`dump_sbnd_opdets.fcl`** — drives the analyzer for one empty event.
3. **`build_semi_analytical_sbnd_json.py`** — small recursive-descent
   parser for `fhicl-dump` output that extracts `VUVHits` / `VISHits`,
   plus a CSV reader for the OpDets, and writes the combined JSON.

End-to-end: see the "How to regenerate" section below.

## How to build and run

### Build (SL7 + sbndcode env)

```bash
in-gpvm-sl7.sh bash -c '
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
  setup sbndcode v10_14_02_03 -q e26:prof
  cd /exp/sbnd/app/users/yuhw/wire-cell-toolkit
  ./wcb configure --prefix=/exp/sbnd/app/users/yuhw/opt ...
  ./wcb -p --notests build install
'
```

### Run the standalone test

```bash
source /exp/sbnd/app/users/yuhw/wcp-porting-img/sbnd/setup-local-opt.sh
cd /exp/sbnd/app/users/yuhw/wcp-porting-img/sbnd/standalone-sample
wire-cell -l stdout -L info \
  -V reality=sim -V DL=6.2 -V DT=9.8 -V lifetime=6 -V driftSpeed=1.565 \
  -V input=. \
  -V semimodel_file=sbnd/photodet/semi-analytical-sbnd.json \
  -c wct-clus-matching-standalone.jsonnet
```

Inputs expected in cwd: `icluster-apa{0,1}-{active,masked}.npz` plus
`opflash_apa{0,1}.tar.gz`. Outputs `data-sep/0/0-{img,op}-apa{0,1}.json`
which `bee-upload.sh` then merges with the cluster JSONs.

### How to regenerate `semi-analytical-sbnd.json`

```bash
# 1. build the dumper inside larsoft (one-off; only needed when geometry
#    or correction tables change)
in-gpvm-sl7.sh bash -c '
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
  setup sbndcode v10_14_02_03 -q e26:prof
  source $LARSOFT_LOCAL/setup
  mrbsetenv
  cd $MRB_BUILDDIR && mrb b -j4 && mrb i -j4
'

# 2. dump OpDets + geometry
in-gpvm-sl7.sh bash -c '
  source $LARSOFT_LOCAL/setup; mrbslp
  lar -c $LARWIRECELL/larwirecell/qlmatch/dump_sbnd_opdets.fcl
' 2>&1 > dump.log
grep "^OPDET:" dump.log | sed "s/^OPDET://" > sbnd_opdets.csv
grep "^GEOM:"  dump.log > sbnd_geom.txt

# 3. parse FCL dump + opdets CSV into one JSON
./build_semi_analytical_sbnd_json.py \
  --fcl   /path/to/semimodel_sbnd-dump.fcl \
  --opdets sbnd_opdets.csv \
  --active-center-y 0 --active-center-z 250.5 \
  --active-size-y 407.465 --active-size-z 501 \
  --cathode-x 0 --vuv-absorption-length 85 \
  --out /exp/sbnd/app/users/yuhw/wire-cell-data/sbnd/photodet/semi-analytical-sbnd.json
```

## Verified

- **Build**: `WireCellMatch` builds cleanly with the WCT-master `wcb`
  toolchain under sbndcode `v10_14_02_03 / e26:prof`. Installed
  `libWireCellMatch.so` to `/exp/sbnd/app/users/yuhw/opt/lib/`.
- **Factory**: `wire-cell` loads `WireCellMatch`, `Factory::find_tn` for
  `QLMatching` resolves, `configure()` reads the JSON and instantiates
  the `SemiAnalyticalModel`.
- **End-to-end**: ran `wct-clus-matching-standalone.jsonnet` on the
  reference cluster + opflash inputs in
  `wcp-porting-img/sbnd/standalone-sample/`, produced
  `data-sep/0/0-{img,op}-apa{0,1}.json`, packaged with `bee-upload.sh`
  and uploaded to BEE for visual inspection.
- **Numerical agreement with legacy `lar` job**: per-(flash, cluster)
  predicted PE agrees to ~0.01 PE / ~5e-5 relative on every common
  bundle (APA1: bit-perfect on all 6 commons; APA0: 10/10 commons
  match, one bundle off by 5% downstream of `add_bundle` merge
  logic). Which marginal bundles survive the LASSO 0.05 threshold
  differs slightly because the threshold is just above the smallest
  pred-PE rows, so any 1e-5 jitter flips inclusion — that's pre-
  existing in the algorithm, not introduced by the port.

## Notes / gotchas

- **`vuv_absorption_length`**: this is *not* a guess. Take it from
  `lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum()`
  interpolated at 9.7 eV (Ar VUV peak). For SBND `v10_14_02_03` /
  `standard_properties` the value is **2000 cm**. The first pass of
  this port used a 85 cm placeholder and shipped half the predicted
  light — visible only when comparing to the legacy `lar` job. The
  `SBNDOpDetDumper` analyzer now emits the correct value into
  `sbnd_geom.txt` so it can't drift.
- **Boost 1.82 + gcc 12 + `-Werror=deprecated-declarations`**: the
  ellint / multi_array headers transitively include
  `boost/functional.hpp` which uses the deprecated
  `std::unary_function` / `binary_function`. Each `#include` that pulls
  these in is wrapped in `#pragma GCC diagnostic push / pop` in
  `SemiAnalyticalModel.cxx`, `Opflash.cxx`, `Util.cxx`.
- **`units::pi` is not `constexpr`**: WCT's `WireCellUtil/Units.h`
  defines `pi` as `static const double`, not `constexpr`. Use `M_PI`
  in `constexpr` contexts.
- **Cathode at X=0**: SBND's two-TPC layout has the cathode at X=0
  (shared between TPCs). `cathode_x = |0| = 0` in the JSON; the sign
  flip for the VIS hotspot still works because the model evaluates
  `plane_depth = scintPoint.x() < 0 ? -plane_depth : plane_depth`.
- **Same-TPC visibility**: replicated by an inline X-sign check; the
  SBND `phot::SBNDOpticalPath` tool was a 5-line "X sign match", so a
  separate pluggable interface was overkill for this scope.

## What's still to do (out of scope of this PR)

- Re-wire the larwirecell-side `QLMatching` (still living under
  `larwirecell/qlmatch/`) to use the WCT-side implementation, so the
  original `wct-clus-matching.fcl` larsoft job picks up the new code
  too. Currently the larwirecell tree still builds the legacy
  larsoft-dependent component; standalone wire-cell uses the new one.
- Port the remaining `SemiAnalyticalModel` branches if a future
  detector (DUNE-VD, ICARUS) needs them.
