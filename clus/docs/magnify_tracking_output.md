# Magnify Tracking Output: ROOT File Format and Conversion

**Source files:**
- `root/apps/wire-cell-uboone-magnify-tracking-convert.cxx` — the conversion tool
- `qlport/uboone-mabc.jsonnet` — `tracking_visitor` configuration (line ~1176+)

**External viewer:** https://github.com/BNLIF/Magnify-tracking

---

## End-to-End Data Flow

```
TaggerCheckNeutrino::visit()
        │  PR graph with fitted segments, kinematics
        ▼
[tracking_visitor]  (UbooneMagnifyTrackingVisitor)
        │  writes ROOT file per event
        ▼
tracking_<run>_<subrun>_<event>.root  (or track_com_*.root)
        │
        ▼
wire-cell-uboone-magnify-tracking-convert  -b input.root -o output.root -f2
        │
        ▼
temp.root   →  opened by Magnify-tracking ROOT macro viewer
```

---

## The Tracking Visitor in jsonnet

In `uboone-mabc.jsonnet`, the `tracking_visitor` is configured at the end of
the pipeline (after `tagger_check_neutrino`):

```jsonnet
local tracking_visitor = {
    type: "UbooneMagnifyTrackingVisitor",
    data: { ... }
};

// Appended to cm_pipeline only when tracking_output != "":
] + (if tracking_output != "" then [tracking_visitor] else []);
```

`tracking_output` is a parameter passed to the jsonnet (e.g., set on the
command line via `-V tracking_output=tracking`). When empty, no tracking
file is produced, skipping this potentially slow step.

---

## Input ROOT File: `T_rec_charge` Tree

The tracking visitor writes one flat-ntuple entry per fitted point on each
segment. Branches:

| Branch | Type | Description |
|---|---|---|
| `x`, `y`, `z` | `Double_t` | 3D position of the fitted point (cm) |
| `q` | `Double_t` | Charge at this point (`dQ`, electrons) |
| `nq` | `Double_t` | Number of hits ("`dx`" in the conversion code, hit count) |
| `ndf` | `Double_t` | Cluster ID (used to group points per track) |
| `pu`, `pv`, `pw`, `pt` | `Double_t` | Projected charges on U, V, W planes and time plane |
| `reduced_chi2` | `Double_t` | Track fit quality (chi² / ndf) |
| `flag_vertex` | `Int_t` | 1 if this point is at the primary vertex |
| `sub_cluster_id` | `Int_t` | Sub-cluster index within the cluster |

Other trees cloned from the input:
- `T_bad_ch`: dead channel list (copied as-is to the output)
- `T_proj`, `T_proj_data`: 2D wire projection data
- `Trun`: run-level calibration constants including `dQdx_scale`, `dQdx_offset`

---

## The Conversion Tool

### Command

```bash
wire-cell-uboone-magnify-tracking-convert \
    -b track_com_5384_6520.root \   # reconstructed input
    -o temp.root \                  # converted output
    -f2                             # file type: 2=data, 1=MC
```

All options:

| Flag | Default | Meaning |
|---|---|---|
| `-b<file>` | `tracking_0_0_0.root` | Reconstructed ROOT file |
| `-o<file>` | `track_com.root` | Output file |
| `-f<1\|2>` | `1` | `1`=MC (with truth comparison), `2`=data (reco only) |
| `-a<file>` | `mcs-tracks.root` | MC truth file (only used with `-f1`) |
| `-t<tree>` | `T_rec_charge` | Reco tree name in `-b` file |
| `-n<tree>` | `T` | Truth tree name in `-a` file |

### Calibration loading

At startup, the tool reads `dQdx_scale` and `dQdx_offset` from the `Trun` tree
in the input file (if it exists):

```cpp
TTree *Trun = (TTree*)file1->Get("Trun");
if (Trun) {
    Trun->SetBranchAddress("dQdx_scale", &dQdx_scale);
    Trun->SetBranchAddress("dQdx_offset", &dQdx_offset);
    Trun->GetEntry(0);
}
```

These factors calibrate the raw charge values before they are processed.

---

## MC Mode (`-f1`): Truth Comparison

When `-f1` is specified, the tool also reads the truth MC track file
(`-a`). It applies a drift velocity correction to the truth x-coordinates:

```cpp
x_corrected = (x_raw + 0.6) / 1.098 * 1.1009999 - 0.1101
```

This converts from the MC simulation's drift velocity (1.098 mm/µs) to the
reconstructed drift velocity (1.1009999 mm/µs) used in reconstruction.

The corrected truth points are loaded into an `NFKDVec::Tree<double>` k-d tree
for fast nearest-neighbor queries.

For each reconstructed point, the tool finds the nearest truth point:

```cpp
auto [dist, closest_pt] = get_closest_point(pcloud, reco_pt);
```

This populates the comparison statistics per track.

---

## Output File: `T_rec` Tree

The output tree groups points by cluster ID (the `ndf` branch) into nested
vectors. Each entry in `T_rec` corresponds to one event (one fill per event):

| Branch | Type | Description |
|---|---|---|
| `rec_x`, `rec_y`, `rec_z` | `vector<vector<double>>` | Reco points grouped by cluster ID |
| `rec_dQ` | `vector<vector<double>>` | Charge per point per track |
| `rec_dx` | `vector<vector<double>>` | Hit count per point |
| `rec_L` | `vector<vector<double>>` | Cumulative arc length along track (cm) |
| `rec_N` | `vector<int>` | Number of points per track |
| `stat_total_L` | `vector<double>` | Total track length per cluster (cm) |
| `stat_max_dis` | `vector<double>` | Max distance from reco to truth point (cm) (MC only) |
| `stat_beg_dis` | `vector<double>` | Distance at the track start (cm) (MC only) |
| `stat_end_dis` | `vector<double>` | Distance at the track end (cm) (MC only) |
| `true_x`, `true_y`, `true_z` | `vector<double>` | Nearest truth point per reco point (MC only) |
| `true_dQ` | `vector<double>` | Nearest truth charge (MC only) |
| `com_dis` | `vector<vector<double>>` | Distance reco↔truth per point (MC only) |
| `com_dtheta` | `vector<vector<double>>` | Angular deviation reco↔truth per point (MC only) |

Additionally:
- `T_bad_ch`: cloned from input (dead channel list)
- `T_proj`, `T_proj_data`: cloned from input (2D projections)
- `T_true`: truth points after drift correction (MC only)

---

## Magnify-tracking Viewer

The output `temp.root` (or whatever `-o` specifies) is consumed by the
Magnify-tracking ROOT macro viewer at:

**https://github.com/BNLIF/Magnify-tracking**

The viewer reads `T_rec` and displays:
- 3D track trajectories (`rec_x/y/z` vectors per cluster).
- 2D projections from `T_proj` / `T_proj_data`.
- Dead channel overlays from `T_bad_ch`.
- In MC mode, truth vs. reco comparison with `com_dis` per-point coloring.

The typical workflow for validating reconstruction quality:

```bash
# Run reconstruction to produce tracking root file
wire-cell -c uboone-mabc.jsonnet -V tracking_output=tracking ...

# Convert to Magnify format (data mode)
wire-cell-uboone-magnify-tracking-convert -btracking_5384_6520.root -otemp.root -f2

# Open in ROOT
root temp.root
# Then load the Magnify-tracking macro from the GitHub repo
```

---

## Typical Validation Use

The `stat_max_dis` and `com_dis` quantities are the primary metrics for
reconstruction quality:
- `stat_max_dis < 1 cm`: the reconstructed track stays within 1 cm of truth
  along its entire length — indicates good reconstruction.
- Large `stat_beg_dis` or `stat_end_dis`: the track endpoints are
  misreconstructed (vertex or stopping point is wrong).
- `com_dtheta > 10°`: the local track direction deviates significantly from
  truth — may indicate a kink artifact or a poorly fitted segment.
