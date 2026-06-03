# SBND PMT VUV/VIS efficiencies and the "group of five"

How the per-OpDet `VUVEfficiency` / `VISEfficiency` values in `QLMatching` map onto
the physical SBND PMT layout — specifically the side-vs-central PMTs of each
quintuplet.

## The quintuplet geometry

Behind each SBND anode wall there are **120 PMTs**, tiled as **24 quintuplets**.
Each quintuplet is *four side PMTs around one central PMT*:

- **1 central PMT** — TPB-**uncoated**, on an "offset" ring (e.g. ch36 at y=−135, z=53.67 cm).
- **4 side PMTs** — TPB-**coated**, its four nearest neighbours, each **50 cm** away on
  the two adjacent "main" rings (e.g. ch6/ch8/ch60/ch62 around ch36).

The count follows directly: 24 × (1 + 4) = **24 uncoated + 96 coated = 120 PMTs**.
The 96:24 = 4:1 coated:uncoated ratio is exactly why every central uncoated PMT has
four equidistant coated neighbours.

Verified geometrically from the dumped OpDet table (`sbnd_xin/work/.../calib-evt2.json`,
fields `{ch,x,y,z,type,apa}`) and the PMT positions in
`wire-cell-bee3/events/static/js/bee/physics/experiment.js` (`SBND.updateOPLocation`):
for each uncoated PMT the four nearest PMTs are all coated, all at d = 50 cm.

## The efficiencies

Per-OpDet defaults baked into `match/inc/WireCellMatch/QLMatching.h`
(`m_VUVEfficiency` / `m_VISEfficiency`, indexed by OpDet over all 312 channels).
They are overridable via the `"VUVEfficiency"` / `"VISEfficiency"` config arrays, but
**no SBND jsonnet currently overrides them** — these defaults are what runs.

| PMT | role in the group | TPB | `VUVEfficiency` | `VISEfficiency` | count |
|-----|-------------------|-----|-----------------|-----------------|-------|
| **4 side** PMTs   | coated   | yes | **0.0392** | **0.0260** | 96 |
| **1 central** PMT | uncoated | no  | **0.0357** | **0.0357** | 24 |

### The other 192 OpDets — X-ARAPUCAs

X-ARAPUCAs are `type` 0 (PMTs are `type` 1). They come in two variants — **VUV** and
**VIS** — distinguished by the dichroic filter window, **not** by a TPB coating (the
coated/uncoated split is a PMT concept and does not apply here):

| X-ARAPUCA | filter sees | `VUVEfficiency` | `VISEfficiency` | count |
|-----------|-------------|-----------------|-----------------|-------|
| **VUV-XA** | direct 128 nm | **0.01752** | 0.00271 | 96 |
| **VIS-XA** | reflected ~430 nm | **0.00000** | **0.01264** | 96 |

Reading the two efficiency arrays by hand:

- **`VUVEfficiency`** — the value **0.00000** is the **VIS-XA** (its dichroic filter
  rejects 128 nm, so it is blind to direct VUV); **0.01752** is the **VUV-XA**
  (sensitive to the direct scintillation).
- **`VISEfficiency`** — the values **0.01264** and **0.00271** are *both* X-ARAPUCA
  (variant, not coating): **0.01264** = VIS-XA (its main channel, reflected visible),
  **0.00271** = VUV-XA (small residual visible response). The PMT VIS values are the
  *other* two entries in the array: 0.0260 (coated) and 0.0357 (uncoated).

So each array carries four distinct levels — two PMT (coated/uncoated) and two
X-ARAPUCA (VUV/VIS):

| level meaning | `VUVEfficiency` | `VISEfficiency` |
|---------------|-----------------|-----------------|
| coated PMT (side)     | 0.03920 | 0.02600 |
| uncoated PMT (central)| 0.03570 | 0.03570 |
| VUV X-ARAPUCA         | 0.01752 | 0.00271 |
| VIS X-ARAPUCA         | 0.00000 | 0.01264 |

Geometrically the X-ARAPUCAs are not part of the PMT quintuplet: VUV-XA and VIS-XA sit
on **separate z-rings** (e.g. TPC0 VIS-XA at z≈16.05 cm, VUV-XA at z≈31.29 cm), three
per ring at y = −135 / 0 / +135 cm — interleaved bands rather than a 4-around-1 cluster.

## Why side and central differ

The semi-analytical model predicts per-channel light as

```
pred_PE_ch = VUVeff_ch · (direct visibility) + VISeff_ch · (reflected visibility)
```

where **VUV** = the direct 128 nm scintillation, and **VIS** = the reflected
component (light that bounced off the TPB-coated cathode/foils and arrives already
wavelength-shifted to ~430 nm).

- **Side (coated) PMTs** carry TPB on the window, so they convert and detect the
  **direct 128 nm** efficiently (VUV 0.0392) plus the reflected visible (VIS 0.0260)
  — VUV-dominant.
- **Central (uncoated) PMT** has no shifter, so it cannot enhance the direct VUV the
  way a coated PMT does; its two efficiencies are set **equal (0.0357)** — the model
  gives it the same response to both visibility components rather than a TPB-boosted
  VUV term.

So summing predicted light over a quintuplet, the four coated PMTs are the
VUV-sensitive "eyes" and the central uncoated one responds the same to direct and
reflected light.

## Where this lives in code

- Arrays: `match/inc/WireCellMatch/QLMatching.h` (`m_VUVEfficiency`, `m_VISEfficiency`).
- Applied in: `match/src/SemiAnalyticalModel.*` via `QLMatching` (see
  [`semi-analytical-model.md`](semi-analytical-model.md)).
- OpDet `type` (1 = PMT, 0 = (X)Arapuca) loaded from the semimodel JSON
  (`m_semimodel_file`, default `sbnd/photodet/semi-analytical-sbnd.json`); the
  coated/uncoated split is read off the efficiency values, not a separate flag.
