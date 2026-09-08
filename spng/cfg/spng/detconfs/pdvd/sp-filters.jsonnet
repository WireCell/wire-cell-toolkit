// SPNG "mirror OSP" filter components for ProtoDUNE-VD.
//
// This is a faithful copy of the official
//   reference/dunereco/dunereco/DUNEWireCell/protodunevd/sp-filters.jsonnet
// (byte-identical to toolkit/cfg/pgrapher/experiment/protodunevd/sp-filters.jsonnet).
//
// WARNING: the OSP C++ hard-codes these filter instance names.  PDVD splits
// every filter into "_b" (bottom drift, anode ident 0..3) and "_t" (top drift,
// ident 4..7) instances; sp.jsonnet selects the suffix via (anode.data.ident<4).
// The values are currently identical for _b and _t but are kept separate so
// they can diverge without code changes.

local wc = import 'wirecell.jsonnet';

local lf(name, data={}) = {
  type: 'LfFilter',
  name: name,
  data: {
    max_freq: 1 * wc.megahertz,
    tau: 0.0 * wc.megahertz,
  } + data,
};
local hf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1 * wc.megahertz,
    sigma: 0.0 * wc.megahertz,
    power: 2,
    flag: true,
  } + data,
};
// All "wire" filters are Hf with different base values.
local wf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1,  // warning: units
    power: 2,
    flag: false,
    sigma: 0.0,  // caller should provide
  } + data,
};

[
  lf('ROI_tight_lf_b',   { tau: 0.014 * wc.megahertz }),
  lf('ROI_tight_lf_t',   { tau: 0.014 * wc.megahertz }),
  lf('ROI_tighter_lf_b', { tau: 0.06  * wc.megahertz }),
  lf('ROI_tighter_lf_t', { tau: 0.06  * wc.megahertz }),
  lf('ROI_loose_lf_b',   { tau: 0.003 * wc.megahertz }),
  lf('ROI_loose_lf_t',   { tau: 0.003 * wc.megahertz }),

  hf('Gaus_tight'),
  hf('Gaus_wide_b', { sigma: 0.12 * wc.megahertz }),
  hf('Gaus_wide_t', { sigma: 0.12 * wc.megahertz }),

  hf('Wiener_tight_U_b', { sigma: 0.148788  * wc.megahertz, power: 3.76194 }),
  hf('Wiener_tight_U_t', { sigma: 0.148788  * wc.megahertz, power: 3.76194 }),
  hf('Wiener_tight_V_b', { sigma: 0.1596568 * wc.megahertz, power: 4.36125 }),
  hf('Wiener_tight_V_t', { sigma: 0.1596568 * wc.megahertz, power: 4.36125 }),
  hf('Wiener_tight_W_b', { sigma: 0.13623   * wc.megahertz, power: 3.35324 }),
  hf('Wiener_tight_W_t', { sigma: 0.13623   * wc.megahertz, power: 3.35324 }),

  hf('Wiener_wide_U_b',  { sigma: 0.186765  * wc.megahertz, power: 5.05429 }),
  hf('Wiener_wide_U_t',  { sigma: 0.186765  * wc.megahertz, power: 5.05429 }),
  hf('Wiener_wide_V_b',  { sigma: 0.1936    * wc.megahertz, power: 5.77422 }),
  hf('Wiener_wide_V_t',  { sigma: 0.1936    * wc.megahertz, power: 5.77422 }),
  hf('Wiener_wide_W_b',  { sigma: 0.175722  * wc.megahertz, power: 4.37928 }),
  hf('Wiener_wide_W_t',  { sigma: 0.175722  * wc.megahertz, power: 4.37928 }),

  wf('Wire_ind_b', { sigma: 1.0 / wc.sqrtpi * 5.0 }),
  wf('Wire_ind_t', { sigma: 1.0 / wc.sqrtpi * 5.0 }),
  wf('Wire_col_b', { sigma: 1.0 / wc.sqrtpi * 10.0 }),
  wf('Wire_col_t', { sigma: 1.0 / wc.sqrtpi * 10.0 }),
]
