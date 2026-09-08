// PDVD OSP filter components (nominal, bottom drift).  The OSP C++ hard-codes
// these instance names.  PDVD's production config splits every filter into
// "_b"/"_t" instances, but the single-APA layers/spdir job uses one bottom
// anode, so here we register the bottom-drift values under the bare names the
// OSP defaults expect.  Values from the dunereco protodunevd sp-filters.

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
local wf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1,  // warning: units
    power: 2,
    flag: false,
    sigma: 0.0,
  } + data,
};

[
  lf('ROI_tight_lf', { tau: 0.014 * wc.megahertz }),
  lf('ROI_tighter_lf', { tau: 0.06 * wc.megahertz }),
  lf('ROI_loose_lf', { tau: 0.003 * wc.megahertz }),

  hf('Gaus_tight'),
  hf('Gaus_wide', { sigma: 0.12 * wc.megahertz }),

  hf('Wiener_tight_U', { sigma: 0.148788 * wc.megahertz, power: 3.76194 }),
  hf('Wiener_tight_V', { sigma: 0.1596568 * wc.megahertz, power: 4.36125 }),
  hf('Wiener_tight_W', { sigma: 0.13623 * wc.megahertz, power: 3.35324 }),

  hf('Wiener_wide_U', { sigma: 0.186765 * wc.megahertz, power: 5.05429 }),
  hf('Wiener_wide_V', { sigma: 0.1936 * wc.megahertz, power: 5.77422 }),
  hf('Wiener_wide_W', { sigma: 0.175722 * wc.megahertz, power: 4.37928 }),

  wf('Wire_ind', { sigma: 1.0 / wc.sqrtpi * 5.0 }),
  wf('Wire_col', { sigma: 1.0 / wc.sqrtpi * 10.0 }),
]
