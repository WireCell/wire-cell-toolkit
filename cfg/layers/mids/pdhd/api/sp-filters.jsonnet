// PDHD OSP filter components (nominal).  The OSP C++ hard-codes these instance
// names.  Values from the dunereco pdhd sp-filters (nominal, non-APA1).

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
    sigma: 0.0,
  } + data,
};

[
  lf('ROI_loose_lf', { tau: 0.003 * wc.megahertz }),
  lf('ROI_tight_lf', { tau: 0.016 * wc.megahertz }),
  lf('ROI_tighter_lf', { tau: 0.08 * wc.megahertz }),

  hf('Gaus_tight'),
  hf('Gaus_wide', { sigma: 0.12 * wc.megahertz }),

  hf('Wiener_tight_U', { sigma: 0.221933 * wc.megahertz, power: 6.55413 }),
  hf('Wiener_tight_V', { sigma: 0.222723 * wc.megahertz, power: 8.75998 }),
  hf('Wiener_tight_W', { sigma: 0.225567 * wc.megahertz, power: 3.47846 }),

  hf('Wiener_wide_U', { sigma: 0.186765 * wc.megahertz, power: 5.05429 }),
  hf('Wiener_wide_V', { sigma: 0.1936 * wc.megahertz, power: 5.77422 }),
  hf('Wiener_wide_W', { sigma: 0.175722 * wc.megahertz, power: 4.37928 }),

  wf('Wire_ind', { sigma: 1.0 / wc.sqrtpi * 0.75 }),
  wf('Wire_col', { sigma: 1.0 / wc.sqrtpi * 10.0 }),
]
