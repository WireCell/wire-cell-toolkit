// Drift-distance closure check of the OmnibusSigProc collection time offset
// (ctoffset) for FD-VD: simulation -> signal processing.
//
// Straight lines of charge along z at fixed y are placed at known distances
// from the collection (W) strip plane of one FD-VD CRU, simulated with
// DepoTransform and signal-processed with OmnibusSigProc.
// check-ctoffset.py converts each collection channel's peak time to a
// distance with the simulated drift speed and compares with the truth.
//
// The field response starts at response_plane = 181 mm from the W strips and
// assumes fr_speed = 1.53 mm/us.  When the simulated drift speed v differs,
// the collection signal arrives at a different time than the field response
// assumes, and ctoffset must be
//
//     ctoffset = response_plane / v - response_plane / fr_speed
//
// With this value the reconstructed distance is unbiased (up to a fixed
// ~-0.1 cm from tick rounding and peak finding) at every distance.
//
// The configuration is the FD-VD (dune10kt-vd) one from dunereco: v6 wires,
// PD-VD CRP field response, 14 mV/fC cold electronics (bottom) or top drift
// electronics (top), PD-VD noise spectra, its sigproc parameters.  All files
// are found via WIRECELL_PATH (wire-cell-data).
//
// TLAs (wire-cell -A name=value; all values may be strings):
//   drift   : "bottom" (anode 0) or "top" (anode 160)
//   ctoffset: "formula" (above), or a number in us (e.g. 4, the PD-VD data value)
//   noise   : "false" or "true", add electronics noise
//   speed   : simulated drift speed in mm/us.  Default 1.6056, the LArSoft
//             E-field parametrization at 0.5 kV/cm and 87 K (FD-VD nominal).
//   outname : FrameFileSink output, holds the "gauss<ident>" frame

local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';

// Distances of the lines from the W plane [cm].  Lines under 18.1 cm start
// between the anode and the response plane and are drifted backward.
// check-ctoffset.py must use the same list and z segments.
local distances = [0.5, 5, 15, 50, 150, 300, 500, 640];
local seg_z0 = 20 * wc.mm;
local seg_dz = 180 * wc.mm;
local seg_len = 160 * wc.mm;
local line_y = -4200 * wc.mm;  // face 0 of anodes 0 and 160: y in [-5046, -3369] mm

local tick = 0.5 * wc.us;
local nticks = 8500;
local response_plane = 181 * wc.mm;
local fr_speed = 1.53 * wc.mm / wc.us;

// Drift planes of the v6 geometry: anode at the W strips, response plane
// 181 mm into the drift, cathode at the TPC active-volume face.
local bottom_face = { anode: -6519.9 * wc.mm, response: (-6519.9 + 181) * wc.mm, cathode: -19.2 * wc.mm };
local top_face = { anode: 6521.5 * wc.mm, response: (6521.5 - 181) * wc.mm, cathode: 20.8 * wc.mm };

local tobool(x) = x == true || x == 'true';
local tonum(x) = if std.isNumber(x) then x else std.parseJson(x);

function(drift='bottom', ctoffset='formula', noise='false', speed=1.6056288906533874,
         outname='check-ctoffset.tar.bz2')

  local is_top = if drift == 'top' then true
                 else if drift == 'bottom' then false
                 else error 'drift must be "bottom" or "top", got "%s"' % drift;
  local ident = if is_top then 160 else 0;
  local face = if is_top then top_face else bottom_face;
  local drift_speed = tonum(speed) * wc.mm / wc.us;
  local sp_ctoffset = if ctoffset == 'formula'
                      then response_plane / drift_speed - response_plane / fr_speed
                      else tonum(ctoffset) * wc.us;

  // ---- services -----------------------------------------------------------
  local random = { type: 'Random', data: { generator: 'default', seeds: [0, 1, 2, 3, 4] } };
  local dft = { type: 'FftwDFT' };
  local wires = { type: 'WireSchemaFile', data: { filename: 'dunevd10kt_3view_30deg_v6_refactored.json.bz2' } };
  local anode = {
    type: 'AnodePlane',
    name: 'crm%d' % ident,
    data: { ident: ident, nimpacts: 10, wire_schema: wc.tn(wires), faces: [face, face] },
    uses: [wires],
  };
  local field = { type: 'FieldResponse', data: { filename: 'protodunevd_FR_imbalance3p_260501.json.bz2' } };

  // Response objects are padded beyond the readout as in the production config.
  local resp_nticks = nticks + 225;
  local gain = 14 * wc.mV / wc.fC;
  local elec = if is_top then {
    type: 'JsonElecResponse',
    data: { filename: 'dunevd-coldbox-elecresp-top-psnorm_400.json.bz2', gain: gain,
            shaping: 2 * wc.us, postgain: 1.36, tick: tick, nticks: resp_nticks },
  } else {
    type: 'ColdElecResponse',
    data: { filename: '', gain: gain, shaping: 2.2 * wc.us, postgain: 1.0,
            tick: tick, nticks: resp_nticks },
  };
  local pirs = [{
    type: 'PlaneImpactResponse',
    name: 'PIRplane%d' % plane,
    data: {
      plane: plane,
      dft: wc.tn(dft),
      field_response: wc.tn(field),
      short_responses: [wc.tn(elec)],
      overall_short_padding: 0.2 * wc.ms,
      long_responses: [],
      long_padding: 1.5 * wc.ms,
      tick: tick,
      nticks: resp_nticks,
    },
    uses: [dft, field, elec],
  } for plane in [0, 1, 2]];

  // ---- simulation ---------------------------------------------------------
  local wplane = face.anode;
  local sign = if is_top then -1 else 1;
  local depos = g.pnode({
    type: 'TrackDepos',
    data: {
      step_size: 1 * wc.mm,
      tracks: [
        local x = wplane + sign * distances[i] * wc.cm;
        local z0 = seg_z0 + seg_dz * i;
        {
          time: 0,
          charge: -5000,  // electrons per step
          ray: wc.ray(wc.point(x, line_y, z0, 1), wc.point(x, line_y, z0 + seg_len, 1)),
        }
        for i in std.range(0, std.length(distances) - 1)
      ],
    },
  }, nin=0, nout=1);

  local drifter = g.pnode({
    type: 'Drifter',
    data: {
      rng: wc.tn(random),
      xregions: [bottom_face, top_face],
      time_offset: 0,
      drift_speed: drift_speed,
      fluctuate: true,
      DL: 6.2e-9 * wc.cm2 / wc.ns,
      DT: 1.63e-8 * wc.cm2 / wc.ns,
      lifetime: 10.4 * wc.ms,
    },
  }, nin=1, nout=1, uses=[random]);

  // The ductor opens response_plane/v early; the Reframer cuts that off again.
  local response_time_offset = response_plane / drift_speed;
  local response_nticks = wc.roundToInt(response_time_offset / tick);
  local ductor_start = -response_time_offset;
  local ductor_readout = (nticks + response_nticks) * tick;

  local bagger = g.pnode({
    type: 'DepoBagger',
    data: { gate: [ductor_start, ductor_start + ductor_readout] },
  }, nin=1, nout=1);

  local ductor = g.pnode({
    type: 'DepoTransform',
    name: 'ductor%d' % ident,
    data: {
      rng: wc.tn(random),
      dft: wc.tn(dft),
      anode: wc.tn(anode),
      pirs: [wc.tn(p) for p in pirs],
      fluctuate: true,
      drift_speed: drift_speed,
      first_frame_number: 100,
      readout_time: ductor_readout,
      start_time: ductor_start,
      tick: tick,
      nsigma: 3,
    },
  }, nin=1, nout=1, uses=[random, dft, anode] + pirs);

  local reframer = g.pnode({
    type: 'Reframer',
    name: 'reframer%d' % ident,
    data: { anode: wc.tn(anode), tags: [], fill: 0.0, tbin: response_nticks, toffset: 0, nticks: nticks },
  }, nin=1, nout=1, uses=[anode]);

  local noise_model = {
    type: 'EmpiricalNoiseModel',
    name: 'noise%d' % ident,
    data: {
      anode: wc.tn(anode),
      dft: wc.tn(dft),
      chanstat: '',
      spectra_file: if is_top then 'pdvd-top-noise-spectra-v3.json.bz2'
                    else 'pdvd-bottom-noise-spectra-14mVfC-v1.json.bz2',
      nsamples: nticks,
      period: tick,
      wire_length_scale: 1.0 * wc.cm,
    },
    uses: [anode, dft],
  };
  local addnoise = g.pnode({
    type: 'AddNoise',
    name: 'addnoise%d' % ident,
    data: { rng: wc.tn(random), dft: wc.tn(dft), model: wc.tn(noise_model),
            nsamples: nticks, replacement_percentage: 0.02 },
  }, nin=1, nout=1, uses=[random, dft, noise_model]);

  local adc = if is_top then {
    baselines: [1.0 * wc.volt, 1.0 * wc.volt, 1.0 * wc.volt],
    fullscale: [0.0 * wc.volt, 2.0 * wc.volt],
  } else {
    baselines: [1003.4 * wc.millivolt, 1003.4 * wc.millivolt, 507.7 * wc.millivolt],
    fullscale: [0.2 * wc.volt, 1.6 * wc.volt],
  };
  local resolution = 14;
  local digitizer = g.pnode({
    type: 'Digitizer',
    name: 'digitizer%d' % ident,
    data: adc { anode: wc.tn(anode), gain: 1.0, resolution: resolution, frame_tag: 'orig%d' % ident },
  }, nin=1, nout=1, uses=[anode]);

  // ---- signal processing --------------------------------------------------
  // Filters as dune10kt-vd sp-filters.jsonnet; OmnibusSigProc finds them by name.
  local lf(name, tau) = { type: 'LfFilter', name: name, data: { max_freq: 1 * wc.megahertz, tau: tau } };
  local hf(name, sigma, power, flag=true, max_freq=1 * wc.megahertz) = {
    type: 'HfFilter', name: name, data: { max_freq: max_freq, sigma: sigma, power: power, flag: flag },
  };
  local filters = [
    lf('ROI_tight_lf', 0.014 * wc.megahertz),
    lf('ROI_tighter_lf', 0.06 * wc.megahertz),
    lf('ROI_loose_lf', 0.002 * wc.megahertz),
    hf('Gaus_tight', 0.0, 2),
    hf('Gaus_wide', 0.12 * wc.megahertz, 2),
    hf('Wiener_tight_U', 0.148788 * wc.megahertz, 3.76194),
    hf('Wiener_tight_V', 0.1596568 * wc.megahertz, 4.36125),
    hf('Wiener_tight_W', 0.13623 * wc.megahertz, 3.35324),
    hf('Wiener_wide_U', 0.186765 * wc.megahertz, 5.05429),
    hf('Wiener_wide_V', 0.1936 * wc.megahertz, 5.77422),
    hf('Wiener_wide_W', 0.175722 * wc.megahertz, 4.37928),
    hf('Wire_ind', 0.75 / std.sqrt(std.pi), 2, false, 1.0),
    hf('Wire_col', 3.0 / std.sqrt(std.pi), 2, false, 1.0),
  ];

  local fullscale = adc.fullscale[1] - adc.fullscale[0];
  local sigproc = g.pnode({
    type: 'OmnibusSigProc',
    name: 'sigproc%d' % ident,
    data: {
      anode: wc.tn(anode),
      dft: wc.tn(dft),
      field_response: wc.tn(field),
      elecresponse: wc.tn(elec),
      ftoffset: 0.0,
      ctoffset: sp_ctoffset,
      per_chan_resp: '',
      fft_flag: 0,
      postgain: 1.0,
      ADC_mV: ((1 << resolution) - 1) / fullscale,
      troi_col_th_factor: 5.0,
      troi_ind_th_factor: 3.0,
      lroi_rebin: 6,
      lroi_th_factor: 3.5,
      lroi_th_factor1: 0.7,
      lroi_jump_one_bin: 1,
      r_th_factor: 3.0,
      r_fake_signal_low_th: 375,
      r_fake_signal_high_th: 750,
      r_fake_signal_low_th_ind_factor: 1.0,
      r_fake_signal_high_th_ind_factor: 1.0,
      r_th_peak: 3.0,
      r_sep_peak: 6.0,
      r_low_peak_sep_threshold_pre: 1200,
      wiener_tag: 'wiener%d' % ident,
      gauss_tag: 'gauss%d' % ident,
      use_roi_debug_mode: false,
      use_multi_plane_protection: false,
    },
  }, nin=1, nout=1, uses=[anode, dft, field, elec] + filters);

  local sink = g.pnode({
    type: 'FrameFileSink',
    name: outname,
    data: { outname: outname, tags: ['gauss%d' % ident], digitize: false },
  }, nin=1, nout=0);

  local graph = g.pipeline([depos, drifter, bagger, ductor, reframer]
                           + (if tobool(noise) then [addnoise] else [])
                           + [digitizer, sigproc, sink]);

  local app = { type: 'Pgrapher', data: { edges: g.edges(graph) } };
  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellSio', 'WireCellSigProc', 'WireCellAux'],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + g.uses(graph) + [app]
