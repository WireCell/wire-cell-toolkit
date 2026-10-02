// Counterpart of fdvd_sim/cfg/params-fdvd.jsonnet (wcp WireCell/wcp-porting-validation), promoted
// unchanged for the in-toolkit FD-VD low-energy reconstruction (fdvd_sim doc 16).
//
// fdvd_sim doc 02: DUNE FD-VD 1x8x14 (3view 30deg, v7 wires, top drift only) simulation
// parameters, SYNCED TO THE PDVD SIMULATION.
//
// Base = pgrapher/experiment/protodunevd/simparams.jsonnet (PDVD sim, which layers on
// protodunevd/params.jsonnet), so every PDVD sim lesson is inherited unless overridden here:
// 14-bit ADC, top electronics JsonElecResponse dunevd-coldbox-elecresp-top-psnorm_400 x 1.36,
// field response protodunevd_FR_imbalance3p_260501 (18.1 cm response plane), top noise
// pdvd-top-noise-spectra-v3, no RC layers (rc_layers 0), overall_short_padding 0.2 ms,
// DL/DT = 4.0/8.8 cm2/s.
//
// Overridden for FD-VD (owner ruling 2026-09-28, doc 02 sec 1):
//   * geometry: 112 CRMs, all TOP drift (v7 wires: every collection plane at x = +3250.7 mm),
//     cathode at x = -325 cm;
//   * transport: drift speed 1.60563 mm/us and lifetime 10.4 ms = the LArSoft FD-VD sample's
//     values (so light t0 <-> x is consistent); DL/DT kept at the PDVD sim pair 4.0/8.8;
//   * every anode uses the PDVD TOP electronics chain (elecs[0] = elecs[1] = top, both noise
//     entries = top v3, both FR entries the same file) -- the PDVD sim.jsonnet/sp.jsonnet
//     select bottom settings for ident < 4, which is why sim-fdvd/sp-fdvd are forks;
//   * readout: 8500 ticks x 0.5 us = [0, 4.25] ms, tick 0 at t = 0 (the official FD-VD
//     TPC window, DefaultTrigTime 0; PDVD's -250 us pre-trigger is a PDVD-data convention).

local wc = import 'wirecell.jsonnet';
local base = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';

base {
    det: {
        // v7 wires: U/V/W at x = 3250.3/3250.5/3250.7 mm (LArSoft 0.2 mm-step convention).
        local w_x = 3250.7 * wc.mm,
        // PDVD lesson (protodunevd/params.jsonnet:41-52): the drift-facing boundary of the
        // active LAr is the CRP shield plane, 16.4 mm below W (3.2 + 10 + 3.2 mm stack).
        // FD-VD uses the same CRPs.
        local apa_plane = 16.4 * wc.mm,
        // FR origin: 18.1 cm from the collection plane (protodunevd_FR_imbalance3p_260501).
        local resp = 18.1 * wc.cm,
        response_plane: resp,
        local cathode_x = -325.0 * wc.cm,

        local face = {
            anode: w_x - apa_plane,
            response: w_x - resp,
            cathode: cathode_x,
        },
        volumes: [
            {
                wires: n,
                name: 'crm%d' % n,
                faces: [face, face],   // one-face CRM; same convention as dune-vd/params.jsonnet
            } for n in std.range(0, 111)
        ],
        bounds: {
            tail: wc.point(-3.25, -6.751, 0.0085, wc.m),
            head: wc.point(3.251, 6.751, 21.0033, wc.m),
        },
    },

    lar: super.lar {
        drift_speed: 1.60563 * wc.mm / wc.us,   // FD-VD LArSoft (wirecell_dune.fcl driftSpeed)
        lifetime: 10.4 * wc.ms,                 // FD-VD LArSoft detsim value
        DL: 4.0 * wc.cm2 / wc.s,                // PDVD sim pair (protodunevd/simparams.jsonnet)
        DT: 8.8 * wc.cm2 / wc.s,
    },

    daq: super.daq {
        nticks: 8500,
    },

    adc: super.adc {
        resolution: 14,
        // PDVD top (TDE) digitizer: ~1 V baselines, 0-2 V full scale (protodunevd/sim.jsonnet:32-37)
        baselines: [1.0 * wc.volt, 1.0 * wc.volt, 1.0 * wc.volt],
        fullscale: [0.0 * wc.volt, 2.0 * wc.volt],
    },

    local top_elec = super.elecs[1],
    elecs: [top_elec, top_elec],
    elec: top_elec,

    sim: super.sim {
        local tick0_time = 0.0 * wc.us,
        local response_time_offset = $.det.response_plane / $.lar.drift_speed,
        local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),
        ductor: {
            nticks: $.daq.nticks + response_nticks,
            readout_time: self.nticks * $.daq.tick,
            start_time: tick0_time - response_time_offset,
        },
        reframer: {
            tbin: response_nticks,
            nticks: $.daq.nticks,
        },
    },

    files: {
        wires: 'dunevd10kt_3view_30deg_v7_refactored_1x8x14.json.bz2',
        fields: [
            'protodunevd_FR_imbalance3p_260501.json.bz2',
            'protodunevd_FR_imbalance3p_260501.json.bz2',
        ],
        noises: [
            'pdvd-top-noise-spectra-v3.json.bz2',
            'pdvd-top-noise-spectra-v3.json.bz2',
        ],
        chresp: null,
    },
}
