// A wire-cell job to measure how trigger activities (TAs) select supernova
// neutrino interactions over radiological backgrounds in one DUNE FD HD APA
// or VD CRU (see detectors.jsonnet).
//
// Per frame (an independent readout window of nticks):
//
//   LowEDepoSource ─▶ DepoSetFanout ─┬─▶ DepoSetDrifter ─▶ DepoTransform ─▶ Reframer
//    (MARLEY SNB +                   │   ─▶ AddNoise ─▶ Digitizer ─▶ TPFinder ─▶ TAFinder ─┐
//     radiologicals)                 │                                                     ▼
//                                    └───────────── depos (truth) ──────────────────▶ TATruth ─▶ TensorFileSink
//
// The SNB rate is set far above any physical one (snb_rate_hz) to spend the
// compute on signal; the metadata records it with a physical reference rate
// (snb_physical_rate_hz) so the analysis can rescale (trig/pytools).
//
// Top-level arguments (strings, as from "wire-cell -A name=value"):
//
// - detector :: "hd" (default) or "vd".
// - anode :: which APA (hd, default 0) or CRU (vd, default 160, top drift).
// - marley :: MARLEY JSON events file, or "" for no SNB (default).
// - outfile :: output .npz of truth tables.
// - nframes, seed :: number of frames, random seed.
// - nticks :: readout window in ticks of 500 ns (default 6000).
// - lifetime_ms :: electron lifetime, default that of the detector params.
// - snb_rate_hz :: simulated SNB interaction rate in the anode's volume.
// - snb_physical_rate_hz_per_t :: the reference physical rate per tonne of
//   argon, default 0.0075 Hz/t: about 3000 events in 40 kt over 10 s, a
//   10 kpc burst (model dependent).  The analysis takes the boost as
//   snb_rate_hz over this times the simulated mass.
// - radiologicals :: comma separated names from radiologicals.jsonnet,
//   where "all" and "stress" name lists, or "none".
// - source_boosts :: comma separated name:factor multiplying the rate of
//   background sources (eg "K42:50,Bi214:10"), to measure rare tails.  A
//   chained daughter inherits its parent's boost.  Recorded in the
//   metadata as source_boosts; the analysis divides rates by them.
// - noise :: "true" or "false".
// - thresholds :: TPG thresholds in ADC for U, V, W, comma separated; a
//   negative one disables the plane (default per detector, as dunetrigger).
// - ta_algorithm :: "Supernova" or "SWIFT".
// - plane_offsets_us :: response plane to TP start delay per plane, for
//   truth matching (default per detector, calibrate with "offset_est").
// - time_tolerance_us :: truth matching tolerance (default 10).
// - depo_out, frame_out :: if given, also save the depos (before drift) and
//   the digitized frames to these files (.tar.bz2), for displays.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local radio = import "radiologicals.jsonnet";
local detectors = import "detectors.jsonnet";

local split(s) = std.map(function(x) std.stripChars(x, " "), std.split(s, ","));
local numbers(s) = std.map(std.parseJson, split(s));

function(detector="hd", anode="", marley="", outfile="snb-ta-study.npz", nframes="10", seed="1",
         nticks="6000", lifetime_ms="", snb_rate_hz="1000", snb_physical_rate_hz_per_t="0.0075",
         radiologicals="all", noise="true", thresholds="", ta_algorithm="Supernova",
         plane_offsets_us="", time_tolerance_us="10", depo_out="", frame_out="",
         source_boosts="")

local det = detectors[detector](
    anode=if anode == "" then (if detector == "vd" then 160 else 0) else std.parseJson(anode),
    nticks=std.parseJson(nticks),
    lifetime=if lifetime_ms == "" then null else std.parseJson(lifetime_ms)*wc.ms,
    noise=noise == "true");
local params = det.params;
local tick = params.daq.tick;
local readout_time = params.daq.nticks * tick;
local thr = if thresholds == "" then det.thresholds else numbers(thresholds);
local offsets = if plane_offsets_us == "" then det.plane_offsets_us else numbers(plane_offsets_us);
local snb_rate = std.parseJson(snb_rate_hz);
local snb_ref_per_t = std.parseJson(snb_physical_rate_hz_per_t);

local given_boosts = if source_boosts == "" then {} else {
    [std.split(kv, ":")[0]]: std.parseJson(std.split(kv, ":")[1]) for kv in split(source_boosts)
};
local radio_unboosted =
    if radiologicals == "none" then []
    else std.flattenArrays([
        if name == "all" || name == "stress" then radio[name] else [radio[name]]
        for name in split(radiologicals)]);
// daughters take their parent's boost
local boosts = given_boosts + {
    [src.daughter.name]: given_boosts[src.name]
    for src in radio_unboosted
    if std.objectHas(src, "daughter") && std.objectHas(given_boosts, src.name)
};
local boosted(src) =
    if !std.objectHas(boosts, src.name) then src
    else if std.objectHas(src, "activity_bq_per_kg")
    then src { activity_bq_per_kg: super.activity_bq_per_kg * boosts[src.name] }
    else src { rate_hz: super.rate_hz * boosts[src.name] };
local radio_sources = [boosted(src) for src in radio_unboosted];
local sources = (if marley == "" then []
                 else [{name: "snb", kind: "marley", file: marley, rate_hz: snb_rate, signal: true}])
                + radio_sources;

local recomb = {type: "BoxRecombination", name: "", data: {}};

local source = g.pnode({
    type: "LowEDepoSource",
    data: {
        anode: wc.tn(det.anode),
        recombination: wc.tn(recomb),
        seed: std.parseJson(seed),
        nframes: std.parseJson(nframes),
        start_time: 0,
        readout_time: readout_time,
        // the longest drift plus a margin
        pre_time: (det.drift_length + 10*wc.cm) / params.lar.drift_speed,
        sources: sources,
    },
}, nin=0, nout=1, uses=[det.anode, recomb]);

local ndepo_out = if depo_out == "" then 2 else 3;
local fanout = g.pnode({
    type: "DepoSetFanout",
    data: {multiplicity: ndepo_out},
}, nin=1, nout=ndepo_out);

local setdrifter = g.pnode({
    type: "DepoSetDrifter",
    data: {drifter: wc.tn(det.drifter)},
}, nin=1, nout=1, uses=[det.drifter]);

local tpfinder = g.pnode({
    type: "TPFinder",
    data: {
        period: tick,
        anode: wc.tn(det.anode),
        streaming: false,
        tpg: {pedestal_init: "mode", filter: "none", threshold: thr[2]},
        planes: [{threshold: t} for t in thr],
    },
}, nin=1, nout=1, uses=[det.anode]);

local ta_configs = {
    Supernova: {both: true},
    SWIFT: {},
};
local tafinder = g.pnode({
    type: "TAFinder",
    data: {
        algorithm: ta_algorithm,
        config: ta_configs[ta_algorithm],
        anode: wc.tn(det.anode),
        period: tick,
        streaming: false,
        store_inputs: true,
    },
}, nin=1, nout=1, uses=[det.anode]);

local truth = g.pnode({
    type: "TATruth",
    data: {
        anode: wc.tn(det.anode),
        drift_speed: params.lar.drift_speed,
        lifetime: params.lar.lifetime,
        period: tick,
        plane_offsets: [x*wc.us for x in offsets],
        time_tolerance: std.parseJson(time_tolerance_us)*wc.us,
        channel_tolerance: 1,
        sources: sources,
        metadata: {
            detector: detector,
            geometry: det.geometry,
            snb_rate_hz: snb_rate,
            snb_physical_rate_hz_per_t: snb_ref_per_t,
            source_boosts: boosts,
            readout_time: readout_time,
            tick: tick,
            thresholds: thr,
            ta_algorithm: ta_algorithm,
            noise: noise == "true",
            nframes: std.parseJson(nframes),
            lifetime: params.lar.lifetime,
            drift_speed: params.lar.drift_speed,
            adc_resolution: params.adc.resolution,
        },
    },
}, nin=2, nout=1, uses=[det.anode]);

local sink = g.pnode({
    type: "TensorFileSink",
    data: {outname: outfile, prefix: ""},
}, nin=1, nout=0);

// Optional taps for displays.
local depo_sink = g.pnode({
    type: "DepoFileSink",
    data: {outname: depo_out},
}, nin=1, nout=0);
local frame_fanout = g.pnode({
    type: "FrameFanout",
    name: "frametap",
    data: {multiplicity: 2},
}, nin=1, nout=2);
local frame_sink = g.pnode({
    type: "FrameFileSink",
    data: {outname: frame_out, tags: [], digitize: true, masks: false},
}, nin=1, nout=0);

local simchain =
    if frame_out == ""
    then g.pipeline([setdrifter, det.simpipe, tpfinder, tafinder])
    else g.intern(innodes=[setdrifter], outnodes=[tafinder],
                  centernodes=[det.simpipe, frame_fanout, tpfinder, frame_sink],
                  edges=[g.edge(setdrifter, det.simpipe),
                         g.edge(det.simpipe, frame_fanout),
                         g.edge(frame_fanout, tpfinder, 0, 0),
                         g.edge(frame_fanout, frame_sink, 1, 0),
                         g.edge(tpfinder, tafinder)]);

local graph = g.intern(
    innodes=[source],
    outnodes=[sink],
    centernodes=[fanout, simchain, truth] + (if depo_out == "" then [] else [depo_sink]),
    edges=[
        g.edge(source, fanout),
        g.edge(fanout, simchain, 0, 0),
        g.edge(fanout, truth, 1, 0),
        g.edge(simchain, truth, 0, 1),
        g.edge(truth, sink),
    ] + (if depo_out == "" then [] else [g.edge(fanout, depo_sink, 2, 0)]));

local app = {
    type: "Pgrapher",
    data: {edges: g.edges(graph)},
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellAux", "WireCellGen", "WireCellSigProc", "WireCellPgraph", "WireCellSio",
                  "WireCellTrig"],
        apps: ["Pgrapher"],
    },
};

[cmdline] + g.uses(graph) + [app]
