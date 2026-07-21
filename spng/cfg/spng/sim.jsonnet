// Build chunks of sim.   See also drift.jsonnet and detsim.jsonnet

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(tpc, control) {
        
    anode: tpc.anode,

    depo_transform: pg.pnode({
        type: "DepoTransform",
        name: tpc.name,
        data: {
            anode: wc.tn(tpc.anode),
            dft: control.dft,
            rng: control.rng,

            drift_speed: tpc.lar.drift_speed,

            // FIXME: expose to user.
            first_frame_number: 0,
            flucutate: true,
            nsigma: 3,

            pirs: [wc.tn(pir) for pir in tpc.pirs],
            tick: tpc.ductor.tick,
            readout_time: tpc.ductor.readout_duration,
            start_time : tpc.ductor.start_time,

        },            
    }, nin=1, nout=1, uses=[tpc.anode]+tpc.pirs),

    // reframer
    reframer(extra_name=""):: pg.pnode({
        type: "Reframer",
        name: tpc.name + extra_name,
        data: {
            anode: wc.tn(tpc.anode),
            nticks: tpc.adc.readout_nticks,
        },
    }, nin=1, nout=1, uses=[tpc.anode]),
    

    empirical_model: {
        type: "EmpiricalNoiseModel",
        name: tpc.name,
        data: {
            anode: wc.tn(tpc.anode),
            dft: control.dft,
            rng: control.rng,
            nsamples: tpc.adc.readout_nticks,
            period: tpc.adc.tick,
            
            wire_length_scale: tpc.noise.empirical.wire_length_scale,
            spectra_file : tpc.noise.empirical.spectra_file,
            chanstat: tpc.noise.empirical.chanstat,
        },
        uses: [tpc.anode]
    },

    addnoise_empirical: pg.pnode({
        type: "AddNoise",
        name: tpc.name,
        data: {
            dft: control.dft,
            rng: control.rng,
            model: wc.tn($.empirical_model),
            nsamples: tpc.adc.readout_nticks,
            replacement_percentage: tpc.noise.empirical.replacement_percentage,
        }
    }, nin=1, nout=1, uses=[$.empirical_model]),

    // fixme: add other noise models, each getting an attribute in tpc.noise.{}.

    digitizer: pg.pnode({
        type: "Digitizer",
        name: tpc.name,
        data: {
            anode: wc.tn(tpc.anode),
            baselines: tpc.adc.baselines,
            fullscale: tpc.adc.fullscale,
            gain: tpc.adc.gain,
            resolution: tpc.adc.resolution,
        },
    }, nin=1, nout=1, uses=[tpc.anode]),

    splat: pg.pnode({
        type: "DepoFluxSplat",
        name: tpc.name,
        data: {
            anode: wc.tn(tpc.anode),
            field_response: wc.tn(tpc.fr),
            sparse: true,
            window_start: tpc.ductor.start_time,
            window_duration: tpc.ductor.readout_duration,
            tick: tpc.adc.tick,
            reference_time: 0.0,
            // Run wirecell-gen morse-* to find these numbers that match the extra
            // spread the sigproc induces.
        } + tpc.splat,
    }, nin=1, nout=1, uses=[tpc.anode, tpc.fr]),
}


