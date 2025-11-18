// This provides various "full-graph" constructions.
//
// Each function takes a "system" object, see system.jsonnet and par-function arguments and options.
//
// A top level Jsonnet may then define a system and call job function.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local detsim = import "detsim.jsonnet";

{
    depos_to_adc(system, input_filename, output_filename_pattern)::
        local sim = detsim(system.scope.det, system.scope.tpcs, system.control);
        local nout = std.length(sim.oports);
        local output_filenames = [output_filename_pattern % ind for ind in wc.iota(nout)];
        local depos_in = system.io.depo_source(input_filename);
        local adc_outs = pg.crossline([system.io.frame_sink(ofn, digitize=true) for ofn in output_filenames]);
        local head = pg.pipeline([depos_in, sim]);
        pg.shuntline(head, adc_outs)
        // local tail = pg.shuntline(sim, adc_outs);
        // pg.components([head, tail]),
}
