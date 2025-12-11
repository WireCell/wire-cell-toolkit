// This provides various "full-sub-graph" constructions.  Each function is
// provided with an upstream source node and a downstream sink node of the
// proper multiplicity.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local detsim = import "detsim.jsonnet";
local frame_mod = import "frame.jsonnet";

function(control={})
{
    local frame = frame_mod(control),

    /// Connect depo source to adc frame sink with drift and detsim.  Sink is
    /// expected to have one input port per tpc.
    depos_to_adc_frame(det, depos_source, adc_sink)::
        local sim = detsim(det, control);
        local nout = std.length(sim.oports);
        local head = pg.pipeline([depos_source, sim]);
        pg.shuntline(head, adc_sink),

    /// Connect ADC frame source to ADC tensorset sink as frame to TDM tensor
    /// converters in between.  Source and sink expected to have per-TPC ports.
    frame_to_tensorset(det, adc_frame_source, adc_tensorset_sink)::
        local converters = pg.crossline([
            frame.to_tdm(tpc)
            for tpc in det.tpcs]);
        pg.shuntlines([adc_frame_source, converters, adc_tensorset_sink]),



}
