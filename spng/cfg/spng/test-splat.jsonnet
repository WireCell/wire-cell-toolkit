// This produces a final graph to input depos and output ADC waveforms.

local det_mod = import "spng/det.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";


local detconf = import "spng/detconf.jsonnet";
local control_mod = import "spng/control.jsonnet";

/// Top-level arguments:
///
/// input :: file name providing depos
/// output :: if given, a file name with '%d' to receive anode ID number
function(input, output="", detname='pdhd', tpcids=[], engine='Pgrapher', device='cpu', verbosity=0)
    
    local controls = control_mod(device=device, verbosity=verbosity);
    local det = detconf.get(detname, tpcids);

    // make source node
    local source = io.depo_source(input);

    // make sink node
    local onames = if output == ""
                   then ["out"+std.toString(tpc.ident) for tpc in det.tpcs]
                   else [output%tpc.ident for tpc in det.tpcs];
    local sinkf = if output == ""
                  then io.frame_null_sink
                  else io.frame_array_sink;
    local sinks =[sinkf(name, digitize=true) for name in onames];
    local sink = pg.crossline(sinks);

    local splat = det_mod(det, controls.config).depos_to_splat;
    local graph = pg.shuntlines([source, splat, sink]);

    pg.main(graph, app=engine, plugins=["WireCellSpng"], uses=controls.uses)

    
