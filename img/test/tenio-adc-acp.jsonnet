// round trip

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         infiles="frames-adc-%(anode)s.npz",
         outfiles="frames-acp-%(anode)s.npz",
         frame_mode="sparse")

    local mid = high.mid(detector, variant, options={sparse:false});

    local anodes = mid.anodes();
    local nanodes = std.length(anodes);
    local anode_iota = std.range(0, nanodes-1);

    local components = [
        local anode = anodes[aid];
        local acfg={anode: anode.data.ident};
        pg.pipeline([
            high.fio.frame_tensor_file_source(std.format(infiles, acfg)),
                
            high.fio.frame_tensor_file_sink(std.format(outfiles, acfg),digitize=true),

        ]) for aid in anode_iota];

    local graph = pg.components(components);
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

