// round trip

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         infiles="frames-adc-%(anode)s.npz",
         outfiles="frames-acp-%(anode)s.npz",
         frame_mode="sparse",anode_iota=null)

    local params = high.params(detector, variant);  
    local mid = high.api(detector, params, options={sparse:false});

    local anodes = mid.anodes();
    local iota = if std.type(anode_iota) == "null" then std.range(0, std.length(anodes)-1) else anode_iota;

    local components = [
        local anode = anodes[aid];
        local acfg={anode: anode.data.ident};
        pg.pipeline([
            high.fio.frame_tensor_file_source(std.format(infiles, acfg)),
                
            high.fio.frame_tensor_file_sink(std.format(outfiles, acfg),digitize=true),

        ]) for aid in iota];

    local graph = pg.components(components);
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

