local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         infiles="frames-sig-%(anode)s.npz",
         outfiles="clusters-img-%(anode)s.npz")

    local mid = high.mid(detector, variant, options={sparse:false});

    local anodes = mid.anodes();
    local nanodes = std.length(anodes);
    local anode_iota = std.range(0, nanodes-1);

    local components = [
        local anode = anodes[aid];
        local acfg={anode: anode.data.ident};
        pg.pipeline([
            high.fio.frame_tensor_file_source(std.format(infiles, acfg)),

            mid.img(anode),

            high.fio.cluster_file_sink(std.format(outfiles, acfg))

        ]) for aid in anode_iota];

    local graph = pg.components(components);
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

