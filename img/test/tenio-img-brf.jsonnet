// Read in cluster graph, reframe it, write out the frame.

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         infiles="clusters-img-%(anode)s.npz",
         outfiles="frames-brf-%(anode)s.npz",
         anode_iota=null)

    local params = high.params(detector, variant);  
    local mid = high.api(detector, params, options={sparse:false});

    local anodes = mid.anodes();
    local iota = if std.type(anode_iota) == "null" then std.range(0, std.length(anodes)-1) else anode_iota;

    local components = [
        local anode = anodes[aid];
        local said = std.toString(anode.data.ident);
        local acfg={anode: said};

        pg.pipeline([
            high.fio.cluster_file_source(std.format(infiles, acfg), anodes),

            pg.pnode({
                type:"BlobDeclustering",
                name:anode.data.ident,
            }, nin=1,nout=1),

            pg.pnode({
                type:"BlobSetReframer",
                name:anode.data.ident,
                data: {
                    frame_tag: "brf" + said,
                    measure: "value",
                    // source: "blob",
                    source: "activity",
                },
            }, nin=1,nout=1),

            high.fio.frame_file_sink(std.format(outfiles, acfg),digitize=false,
                                     dense={
                                         chbeg:2560*anode.data.ident,
                                         chend:2560*(anode.data.ident+1),
                                         tbbeg:0, tbend:6000}),

        ]) for aid in iota];

    local graph = pg.components(components);
    // local executor = "TbbFlow";
    local executor = "Pgrapher";
    high.main(graph, executor)

    
