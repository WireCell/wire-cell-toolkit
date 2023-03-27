// Read in cluster graph, produce point cloud, write out.

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         infiles="clusters-img-%(anode)s.npz",
         outfiles="pointcloud-img-%(anode)s.npz")

    local mid = high.mid(detector, variant, options={sparse:false});

    local anodes = mid.anodes();
    local nanodes = std.length(anodes);
    local anode_iota = std.range(0, nanodes-1);


    local components = [
        local anode = anodes[aid];
        local acfg={anode: anode.data.ident};
        local bs = {
            type: "BlobSampler",
            name: anode.data.ident,
            data: {
                strategy: ["center","corner","edge"]
                #strategy: "center"
            }};

        pg.pipeline([
            high.fio.cluster_file_source(std.format(infiles, acfg), anodes),

            pg.pnode({
                type:"BlobDeclustering",
                name:anode.data.ident,
            }, nin=1,nout=1),

            pg.pnode({
                type:"BlobSampling",
                name:anode.data.ident,
                data: {
                    sampler: wc.tn(bs),
                },
            }, nin=1, nout=1, uses=[bs]),

            // output
            high.fio.tensor_file_sink(std.format(outfiles, acfg)),


        ]) for aid in anode_iota];

    local graph = pg.components(components);
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

    
