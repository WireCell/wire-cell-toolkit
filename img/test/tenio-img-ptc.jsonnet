// Read in cluster graph, produce point cloud, write out.

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         infiles="clusters-img-%(anode)s.npz",
         outfiles="pointcloud-img-%(anode)s.npz",
         anode_iota=null)

    local params = high.params(detector, variant);  
    local mid = high.api(detector, params, options={sparse:false});

    local anodes = mid.anodes();
    local iota = if std.type(anode_iota) == "null" then std.range(0, std.length(anodes)-1) else anode_iota;

    local components = [
        local anode = anodes[aid];
        local acfg={anode: anode.data.ident};

        // Note, the "sampler" must be unique to the "sampling".
        local bs = {
            type: "BlobSampler",
            name: anode.data.ident, 
            data: {
                strategy: [
                    "center","corner","edge","bounds","stepped",
                    {name:"grid", step:1, planes:[0,1]},
                    {name:"grid", step:1, planes:[1,2]},
                    {name:"grid", step:1, planes:[2,0]},
                    {name:"grid", step:2, planes:[0,1]},
                    {name:"grid", step:2, planes:[1,2]},
                    {name:"grid", step:2, planes:[2,0]},
                ],
                extra: [".*"] // want all the extra
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
                    samplers: {samples: wc.tn(bs)},
                },
            }, nin=1, nout=1, uses=[bs]),

            // output
            high.fio.tensor_file_sink(std.format(outfiles, acfg)),


        ]) for aid in iota];

    local graph = pg.components(components);
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

    
