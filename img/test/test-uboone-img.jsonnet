// This runs uboone sim+sigproc+img with a side pipeline of
// BlobDepoFill to produce blobs with "true" charge based on depos.
//
// Use like:
// wire-cell -A depos=depos.npz -A outimg=blobs-img.npz -A outtru=blobs-tru.npz -c test-uboone-img.jsonnet


local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);
local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";
local nf_maker = import "pgrapher/experiment/uboone/nf.jsonnet";
local chndb_maker = import "pgrapher/experiment/uboone/chndb.jsonnet";
local sp_maker = import "pgrapher/experiment/uboone/sp.jsonnet";


local anodes = tools.anodes;
local anode = anodes[0];
local aname = anode.data.ident;

// CHECK THIS FOR MICROBOONE.  ProtoDUNE-SP sim needs 314*wc.us.
local time_offset = 314*wc.us;

local make_random(seeds=[0,1,2,3,4], generator="default") = {
    type: "Random",
    name: generator + "-" +
          std.join('-', std.map(std.toString,seeds)),
    data: {
        generator: generator,
        seeds: seeds,
    }
};

// don't fluctuate to help keep down variations between tests
local make_drifter(random, xregions, lar, fluctuate=false, name="") =
    local d = pg.pnode({
        type: 'Drifter',
        name: name,
        data: lar {
            rng: wc.tn(random),
            xregions: xregions,
            fluctuate: fluctuate,
        },
    }, nin=1, nout=1, uses=[random]);
    pg.pnode({
        type: 'DepoSetDrifter',
        name: name,
        data: { drifter: wc.tn(d) }
    }, nin=1, nout=1, uses=[d]);

local random = make_random();
local xregions = wc.unique_list(std.flattenArrays([v.faces for v in params.det.volumes]));
local drifter = make_drifter(random, xregions, params.lar);

local waveform_map = {
    type: 'WaveformMap',
    name: "",
    data: {
        filename: "microboone-charge-error.json.bz2",
    },
    uses: [],
};

local charge_err = pg.pnode({
    type: 'ChargeErrorFrameEstimator',
    name: "",
    data: {
        intag: "gauss",
        outtag: 'gauss_error',
        anode: wc.tn(anode),
	rebin: 4,  // this number should be consistent with the waveform_map choice
	fudge_factors: [2.31, 2.31, 1.1],  // fudge factors for each plane [0,1,2]
	time_limits: [12, 800],  // the unit of this is in ticks
        errors: wc.tn(waveform_map),
    },
}, nin=1, nout=1, uses=[waveform_map, anode]);


local slicing(tag="", span=4, active_planes=[0,1,2], masked_planes=[], dummy_planes=[]) =
    pg.pnode({
        type: "MaskSlices",
        name: "",
        data: {
            tag: tag,
            tick_span: span,
            anode: wc.tn(anode),
            min_tbin: 0,
            max_tbin: 9592,
            active_planes: active_planes,
            masked_planes: masked_planes,
            dummy_planes: dummy_planes,
        },
    }, nin=1, nout=1, uses=[anode]);


local tiling() = pg.pnode({
    type: "GridTiling",
    name: "",
    data: {
        anode: wc.tn(anode),
        face: 0,
    }
}, nin=1, nout=1, uses=[anode]);

local solving(spans=1.0, threshold=0.0) =
    pg.pipeline([

        pg.pnode({
            type: "BlobClustering",
            name: "",
            data:  { spans : spans }
        }, nin=1, nout=1),

        pg.pnode({
            type: "BlobGrouping",
            name: ""
        }, nin=1, nout=1),

        pg.pnode({
            type: "ChargeSolving",
            name: "",
            data:  {
                weighting_strategies: ["uniform"], //"uniform", "simple"
            }
        }, nin=1, nout=1)
    ], name="");


local blob_sink(outname, fmt="json") = pg.pnode({
    type: "ClusterFileSink",
    name: outname,
    data: {
        outname: outname,
        format: fmt,
    }
}, nin=1, nout=0);



// Catch clusters for depo filling
local make_catcher(outtru) =
    local bsf = pg.pnode({
        type:'ClusterFanout',
        name: "",
        data: {
            multiplicity: 2
        }}, nin=1, nout=2);
    local dbf = pg.pnode({
        type:'BlobDepoFill',
        name:"",
        data: {
            speed: params.lar.drift_speed, // 1.56*wc.mm/wc.us, 
            time_offset: time_offset,
        }}, nin=2, nout=1);
    local cfs = blob_sink(outtru);
    pg.intern([bsf, dbf], centernodes=[cfs],
              edges=[pg.edge(bsf,dbf,1,0),
                     pg.edge(dbf,cfs,0,0)],
              iports=[bsf.iports[0], dbf.iports[1]],
              oports=[bsf.oports[0]],
              name="");

local sim = sim_maker(params, tools);
local chndb = chndb_maker(params, tools).wct("perfect");

local depo_file_source(filename, scale=1.0) =
    pg.pnode({
        type: 'DepoFileSource',
        name: filename,
        data: { inname: filename, scale: scale }
    }, nin=0, nout=1);

local make_depos(depos="depos.npz") =
    local ds = depo_file_source(depos);
    local dsfan = pg.pnode({
        type:'DepoSetFanout',
        name:'',
        data: {
            multiplicity: 2,
        }}, nin=1, nout=2);
    pg.intern(centernodes=[ds,dsfan],
              edges=[pg.edge(ds,dsfan,0,0)],
              oports=dsfan.oports, name="deposource");


local sp = sp_maker(params, tools);

local frametaps(out, tags=[], digitize=false) =
    if std.type(out) == "null" || std.length(out) == 0 then
        []
    else
        local sink = pg.pnode({
            type: 'FrameFileSink',
            name: out,
            data: {
                outname: out,
                tags: tags,
                digitize: digitize,
            },
        }, nin=1, nout=0);
        local fan = pg.pnode({
            type: 'FrameFanout',
            name: out,
            data: { multiplicity: 2 },
        }, nin=1, nout=2);
        [pg.intern(innodes=[fan], outnodes=[fan],
                   centernodes = [sink],
                   edges=[pg.edge(fan, sink, 1, 0)],
                   name=out)];
local depotaps(out) =
    if std.type(out) == "null" || std.length(out) == 0 then
        []
    else
        local sink = pg.pnode({
            type: 'DepoFileSink',
            name: out,
            data: { outname: out, },
        }, nin=1, nout=0);
        local fan = pg.pnode({
            type: 'DepoFanout',
            name: out,
            data: { multiplicity: 2 },
        }, nin=1, nout=2);
        [pg.intern(innodes=[fan], outnodes=[fan],
                   centernodes = [sink],
                   edges=[pg.edge(fan, sink, 1, 0)],
                   name=out)];
                         

// Top level function with TLAs
local make_graph(depos="depos.npz",
                 outimg="blobs-img.npz", outtru="blobs-tru.npz",
                 drifted="", outadc="", outsig="") =
    local ds = make_depos(depos);
    local catcher = make_catcher(outtru);

    local pipeline = [
        ds, 
    ] + depotaps(drifted) + [
        sim.signal_sets,
        sim.add_noise(sim.make_noise_model(anode, sim.miscfg_csdb)),
        sim.digitizer(anode, tag="orig")
    ] + frametaps(outadc, digitize=true) + [
        nf_maker(params, tools, chndb),
        sp,
    ] + frametaps(outsig, digitize=false) + [
        charge_err,
        slicing("gauss", 109),
        tiling(),
        solving(),
        blob_sink(outimg)
    ];

    local pl = pg.pipeline(pipeline, name='mainpipe');
    // Insert the catcher node.
    local pl2 = pg.insert_node(pl, {
        tail: {node: 'BlobGrouping', port: 0},
        head: {node: 'ChargeSolving', port: 0}
    }, catcher, catcher, 0, 0, name='insert_catcher');
    // Connect the ports 1 of depos and catcher
    local graph = pg.intern(centernodes=[pl2],
                            edges=[pg.edge(ds, catcher, 1, 1)],
                            name="graph");
    local engine = 'TbbFlow';
    local app = {
        // type: 'Pgrapher',
        type: engine,
        name: "",
        data: {
            edges: pg.edges(graph),
        },
    };
    local cmdline = {
        type: "wire-cell",
        name: "",
        data: {
            plugins: [
                "WireCellPgraph", "WireCellTbb", "WireCellSio",
                "WireCellGen", "WireCellSigProc", "WireCellImg",
            ],
            apps: [wc.tn(app)]
        }
    };
    [cmdline] + pg.uses(graph) + [app];


make_graph
