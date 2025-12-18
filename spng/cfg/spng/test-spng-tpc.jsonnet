// This covers SPNG from ADC input to signal output.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";

local frame_mod = import "spng/frame.jsonnet";
local spng_mod = import "spng/spng.jsonnet";

local detconf = import "spng/detconf.jsonnet";
local detector = import "spng/detector.jsonnet";
local control_module = import "spng/control.jsonnet";



// An experiment to inject a generic wrapper.
// FIXME: move this into its own jsonnet
local wrapped_pnode(class_wrappers)=
    function(inode, nin=0, nout=0, uses=[], name=null) // pnode
        //local itype=std.trace('inode.type=%s'%inode.type, inode.type);
        local itype=inode.type;
        local w = std.get(class_wrappers, itype, function(inode, pnode) pnode);
        w(inode, pg.pnode(inode, nin=nin, nout=nout, uses=uses, name=name));


/// Top-level arguments:
///
/// input :: name of input file
///
/// output :: name for output file
///
/// view_crossed :: which views should be subject to crossviews
///
/// dump :: if given, a list of tiers to dump.  See pipeline definition below for labels
function(input, output, tpcid=0, view_crossed=[1,1,0],
         dump="",
         detname='pdhd',
         // fixme: add seeds=[1,2,3,4]
         engine='Pgrapher', device='cpu', verbosity=0, semaphore=1)

    local tpcid_num = wc.intify(tpcid);

    local det = detector.subset(detconf[detname], [tpcid_num]);
    local tpc = det.tpcs[tpcid_num];

    local controls = control_module(device=device,
                                    verbosity=wc.intify(verbosity),
                                    semaphore=wc.intify(semaphore));
    local control = controls.config;

    local frame_nodes = frame_mod(control);

    local dump_views_maybe(tier, pnode) =
        local wrap="tensors-%(tier)s-%(tpcid)s.pkl";
        local nodump = std.length(std.findSubstr(tier, dump)) == 0;
        if nodump
        then pnode
        else
            local name = 'wrap_tpc' + std.toString(tpcid) + '_' + tier;
            local filename = wrap % {tier: tier, tpcid: std.toString(tpcid)};
            local sink = tio.pickle_tensor_set(filename);
            local tap = frame_nodes.sink_taps(name, std.length(pnode.oports), sink);
            pg.shuntline(pnode, tap);
            

    local dump_inode_maybe(tier, inode, pnode) =
        local wrap="tensors-%(tier)s-%(tpcid)s-%(itype)s-%(iname)s.pkl";
        local nodump = std.length(std.findSubstr(tier, dump)) == 0;
        if nodump
        then pnode
        else
            local name = tier + '_' + inode.type + '_' + inode.name;
            local filename = wrap % {tier: tier, tpcid: std.toString(tpcid), itype:inode.type, iname:inode.name};
            local sink = tio.pickle_tensor_set(filename);
            local tap = frame_nodes.sink_taps(name, std.length(pnode.oports), sink);
            pg.shuntline(pnode, tap);

    local pnode_wrappers = {
        SPNGCrossViews: function(inode, pnode) dump_inode_maybe("crossviews", inode, pnode),
        SPNGCrossViewsExtract: function(inode, pnode) dump_inode_maybe("mps", inode, pnode),
        SPNGTensorForward: function(inode, pnode) dump_inode_maybe("dnnroi", inode, pnode),
    };
    local wpg = pg + { pnode:: wrapped_pnode(pnode_wrappers) };
            

    local source = io.frame_array_source(input);
    local spng = spng_mod(tpc, control, view_crossed, wpg);
    local sink = io.frame_array_sink(output);

    local graph = pg.pipeline([
        source,
        spng,
        sink
    ]);

    // FIXME: still need to bring original frame tensor set to output to form a
    // full frame tensor set and then feed to TdmToFrame

    pg.main(graph, engine, plugins=["WireCellSpng"], uses=controls.uses)


