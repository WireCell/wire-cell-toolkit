/// Process depo files to dnnroi-fodder (DNNROI training input).
///
/// This is a single-TPC job producing only the "fodder" tier from dnnroi-training.jsonnet.
///
/// Example usage:
/// wire-cell spng/cfg/spng/dnnroi-training-fodder.jsonnet -A input=depos.npz -A output=fodder.npz
///
/// The output contains input to DNNROI of shape (3, nchan, ntick) where the feature
/// dimension of size 3 holds the "dense" (eg looseLF) and the mp2 and mp3 images.
/// The dense image is scaled by 4000 and the mp2/mp3 are boolean.
///
/// To visualize, use: https://github.com/brettviren/teepeesee

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local det_js = import "spng/det.jsonnet";
local drift_js = import "spng/drift.jsonnet";
local subgraphs_js = import "spng/subgraphs.jsonnet";
local control_js = import "spng/control.jsonnet";
local frame_js = import "spng/frame.jsonnet";
local hacks = import "spng/hacks.jsonnet";
local detconf = import "spng/detconf.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param output The output file name.
// @param schema The output file schema ("tensor" or "frame")
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc.
// @param rebin The downsample factor along the time dimension.
// @param scale The amount to scale down the dense arrays (looseLF)
// @param tpcid The TPC ID number.
// @param dump The list of intermediates to dump to file.
// @param verbosity The verbosity level for additional logging.
//
// The only required TLA is "input".
//
// Notes:
//
// - Input depos should be arranged to populate the given tpcid.
//
// - A schema of "frame" may not work.
function(input,
         output="dnnroi-training-fodder.npz",
         schema="tensor",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
         rebin=4,
         scale=4000.0,
         tpcid=3,
         /// U,V,W set to 1 if it has mp2/mp3 produced.
         crossed_views=[1,1,0], //TLA Code, not tla str --> provide via "--tla-code"
         dump="",
         gzip=1,
         verbosity=0)

    
    local cross_view_ids = [vi.index for vi in wc.enumerate(crossed_views) if vi.value == 1];
    local ncrossed = wc.sum(crossed_views);

    # Assure numbers
    local irebin = wc.numberify(rebin);
    local fscale = wc.numberify(scale);
    local itpcid = wc.numberify(tpcid);
    local iverbosity = wc.numberify(verbosity);


    local controls = control_js(device=device, verbosity=iverbosity);
    local control = controls.config;

    // We focus here on just one TPC
    local det = detconf.get(detname, [itpcid]);
    local tpc = det.tpcs[0];

    local frame = frame_js(control);

    ///////////for debugging
    // Wrap a tensor producing pnode with a pickle dumper.
    local dump_tensors(tier, inode, pnode) =
        local wrap="tensors-%(tier)s-%(tpcid)s-%(itype)s-%(iname)s.pkl";
        local name = tier + '_' + inode.type + '_' + inode.name;
        local filename = wrap % {tier: tier, tpcid: std.toString(tpcid), itype:inode.type, iname:inode.name};
        local sink = tio.pickle_tensor_set(filename);
        local tap = frame.sink_taps(name, std.length(pnode.oports), sink);
        pg.shuntline(pnode, tap);

    // Wrap pnode in file dumper for tier if user has tier in dump
    local dump_tensors_maybe(tier, inode, pnode) =
        if std.findSubstr(tier,dump) == []
        then pnode
        else dump_tensors(tier, inode, pnode);

    local have_tier(tier, labels, inode) =
        std.findSubstr(tier, dump) != [] && std.findSubstr(labels[tier], inode.name) != [];

    // Wrap pnode if a key of labels is in dump list and value matches inode.name
    local dump_tensors_matched(labels, inode, pnode) =
        local tiers = [tier
                       for tier in std.objectFields(labels)
                       if have_tier(tier, labels, inode)];
        if tiers == []
        then pnode
        else dump_tensors(tiers[0], inode, pnode);

    local pnode_wrappers = {
        SPNGRebinner: function(inode, pnode)
            dump_tensors_maybe("rebin", inode, pnode),
        SPNGResampler: function(inode, pnode)
            dump_tensors_matched({splat:"splat"}, inode, pnode),
        SPNGTransform: function(inode, pnode)
            dump_tensors_matched({dscale:"dnnroi_scale"}, inode, pnode),
        SPNGKernelConvolve: function(inode, pnode)
            dump_tensors_matched({dnnroi:"dnnroi",wiener:"wiener",decon:"_group"}, inode, pnode),
        SPNGThreshold: function(inode, pnode)
            dump_tensors_matched({wthresh:"wthresh"}, inode, pnode),
        SPNGReduce: function(inode, pnode)
            dump_tensors_matched({stack:"stack"}, inode, pnode),
    };
    ///////////end for debugging
    local wpg = hacks.wrap_pnode(pnode_wrappers);


    // Source of depos
    local source = io.depo_source(input);

    // Common drifter
    local drift = drift_js(det, control).drifter;

    local upstream = pg.pipeline([source, drift]);


    local sg = subgraphs_js(tpc, control, pg=wpg);

    local detmod = det_js(det, control);

    local sim = detmod.inducer;

    local to_tdm = sg.frame_to_tdm();
    local dnnroi_pre = sg.dnnroi_training_preface(crossed_views, rebin, extra_name="_preface");

    // This is a point of collusion between final metadata and the tdm to frame conversion.
    local fodder_tag = "fodder";

    // We have to have a little subgraph just to get the packed tensors into a
    // form that the TdmToFrame can consume.
    local final_metadata = pg.crossline([
        pg.pnode({
            type: 'SPNGTransform',
            name: tpc.name + 'v'+std.toString(view) + 'f' + std.toString(feat.index) + "_final_metadata",
            data: {
                operation: "noop",
                tag: fodder_tag,
                datapath_format: "/frames/{ident}/tags/{tag}/view/%(view)d/feature/%(feat)s/traces"
                                 % {view:view, feat:feat.value},
            }
        }, nin=1, nout=1)
        for view in cross_view_ids
        for feat in wc.enumerate(["dense", "mp2", "mp3"])
    ]);

    local training_pre = pg.shuntlines([
        // 2 tensors each with dim=-3 of size 3.
        sg.expand_dimension(dnnroi_pre, "_fodder"),
        final_metadata,
        sg.tensor_packer(multiplicity=ncrossed*3, extra_name="_fodder")
    ]);

    local fodder = {
        frame: pg.pipeline([
            sg.wrap_bypass(training_pre),
            sg.tdm_to_frame(),  // Warning: this does not yet work
            io.frame_sink(output)
        ]),
        tensor: pg.pipeline([
            training_pre,
            io.ttensors_sink(output,
                             include_rules=[], exclude_rules=[],
                             datapath_pattern="tensorsets/{ident}", gzip=wc.numberify(gzip))
        ]),
    }[schema];


    local graph = pg.pipeline([upstream, sim, to_tdm, fodder]);

    // fixme: strictly, only need HIO if saving to HDF.
    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellGen", "WireCellHio"],
                           uses=controls.uses);
    result


