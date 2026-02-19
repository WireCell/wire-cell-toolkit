/// Process depo files to dnnroi-truth (ground truth ROI masks).
///
/// This is a single-TPC job producing only the "truth" tier from dnnroi-training.jsonnet.
///
/// Example usage:
/// wire-cell spng/cfg/spng/dnnroi-training-truth.jsonnet -A input=depos.npz -A output=truth.npz
///
/// The output contains ROI masks of shape (nchan, ntick) from depo flux splat.
///
/// To visualize, use: https://github.com/brettviren/teepeesee

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
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
// @param scale The amount to scale down the dense arrays (splat)
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
         output="dnnroi-training-truth.npz",
         schema="tensor",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
         rebin=4,
         scale=4000.0,
         tpcid=3,
         dump="",
         gzip=1,
         verbosity=0)

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
            dump_tensors_matched({sscale:"splat_scale"}, inode, pnode),
    };
    ///////////end for debugging
    local wpg = hacks.wrap_pnode(pnode_wrappers);


    // Source of depos
    local source = io.depo_source(input);

    // Common drifter
    local drift = drift_js(det, control).drifter;

    local upstream = pg.pipeline([source, drift]);


    local sg = subgraphs_js(tpc, control, pg=wpg);

    local truth = {
        frame: pg.pipeline([
            sg.splat_frame(ratio=1.0/irebin, scale=1.0/fscale, extra_name="_splat"),
            // sg.splat_cellview_frame(ratio=1.0/irebin, scale=1.0/fscale, extra_name="_splat"),
            io.frame_sink(output)
        ]),
        tensor: pg.pipeline([
            // sg.splat_cellview_tensor(ratio=1.0/irebin, scale=1.0/fscale, extra_name="_splat"),
            sg.splat_tensor(ratio=1.0/irebin, scale=1.0/fscale, extra_name="_splat"),
            io.ttensors_sink(output, gzip=wc.numberify(gzip), control=control),
        ]),
    }[schema];


    local graph = pg.pipeline([upstream, truth]);

    // fixme: strictly, only need HIO if saving to HDF.
    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellGen", "WireCellHio"],
                           uses=controls.uses);
    result


