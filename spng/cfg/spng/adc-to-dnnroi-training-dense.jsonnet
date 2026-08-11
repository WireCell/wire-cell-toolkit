/// Process ADC frame files to dnnroi-fodder (DNNROI training input), dense only.
///
/// This is a single-TPC job producing only the "fodder" tier of
/// dnnroi-training-dense-only.jsonnet.  Unlike that job, it starts from an ADC
/// level frame file instead of from depos.  As there is thus no simulation there
/// is also no "truth" tier -- pair this with a separate splat job if truth is
/// wanted.
///
/// Example usage.
///
/// Simulated ADC, untagged traces:
///
/// wire-cell spng/cfg/spng/adc-to-dnnroi-training-dense.jsonnet \
///   -A input=cosmicsadc.npz -A output=fodder.npz -A tpcid=0
///
/// LArSoft data dump with traces tagged "raw2" holding the channels of APA 2:
///
/// wire-cell spng/cfg/spng/adc-to-dnnroi-training-dense.jsonnet \
///   -A input=data_events.npz -A output=fodder.npz -A intag=raw2 -A tpcid=2
///
/// The output holds one "dense" (eg looseLF) image per view of shape
/// (nchan, ntick) rebinned by "rebin" and scaled down by "scale".
///
/// To visualize, use: https://github.com/brettviren/teepeesee

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local subgraphs_js = import "spng/subgraphs.jsonnet";
local control_js = import "spng/control.jsonnet";
local frame_js = import "spng/frame.jsonnet";
local hacks = import "spng/hacks.jsonnet";
local detconf = import "spng/detconf.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "frame file" format, usually .npz, holding ADC.
// @param output The output file name.
// @param schema The output file schema ("tensor" or "frame")
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc.
// @param rebin The downsample factor along the time dimension.
// @param scale The amount to scale down the dense arrays (looseLF)
// @param tpcid The TPC ID number.
// @param intag The trace tag carried by the input frame file, default "" for
//        untagged.  LArSoft ADC dumps typically use "raw2".
// @param dump The list of intermediates to dump to file.
// @param verbosity The verbosity level for additional logging.
//
// The only required TLA is "input".
//
// Notes:
//
// - The input ADC must hold the channels of the given tpcid.  Eg a pdhd LArSoft
//   dump of channels 5120-7679 requires tpcid=2.  A mismatch shows up as a
//   "ncols=0" in the FrameToTdm debug log.
//
// - intag must name the trace tag in the input file.  FrameFileSource excludes
//   any tagged array whose tag it is not told about and FrameToTdm converts only
//   traces matching the tag of its rule.
//
// - A schema of "frame" may not work.
function(input,
         output="dnnroi-training-fodder-dense.npz",
         schema="tensor",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
         rebin=4,
         scale=4000.0,
         tpcid=3,
         intag="",
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


    // Source of ADC frames.  This replaces the depo source, drifter and
    // inducer of the depo based training jobs.
    local source = io.frame_array_source(
        input, tags=if std.length(intag) == 0 then [] else [intag]);

    local sg = subgraphs_js(tpc, control, pg=wpg);

    local to_tdm = sg.frame_to_tdm(tag=intag);

    local dnnroi_pre = sg.dnnroi_dense_training_preface(
        rebin=irebin, scale=1.0/fscale, tag=intag, extra_name="_preface");

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
        for view in [0,1,2]
        for feat in wc.enumerate(["dense"])
    ]);

    local training_pre = pg.shuntlines([
        dnnroi_pre,
        final_metadata,
        sg.tensor_packer(multiplicity=3, extra_name="_fodder")
    ]);

    local fodder_body = {
        frame: pg.pipeline([
            sg.wrap_bypass(training_pre),
            // Warning: this does not yet work
            sg.tdm_to_frame(traces_tag=fodder_tag, chid_tag=sg.tdm_tagpath(intag)),
            io.frame_sink(output)
        ]),
        tensor: pg.pipeline([
            training_pre,
            io.ttensors_sink(output,
                             include_rules=[], exclude_rules=[],
                             datapath_pattern="tensorsets/{ident}", gzip=wc.numberify(gzip))
        ]),
    }[schema];


    local graph = pg.pipeline([source, to_tdm, fodder_body]);

    // fixme: strictly, only need HIO if saving to HDF.
    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellHio"],
                           uses=controls.uses);
    result
