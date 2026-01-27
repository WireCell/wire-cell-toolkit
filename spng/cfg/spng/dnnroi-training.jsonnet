/// Process depo files to dnnroi-input and dnnroi-truth.
///
/// This is a single-TPC job.
///
/// Example usage:
/// wire-cell spng/cfg/spng/dnnroi-training.jsonnet -A input=depos.npz
///
/// By default it produces dnnroi-training-{truth,fodder}.npz each with contents like:
///
/// $ unzip -v dnnroi-training-truth.npz
/// Archive:  dnnroi-training-truth.npz
///  Length   Method    Size  Cmpr    Date    Time   CRC-32   Name
/// --------  ------  ------- ---- ---------- ----- --------  ----
///        4  Defl:N        9 -125% 01-13-2026 12:31 25cbfc4f  tensorset_0_metadata.json
///      291  Defl:N      175  40% 01-13-2026 12:31 f539e510  tensor_0_0_metadata.json
///  9600128  Defl:N     9418 100% 01-13-2026 12:31 99d05f42  tensor_0_0_array.npy
///      291  Defl:N      180  38% 01-13-2026 12:31 50445c69  tensor_0_1_metadata.json
///  9600128  Defl:N     9418 100% 01-13-2026 12:31 99d05f42  tensor_0_1_array.npy
///      291  Defl:N      180  38% 01-13-2026 12:31 e61e6ce6  tensor_0_2_metadata.json
/// 11520128  Defl:N    11290 100% 01-13-2026 12:31 cab52a81  tensor_0_2_array.npy
/// --------          -------  ---                            -------
/// 30721261            30670 100%                            7 files
///
/// For training DNNROI, only the *_array.npy are useful.  The last letter
/// implies the plane number.  The fodder only has 0 and 2.  Truth are ROI masks
/// of shape (nchan, ntick) from depo flux splat, fodder is a input to DNNROI of
/// shape (3, nchan, ntick) the feature dimension of size 3 holds the "dense"
/// (eg looseLF) and the mp2 and mp3 images.  The dense image is scaled by 4000
/// and the mp2/mp3 are boolean.
///
/// To visualize, use: https://github.com/brettviren/teepeesee

local wc = import "wirecell.jsonnet";
local pg = import 'pgraph.jsonnet';
local io = import "spng/io.jsonnet";
local tio = import "spng/torchio.jsonnet";
local det_js = import "spng/det.jsonnet";
local drift_js = import "spng/drift.jsonnet";
local splatroi_js = import "spng/splatroi.jsonnet";
local subgraphs_js = import "spng/subgraphs.jsonnet";
local control_js = import "spng/control.jsonnet";
local frame_js = import "spng/frame.jsonnet";
local hacks = import "spng/hacks.jsonnet";
local detconf = import "spng/detconf.jsonnet";

// The TLAs:
//
// @param input The name of a file in WCT "depo file" format, usually .npz.
// @param outpat The output file name pattern with format variables.
// @param detname The name of a supported detector, default "pdhd".
// @param engine The name of the graph execution engine, default Pgrapher or TbbFlow.
// @param device The name of the device for SPNG nodes, default "cpu" or "gpu", "gpu1", etc. 
//
// The only required TLA is "input".  
//
// The outpat must include these format variables:
// - %(tier)s will be filled with the label "input" or "truth".
//
// Note, this hard-wires use of TPC ID 0.  Input depos should be arranged to
// populate this TPC.
//
function(input,
         outpat="dnnroi-training-%(tier)s.npz",
         detname='pdhd',
         engine='Pgrapher',
         device='cpu',
         rebin=4,
         scale=4000.0,
         tpcid=3,
         dump="",
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
            dump_tensors_matched({dscale:"dnnroi_scale", sscale:"splat_scale"}, inode, pnode),
        SPNGKernelConvolve: function(inode, pnode)
            dump_tensors_matched({dnnroi:"dnnroi",wiener:"wiener",decon:"_group"}, inode, pnode),
        SPNGThreshold: function(inode, pnode)
            dump_tensors_matched({wthresh:"wthresh"}, inode, pnode),
        SPNGReduce: function(inode, pnode)
            dump_tensors_matched({stack:"stack"}, inode, pnode),
    };
    local wpg = hacks.wrap_pnode(pnode_wrappers);


    // Source of depos
    local source = io.depo_source(input);

    // Common drifter
    local drift = drift_js(det, control).drifter;

    local upstream = pg.pipeline([source, drift]);

    // Fan out depos to the splat and sim+sigproc subgraphs
    local depo_fan = pg.pnode({
        type:'DepoSetFanout',
        name: det.name + "_depo",
        data: {
            multiplicity:2
        }
    },nin=1, nout=2);


    local truth = splatroi_js(tpc, control, rebin=irebin, scale=1.0/fscale, pg=wpg);
    local detmod = det_js(det, control);

    local sim = detmod.inducer;


    //local fodder = roifodder_js(tpc, control, rebin=irebin, scale=1.0/fscale, pg=wpg);
    local sg = subgraphs_js(tpc, control, pg=wpg);
    local to_tdm = sg.frame_to_tdm();
    local fodder = sg.dnnroi_training_preface();
    local simfodder = pg.pipeline([sim, to_tdm, fodder]);

    local body = pg.intern(innodes=[upstream], outnodes=[truth, simfodder],
                           centernodes=[depo_fan],
                           edges=[
                               pg.edge(upstream, depo_fan),
                               pg.edge(depo_fan, truth, 0, 0),
                               pg.edge(depo_fan, simfodder, 1, 0)]);


    /// Connect a subgraph to the source's oports.  Source is interned.
    local multi_sink(source, filename, extra_name="", prefix="") =
        local mult = std.length(source.oports);
        local this_name = tpc.name + extra_name;
        local pack = pg.pnode({
            type: 'SPNGTensorPacker',
            name: this_name,
            data: {
                multiplicity: mult,
            } + controls
        }, nin=mult, nout=1);
        local ttt = pg.pnode({
            type: 'TorchToTensor',
            name: this_name,
            data: {},
        }, nin=1, nout=1);
        local sink = pg.pnode({
            type: 'TensorFileSink',
            name: this_name,
            data: {
                outname: filename,
                prefix: prefix,
            },
        }, nin=1, nout=0);
        pg.intern(centernodes=[source,pack,ttt,sink],
                  edges=[
                      pg.edge(source, pack, port.index, port.index)
                      for port in wc.enumerate(source.oports)
                  ] + [
                      pg.edge(pack, ttt),
                      pg.edge(ttt, sink)
                  ]);
        
    local truth_sink = multi_sink(truth, outpat % {tier:"truth"}, "truth");
    local fodder_sink = multi_sink(simfodder, outpat % {tier:"fodder"}, "fodder");
    

    local graph = pg.components([body, truth_sink, fodder_sink]); 
    //local graph = fodder;
    //local graph = pg.components([sim, fodder]);

    local result = pg.main(graph, app=engine,
                           plugins=["WireCellSpng", "WireCellGen"],
                           uses=controls.uses);
    result

