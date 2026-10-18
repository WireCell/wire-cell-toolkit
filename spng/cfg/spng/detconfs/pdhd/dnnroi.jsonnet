// This produces a function to configure DNN-ROI for one APA given
// anode and torch service (ts) objects.
// 
// The prefix is prepended to all internal node names if uniqueness
// beyond anode ID is needed.  The output_scale allows for an ad-hoc
// scaling of dnnroi output.  The U and W planes will go through
// dnnroi while hte W plane will be shunted.  What comes out will be a
// unified frame with frame tag "dnnspN" where "N" is the anode ID.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";


local detectors = import "detectors.jsonnet";
local detname = "pdhd";
local det = detectors[detname];  // me

// Defaults match pgrapher/experiment/pdhd/dnnroi_pp.jsonnet for the 6-ch
// KD model: rebin 4 (training), pad so the post-rebin length is divisible by
// 2^4 for the UNet stages, 6 input channels with the per-channel
// normalization baked into the .ts (so input_scale 1.0).
function (tpc, prefix="dnnroi", output_scale=1.0, nticks=6000, tick_per_slice=4,
          tick_pad_multiple=null, mask_thresh=0.2, nchan=6, nchunks=1, device="cpu")
    local ts = {
        type: "TorchService",
        name: "dnnroi",
        data: {
            model: "pipe_distill_transformer_6ch.ts", #Set WIRECELL_PATH i.e. to /nfs/data/1/calcuttj/toolkit_testing/dunereco/dunereco/DUNEWireCell/pdhd/
            device: device, // "gpucpu",
            concurrency: 1,
        },
    };
    local anode = tpc.anode;
    local apaid = anode.data.ident;
    local prename = prefix + std.toString(apaid);
    // Order must match the training im_tags.  3-ch models use the first three.
    assert nchan == 3 || nchan == 6 : 'dnnroi nchan must be 3 or 6, got %d' % nchan;
    local intags = ['loose_lf%d'%apaid, 'mp2_roi%d'%apaid,
                    'mp3_roi%d'%apaid, 'tight_lf%d'%apaid,
                    'decon_charge%d'%apaid, 'gauss%d'%apaid][0:nchan];
    local pad_multiple = if tick_pad_multiple == null
                         then tick_per_slice * 16 else tick_pad_multiple;

    local dnn_node(plane, letter) = pg.pnode({
        type: "DNNROIFinding",
        name: prename+letter,
        data: {
            anode: wc.tn(anode),
            plane: plane,
            intags: intags,
            decon_charge_tag: "decon_charge%d" %apaid,
            outtag: "dnnsp%d%s" % [apaid, letter],
            summary_tag: "wiener%d" % apaid,
            // 6-ch models normalize internally; the 3-ch model needs 1/4000.
            input_scale: if nchan == 6 then 1.0 else 1.0 / 4000,
            output_scale: output_scale,
            mask_thresh: mask_thresh,
            forward: wc.tn(ts),
            tick_per_slice: tick_per_slice,
            tick_pad_multiple: pad_multiple,
            nticks: nticks,
            nchunks: nchunks
        }
    }, nin=1, nout=1, uses=[ts, anode]);
    local dnnroi_u = dnn_node(0, "u");
    local dnnroi_v = dnn_node(1, "v");

    local dnnroi_w = pg.pnode({
        type: "PlaneSelector",
        name: prename+"w",
        data: {
            anode: wc.tn(anode),
            plane: 2,
            tags: ["gauss%d"%apaid],
            summary_tags: ["wiener%d"%apaid],
            tag_rules: [{
                frame: {".*":"DNNROIFinding"},
                trace: {["gauss%d"%apaid]:"dnnsp%dw"%apaid},
            }],
        }
    }, nin=1, nout=1, uses=[anode]);

    local dnnpipes = [dnnroi_u, dnnroi_v, dnnroi_w];
    local dnnfanout = pg.pnode({
        type: "FrameFanout",
        name: prename,
        data: {
            multiplicity: 3
        }
    }, nin=1, nout=3);

    local dnnfanin = pg.pnode({
        type: "FrameFanin",
        name: prename,
        data: {
            multiplicity: 3,
            tag_rules: [{
                frame: {".*": "dnnsp%d%s" % [apaid,plane]},
                trace: {".*": "dnnsp%d%s" % [apaid,plane]},
            } for plane in ["u", "v", "w"]]
        },
    }, nin=3, nout=1);

    local retagger = pg.pnode({
        type: "Retagger",
        name: 'dnnroi%d' % apaid,
        data: {
            // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
            tag_rules: [{
                // Retagger also handles "frame" and "trace" like fanin/fanout
                // merge separately all traces like gaussN to gauss.
                frame: {
                    ".*": "dnnsp%d" % apaid
                },
                merge: {
                    ".*": "dnnsp%d" % apaid
                },
                //   merge: {'dnnsp\\d[uvw]' : 'dnnsp%d' %apaid,},
            }],
        },
    }, nin=1, nout=1);
    
    pg.intern(innodes=[dnnfanout],
              outnodes=[retagger],
              centernodes=dnnpipes+[dnnfanin],
              edges=[pg.edge(dnnfanout, dnnpipes[ind], ind, 0) for ind in [0,1,2]] +
              [pg.edge(dnnpipes[ind], dnnfanin, 0, ind) for ind in [0,1,2]] +
              [pg.edge(dnnfanin, retagger, 0, 0)])
