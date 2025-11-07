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


function (anode, ts, prefix="dnnroi", output_scale=1.0) 
    local apaid = anode.data.ident;
    local prename = prefix + std.toString(apaid);
    local intags = ['loose_lf%d'%apaid, 'mp2_roi%d'%apaid,
                     'mp3_roi%d'%apaid];

    local dnnroi_u = pg.pnode({
        type: "DNNROIFinding",
        name: prename+"u",
        data: {
            anode: wc.tn(anode),
            plane: 0,
            intags: intags,
            decon_charge_tag: "decon_charge%d" %apaid,
            outtag: "dnnsp%du"%apaid,
            output_scale: output_scale,
            forward: wc.tn(ts)
        }
    }, nin=1, nout=1, uses=[ts, anode]);
    local dnnroi_v = pg.pnode({
        type: "DNNROIFinding",
        name: prename+"v",
        data: {
            anode: wc.tn(anode),
            plane: 1,
            intags: intags,
            decon_charge_tag: "decon_charge%d" %apaid,
            outtag: "dnnsp%dv"%apaid,
            output_scale: output_scale,
            forward: wc.tn(ts)
        }
    }, nin=1, nout=1, uses=[ts, anode]);
    local dnnroi_w = pg.pnode({
        type: "PlaneSelector",
        name: prename+"w",
        data: {
            anode: wc.tn(anode),
            plane: 2,
            tags: ["gauss%d"%apaid],
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
                frame: {".*": "dnnsp%d" % apaid},
                // trace: {".*": "dnnsp%d" % apaid},
            } for plane in ["u", "v", "w"]]
        },
    }, nin=3, nout=1);
    
    pg.intern(innodes=[dnnfanout],
              outnodes=[dnnfanin],
              centernodes=dnnpipes,
              edges=[pg.edge(dnnfanout, dnnpipes[ind], ind, 0) for ind in [0,1,2]] +
              [pg.edge(dnnpipes[ind], dnnfanin, 0, ind) for ind in [0,1,2]])
