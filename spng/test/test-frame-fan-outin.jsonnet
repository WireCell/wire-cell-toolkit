///  A test job converting an IFrame to tensors and back.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local omnimap = import "omni/all.jsonnet";
local io = import "fileio.jsonnet";


/// Note, different input may require selecting a different anodeid to get any activity.
function(name="pdhd", input="test/data/muon-depos.npz", output="test-frame-fan-outin.npz", paradigm="orig", anodeid="3", device="cpu")
    local source = io.depo_file_source(input);
    local omni = omnimap[name];
    local anode = omni.anodes[std.parseInt(anodeid)];
    local drift = omni.drift(anode);
    local signal = omni.signal(anode);
    local noise = omni.noise(anode);
    local digitize = omni.digitize(anode);
    local verbose = 2;

    local groups = [
        [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
         wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
         wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
        
        [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
    ];
    local group_iota = std.range(0, std.length(groups)-1);

    local all_fans = {
        local ngroups = std.length(groups),

        // This option runs the initial mainline, non-TDM SPNG nodes
        orig: {
            multiplicity: ngroups,
            fanout : pg.pnode({
                type: 'FrameToTorchSetFanout',
                name: "fanout",
                data: {
                    anode: wc.tn(anode),
                    expected_nticks: 6000,
                    output_groups: groups,
                    unsqueeze_output: true,
                },
            }, nin=1, nout=ngroups, uses=[anode]),
            fanin : pg.pnode({
                type: 'TorchTensorSetToFrameFanin',
                name: "",
                data: {
                    anode: wc.tn(anode),
                    input_groups: groups,
                },
            }, nin=ngroups, nout=1, uses=[anode]),

        },

        // This option runs the TDM-compliant nodes
        tdm: {
            multiplicity: ngroups,
            local frametotdm = pg.pnode({
                type: 'SPNGFrameToTdm',
                name: "totdm",
                data: {
                    anode: wc.tn(anode),
                    verbose: verbose,
                    rules: [{
                        tag: "",
                        groups: [{
                            wpids: groups[g],
                        }, for g in group_iota]}],
                },
            }, nin=1, nout=1, uses=[anode]),
            local fanout = pg.pnode({
                type: 'SPNGFanoutNode',
                name: "fanout",
                data: {
                    multiplicity: ngroups,
                    verbose: verbose,
                },
            }, nin=1, nout=ngroups),
            local pipes = [
                pg.pnode({
                    type: 'SPNGFunctionNode',
                    name: "group%s"%g,
                    data: {
                        verbose: verbose,
                        tensor_selection: [
                            { accept: "/frames/\\d+/frame" },
                            { accept: "/frames/\\d+/tags/null/rules/0/groups/%d/.*" % g },
                        ],
                        keep_unselected: false,
                        select_parents: false,
                        combine_policy: "produced_only",
                    },
                }, nin=1, nout=1, uses=[]) for g in group_iota],
            
            fanout: pg.pipeline([frametotdm, 
                                 pg.intern(innodes=[fanout], outnodes=pipes, edges=[
                                     pg.edge(fanout, pipes[n], n, 0) for n in std.range(0, ngroups-1)]),
            ]),
            fanin: pg.pipeline([
                pg.pnode({
                    type: 'SPNGFaninNode',
                    name: "fanin",
                    data: {
                        multiplicity: ngroups,
                        verbose: verbose,
                    },
                }, nin=ngroups, nout=1),
                pg.pnode({
                    type: 'SPNGTdmToFrame',
                    name: "fromtdm",
                    data: {
                        verbose: verbose,
                    },
                }, nin=1, nout=1, uses=[]),
            ]),
        },
    };
    local fans = all_fans[paradigm];
    local body = pg.intern(innodes=[fans.fanout],
                           outnodes=[fans.fanin],
                           edges=[pg.edge(fans.fanout, fans.fanin, n, n)
                                  for n in std.range(0, fans.multiplicity-1)]);
    //local sink = io.frame_file_sink(output);
    local sink = pg.pnode({ type: "DumpFrames" }, nin=1, nout=0);
    local graph = pg.pipeline([source, drift, signal, noise, digitize, body, sink]);
    pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
