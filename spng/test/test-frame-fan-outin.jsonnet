///  A test job converting an IFrame to tensors and back.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local omnimap = import "omni/all.jsonnet";
local io = import "fileio.jsonnet";

/// Note, different input may require selecting a different anodeid to get any activity.
function(name="pdhd", input="test/data/muon-depos.npz", output="test-frame-fan-outin.npz", paradigm="orig", anodeid="3")
    local source = io.depo_file_source(input);
    local omni = omnimap[name];
    local anode = omni.anodes[std.parseInt(anodeid)];
    local drift = omni.drift(anode);
    local signal = omni.signal(anode);
    local noise = omni.noise(anode);
    local digitize = omni.digitize(anode);

    local all_fans = {
        orig: {
            local groups = [
                [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
                 wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
                 wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
            ],
            local ngroups = std.length(groups),
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
        tdm: {
        },
    };
    local fans = all_fans[paradigm];
    local body = pg.intern(innodes=[fans.fanout],
                           outnodes=[fans.fanin],
                           edges=[pg.edge(fans.fanout, fans.fanin, n, n)
                                  for n in std.range(0, fans.multiplicity-1)]);
    local sink = io.frame_file_sink(output);
    local graph = pg.pipeline([source, drift, signal, noise, digitize, body, sink]);
    pg.main(graph, 'TbbFlow', plugins=["WireCellSpng"])
