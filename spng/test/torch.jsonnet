// TODO -- brief descrip

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';


// local wpidU0 = wc.WirePlaneId(wc.Ulayer, 0);
// local wpidV0 = wc.WirePlaneId(wc.Vlayer, 0);
// local wpidW0 = wc.WirePlaneId(wc.Wlayer, 0);

// local wpidU1 = wc.WirePlaneId(wc.Ulayer, 1);
// local wpidV1 = wc.WirePlaneId(wc.Vlayer, 1);
// local wpidW1 = wc.WirePlaneId(wc.Wlayer, 1);

// function make_wpid(tools)

function(tools, override = {}) {
    
    make_fanout(anode, name=null):: g.pnode({
        type: 'FrameToTorchSetFanout',
        name:
            if std.type(name) == 'null'
            then anode.name + 'torchfanout%d' % anode.data.ident
            else name,
        data: {
            anode: wc.tn(anode),
            expected_nticks: 6000,

            output_groups: [
                [wc.WirePlaneId(wc.Ulayer, 0, anode.data.ident),
                 wc.WirePlaneId(wc.Ulayer, 1, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Vlayer, 0, anode.data.ident),
                 wc.WirePlaneId(wc.Vlayer, 1, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Wlayer, 0, anode.data.ident)],
                
                [wc.WirePlaneId(wc.Wlayer, 1, anode.data.ident)],
            ],

        }
    }, nin=1, nout=4, uses=[anode]),

}