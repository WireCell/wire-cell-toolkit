local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fans = import "fans.jsonnet";

{
    /// Apply Threshold on each view making a 3->3 subgraph
    threshold_views(tpc, extra_name=""):: pg.crossline([pg.pnode(
        {
            type: 'SPNGThreshold',
            name: tpc.name + 'v' + std.toString(it.index),
            data: it.value
        }, nin=1, nout=1) for it in wc.enumerate(tpc.crossview_thresholds)]),

    /// One CrossViews node
    crossview(tpc, view_index, extra_name="")::
        pg.pnode({
            type: 'SPNGCrossViews',
            name: tpc.name + 'v' + std.toString(view_index) + extra_name,
            data: {
                anode: tpc.anode,
                target_index: view_index,
                face_idents: tpc.faces,
            },
        }, nin=3, nout=1, uses=[tpc.anode]),

    /// Build the cross views block producing a 3->3 subgraph with fully connected inner.
    crossfan(tpc, extra_name="")::
        local fanouts = [fans.fanout(tpc.name + 'crossfan_v'+std.toString(view_index), 3)
                         for view_index in [0,1,2]];
        local crossviews = [$.crossview(tpc, view_index, extra_name=extra_name)
                            for view_index in [0,1,2]];
        pg.intern(innodes=fanouts, outnodes=crossviews, 
                  edges=[
                      pg.edge(fanouts[fan_index], crossviews[view_index], view_index, fan_index)
                      for view_index in [0,1,2]
                      for fan_index in [0,1,2]
                  ]),
    


}
