local wc = import "wirecell.jsonnet";
local real_pg = import "pgraph.jsonnet";
local util = import "spng/util.jsonnet";
local fans_mod = import "spng/fans.jsonnet";

function(control={}, pg)
{
    local fans = fans_mod(control),

    /// Apply Threshold on each view making a 3->3 subgraph
    threshold_views(tpc, extra_name=""):: pg.crossline([pg.pnode(
        {
            type: 'SPNGThreshold',
            name: tpc.name + 'v' + std.toString(it.index) + extra_name,
            data: it.value + control,
        }, nin=1, nout=1) for it in wc.enumerate(tpc.crossview_thresholds)]),

    /// One CrossViews node
    crossview(tpc, view_index, extra_name="")::
        pg.pnode({
            type: 'SPNGCrossViews',
            name: tpc.name + 'v' + std.toString(view_index) + extra_name,
            data: {
                anode: wc.tn(tpc.anode),
                target_index: view_index,
                face_idents: tpc.faces,
                multiplicity: 3,
            } + control,
        }, nin=3, nout=1, uses=[tpc.anode]),



    /// 3->ncrossed for ncrossed > 1.
    crossfan_many(tpc, view_crossed=[1,1,0], extra_name="")::
        local ncrosses = wc.sum(view_crossed);
        local fanouts = [       // 3
            fans.fanout(tpc.name+'crossfan_v'+std.toString(view_index),
                        ncrosses)
            for view_index in [0,1,2] ];
        local crossviews = [    // 3, sparse
            if view_crossed[view_index] == 1
            then $.crossview(tpc, view_index, extra_name=extra_name)
            else null
            for view_index in [0,1,2] ];
        // less than or equal to 3
        local crossviews_only = [c for c in crossviews if std.type(c) != "null"];
        pg.intern(innodes=fanouts,
                  centernodes=crossviews_only,
                  oports=[
                      if view_crossed[view_index] == 1
                      then crossviews[view_index].oports[0]
                      //else fanouts[view_index].oports[ncrosses]
                      for view_index in [0,1,2]
                  ],
                  edges=[
                      pg.edge(f.value, c.value, c.index, f.index)
                      for f in wc.enumerate(fanouts)
                      for c in wc.enumerate(crossviews_only)
                  ]),


    /// A 3->ncrossed with fans and crossviews in the middle.
    /// Connect a partial crossfan based on marking each view as crossed (1) or
    /// not (0) in view_crossed.
    crossfan(tpc, view_crossed=[1,1,0], extra_name="")::
        local ncrosses = wc.sum(view_crossed);
        if ncrosses == 1        // special case, no explicit fanouts
        then $.crossfan(tpc, [vi.index for vi in wc.enumerate(view_crossed) if vi.value == 1][0],
                            extra_name=extra_name)
        else $.crossfan_many(tpc, view_crossed, extra_name),
    

    /// Build a full cross views block producing a 3->3 subgraph with fully connected inner.
    crossfanfull(tpc, extra_name="")::
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
    


    /// An N+1 -> 1 fanin that inputs N "dense" inputs and one crossviews for
    /// one view and outputs a stacked tensor consisting of N + Ncv where Ncv
    /// are the number of MP numbers to extract from the crossviews.  The N
    /// dense inputs should be connected to the first N iports and the
    /// crossviews to the last iport.
    extract_stack_one(tpc, view_index, ndense=1, mps=["mp2", "mp3"])::
        local this_name = tpc.name + "v" + std.toString(view_index);
        local nex = std.length(mps);
        local ex = pg.pnode({
            type: 'SPNGCrossViewsExtract',
            name: this_name,
            data: {
                extraction: mps
            },
        }, nin=1, nout=nex);
        local op = util.reduce_one("stack", dim=-3, multiplicity=ndense+nex, name=this_name+"_stack");
        pg.intern(outnodes=[op], // single output
                  centernodes=[ex],
                  iports=[      // ndense stack iports exposed followed by extract iports
                      op.iports[ind]
                      for ind in wc.iota(ndense)
                  ] + ex.iports,
                  edges=[       // connect extract oports to last stack iports
                      pg.edge(ex, op, wind, wind+ndense)
                      for wind in wc.iota(nex)
                  ]),
    
    /// Connect 'dense' and crossviews ('cvs') each a source of 3, per-view
    /// oports with the cross views extract ('exs') with 3, per-view iports.
    /// Returns effectively a source with 3, per-view oports.
    ///
    /// This idea is broken because it assumes all views have dense and
    /// crossviews while W tends to skip this.
    connect_extract_stack(tpc, dense, cvs, exs, extra_name="")::
        pg.intern(name=extra_name,
                  oports=exs,
                  edges=[
                      pg.edge(dense, exs[view_index], view_index, 0)
                      for view_index in [0,1,2]
                  ] + [
                      pg.edge(cvs, exs[view_index], view_index, 1)
                      for view_index in [0,1,2]
                  ]),
}
