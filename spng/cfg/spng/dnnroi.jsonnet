local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fans = import "fans.jsonnet";
local util = import "util.jsonnet";

{

    prepare_pipe(name, view_index, decon, crossview, what=["mp2", "mp3"])::
        local this_name = name + "v" + std.toString(view_index);
        local nex = std.length(what);
        local ex = pg.pnode({
            type: 'SPNGCrossViewsExtract',
            name: this_name,
            data: {
                extraction: what
            },
        }, nin=1, nout=nex);
        local op = util.reduce_one("stack", dim=-3, multiplicity=1+nex, name=this_name);
        pg.intern(innodes=[decon, crossview], outnodes=[op], centernodes=[ex],
                  edges=[
                      pg.edge(crossview, ex, view_index, 0), 
                      pg.edge(decon, op, view_index, 0),
                  ] + [
                      pg.edge(ex, op, wind, wind+1)
                      for wind in wc.iota(nex)
                  ]),
        

    /// Produce the layer between deoncs+crossviews for input to inference.
    prepare(name, decons, crossviews, what=["mp2","mp3"])::
        pg.crossline([
            $.prepare_pipe(name, vi, decons, crossviews, what)
            for vi in [0,1,2]]),
    
}

