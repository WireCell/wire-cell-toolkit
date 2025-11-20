local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

{
    fanin(name, multiplicity=2, type='Tensor'):: pg.pnode({
        type: "SPNGFanin"+type+'s',
        name: name,
        data: { multiplicity: multiplicity },
    }, nin=multiplicity, nout=1),

    fanout(name, multiplicity=2, type='Tensor'):: pg.pnode({
        type: "SPNGFanout"+type+'s',
        name: name,
        data: { multiplicity: multiplicity },
    }, nin=1, nout=multiplicity),

    /// A fanout + shuntlines.  Fan out N oports of upstream node to N ports of
    /// each in list of downstream nodes.  
    fanout_shuntline(upstream, downstreams, extra_name="")::
        local Nu = std.length(upstream.oports);
        local Nd = std.length(downstreams);
        local fans = [pg.pnode({
            type: 'SPNGFanoutTensors',
            name: upstream.name + 'fan' + std.toString(ifan) + extra_name,
            data: {
                multiplicity: Nd,
            }
        }, nin=1, nout=Nd) for ifan in wc.iota(Nu)];
        pg.intern(innodes=[upstream], outnodes=downstreams, centernodes=fans,
                  edges=[
                      // edges from upstream to fans
                      pg.edge(upstream, fans[ifan], ifan, 0)
                      for ifan in wc.iota(Nu)
                  ] + [
                      pg.edge(fans[ifan], downstreams[id], id, ifan)
                      for ifan in wc.iota(Nu)
                      for id in wc.iota(Nd)
                  ]),

}
