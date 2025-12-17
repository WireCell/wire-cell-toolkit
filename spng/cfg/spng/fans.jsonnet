local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

function(control={})
{
    fanin(name, multiplicity=2, type='Tensor'):: pg.pnode({
        type: "SPNGFanin"+type+'s',
        name: name,
        data: { multiplicity: multiplicity } + control
    }, nin=multiplicity, nout=1),

    fanout(name, multiplicity=2, type='Tensor'):: pg.pnode({
        type: "SPNGFanout"+type+'s',
        name: name,
        data: { multiplicity: multiplicity } + control
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
            } + control
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

    /// Forward N ports M ways.
    ///
    /// It returns a list like:
    ///
    /// [N-sink, [N-source]*M]
    ///
    ///
    /// Each of N iports of the N-sink should be connected to upstream.
    ///
    /// Each of the N iports to each of the M N-sources should be connected to
    /// downstream.
    ///
    /// N is typically number of detector views in a TPC or number of TPCs in a
    /// detector.  M is whatever fanout number you want from each.
    fanout_cross(name, N, M, type='Tensor')::
        local fans = [pg.pnode({
            type: "SPNGFanout"+type+"s",
            name: name + "f" + std.toString(num),
            data: { multiplicity: M } + control
        }, nin=1, nout=M) for num in wc.iota(N)];
        local sink = pg.intern(innodes=fans); // sets their iports
        local sources = [
            pg.intern(centernodes=fans, // fixme, even give this?
                      oports=[
                          fans[fnum].oports[mnum]
                          for fnum in wc.iota(N)
                      ])
            for mnum in wc.iota(M)
        ];
        [sink, sources],

    /// An N-by-M fanout like fanout_cross() but caller provides function
    /// ifunc to generate an inode that is the fanout for one num in N.  It is
    /// called:
    ///
    ///   ifunc(num, M) -> inode
    ///
    /// with num in half-open range [0,N) and M being the multiplicity that
    /// may be required to configure the returned inode.  The pnode is made with
    /// the result.
    ///
    /// Return [sink, [source, source, ...]]
    fanout_cross_gen(N, M, ifunc)::
        local the_fans = [pg.pnode(ifunc(num, M), nin=1, nout=M) for num in wc.iota(N)];
        local sink = pg.intern(innodes=the_fans); // sets their iports
        local sources = [
            pg.intern(centernodes=the_fans, // fixme, even give this?
                      oports=[
                          the_fans[fnum].oports[mnum]
                          for fnum in wc.iota(N)
                      ])
            for mnum in wc.iota(M)
        ];
        [sink, sources],



    /// Fanout N source ports according to targets.
    ///
    /// Targets is a list of a list of string.  Each list-of-string gives the
    /// target group names for a given source port.  One Pnode for each group
    /// name is returned which has oports connected to one of corresponding
    /// source port.
    ///
    /// An object is returned with "sink" providing N iports and "targets"
    /// providing an object with keys as in targets and values a Pnode.
    ///
    /// For example:
    /// fanout_select("", 3, [["w","g","d"], ["w","g","d"], ["w","g"]])
    ///
    /// returns:
    /// {
    ///   sink: <a 3 iport pnode>,
    ///   targets: {
    ///      w: <a 3 oport pnode>
    ///      g: <a 3 oport pnode>
    ///      d: <a 2 oport pnode>
    ///   }
    /// }
    // The "d" target's 2 ports connect to source port 0 and 1.
    fanout_select(name, N, targets_list, type='Tensor')::
        local tnames = std.set(std.flattenArrays(targets_list));
        local the_fans = [
            local M = std.length(targets_list[num]);
            pg.pnode({
                type: "SPNGFanout"+type+"s",
                name: name + "f" + std.toString(num),
                data: { multiplicity: M} + control
            }, nin=1, nout=M) for num in wc.iota(N)];

        local fan_port(fanind, group) =
                  wc.index_of(targets_list[fanind], group);

        local group_fan_index = {
            [name]:[x.index for x in wc.enumerate(targets_list) if std.member(x.value, name)]
            for name in tnames
        };

        local tobj = {
            [group]: pg.intern(oports=[
                the_fans[fanind].oports[fan_port(fanind, group)]
                for fanind in wc.iota(N) if std.member(targets_list[fanind], group)])
            for group in tnames
        };

        {
            sink: pg.crossline(the_fans),
            targets: tobj,

        },

    /// Like fanout_select but produce the fanouts with a generator function
    /// ifunc.
    fanout_select_gen(name, N, targets_list, ifunc)::
        local tnames = std.set(std.flattenArrays(targets_list));
        local the_fans = [
            local M = std.length(targets_list[num]);
            pg.pnode(ifunc(num, M), nin=1, nout=M)
            for num in wc.iota(N)];

        local fan_port(fanind, group) =
                  wc.index_of(targets_list[fanind], group);

        local group_fan_index = {
            [name]:[x.index for x in wc.enumerate(targets_list) if std.member(x.value, name)]
            for name in tnames
        };

        local tobj = {
            [group]: pg.intern(oports=[
                the_fans[fanind].oports[fan_port(fanind, group)]
                for fanind in wc.iota(N) if std.member(targets_list[fanind], group)])
            for group in tnames
        };

        {
            sink: pg.crossline(the_fans),
            targets: tobj,

        },


    /// Return two nodes: [orig, oport] where "orig" is effectively the input
    /// node and "oport" is a node that is a replicated source of data from the
    /// output port of "node" at the index oport_index.
    forkout_oport(name, node, oport_index, type='Tensor')::
        local fout = $.fanout(name+'forkout', type=type);
        local copy = pg.intern(innodes=[node],
                               centernodes=[fout],
                               oports=[
                                   if op.index == oport_index
                                   then fout.oports[0]
                                   else op.value
                                   for op in wc.enumerate(node.oports)
                               ],
                               edges=[
                                   pg.edge(node, fout, oport_index, 0)
                               ]);
        local out = pg.intern(oports=[fout.oports[1]]);
        [copy, out],

}
