// Deal with things at frame level

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fans = import "fans.jsonnet";

/// Make a configuration for a subgraph that converts from IFrame to TDM tensor set.
{


    /// Make a FrameToTDM config object for a tpc with a single set of grouped
    /// views.  See the class documentation for details.  If more than one of
    /// these per tpc, use extra_name to distinguish.
    to_tdm(tpc, tag="", extra_name="", control={}):: pg.pnode({
        type: 'SPNGFrameToTdm',
        name: tpc.name + extra_name,
        data: {
            anode: wc.tn(tpc.anode),
            rules: [{
                tag: tag,
                groups: [{
                    wpids: g.wpids(tpc.ident),
                }, for g in tpc.view_groups]
            }],                 // just one rule
        } + wc.object_with(control, ["verbosity", "device"]),
    }, nin=1, nout=1, uses=[tpc.anode]),

    /// Convert a TDM frame tensor set back into IFrame
    from_tdm(tpc, extra_name="", control={}):: pg.pnode({
        type: 'SPNGTdmToFrame',
        name: tpc.name + extra_name,
        data: {
            frame: {datapath: "/frames/\\d+/frame"},
            // chmasks: ...
            tagged_traces: [ {
                // eg datapath of /frames/0/tags/gauss/groups/0/traces
                traces: { tag: "gauss" },
                // eg datapath of /frames/0/tags/null/rules/0/groups/0/chids
                chids: { tag: "null" },
            }]
        } + wc.object_with(control, ["verbosity", "device"]),
    }, nin=1, nout=1, uses=[]),

        

    /// Map a tensors in a set by an iteration of a datapath pattern over a set
    /// of a count to a fan out of individual tensors.
    tensorset_unpacker(tpc,
                       extra_name="", // appended to tpc.name
                       // How to locate each tensor, %d is interpolated as an index
                       datapath_pattern="/frames/\\d+/tags/null/rules/0/groups/%d/traces",
                       control={})::
        local ngroups = std.length(tpc.view_groups);
        pg.pnode({
            type: 'SPNGTorchSetUnpacker',
            name: tpc.name + extra_name,
            data: {
                selections: [{datapath:datapath_pattern % group} for group in wc.iota(ngroups) ],
            } + wc.object_with(control, ["verbosity"]),
        }, nin=1, nout=ngroups),

    /// Repack tensors of a view.  FIXME: address asymmetry between unpacking
    /// groups and packing views.  For now, add "_view_" name to the function.
    tensorset_view_repacker(name, multiplicity=3, control={}):: pg.pnode({
        type: 'SPNGTensorPacker',
        name: name,
        data: {
            multiplicity: multiplicity,
        } + wc.object_with(control, ["verbosity"]),
    }, nin=multiplicity, nout=1),


    /// Currently the pattern is to make a frame tensor set, pass that to the
    /// end of the graph and also split out the individual group tensors.  This
    /// returns that IFrame -> set + n tensor subgraph.
    set_plus_groups(tpc, control={})::
        local tdm = $.to_tdm(tpc);
        local fanout = fans.fanout(tdm.name+"frame", control=control);
        local up = $.tensorset_unpacker(tpc, control=control); // 1->4
        // 1->1+4
        pg.intern(innodes=[tdm],
                  outnodes=[fanout, up],
                  // We add open fanout port zero
                  oports=[fanout.oports[0]]+up.oports,
                  edges=[pg.edge(tdm, fanout),
                         pg.edge(fanout, up, 1, 0)]),
        
    /// A "bypass" is a 2->2 node bridges a 1->2 fanout to a 2->1 fanin.  One
    /// finger of the fans are connected (the bypass) and the other fingers are
    /// left open
    ///
    /// The input on iport=0 goes to the fanout and is immediately available on
    /// oport=0.  That is, port zeros act like a pass-through.  The fanout_index
    /// of the fanout is connected to fanin_index of the fanin.  Input on
    /// iport=1 also goes to the fanin.  And final output of the fanin is
    /// available on oport=1.  Ports zero do not change the data while ports one
    /// expect or produce further processed data.
    ///
    /// This sounds illogical and like a cycle is set up but that is a mirage
    /// due to the bypass being a subgraph:
    ///
    ///                        the bypass
    /// (input0) ---> [fanout] ----->---- [fanin] --> (output1)
    ///                 |                   ^
    ///                 |                   |
    /// (input1)--------|-------------------+
    ///                 |
    ///                 +---------------------------> (output0 == source(input0))
    ///
    /// The bypass is useful when one needs to combine processed input with
    /// original input.  Eg:
    ///
    /// (input0) -> [0:bypass:0] -> [processing] -> [1:bypass:1] -> (output1)
    ///
    /// The fan in/out indices are determine the bypass connection if order
    /// matters to the fans.
    bypass(name, fanout_index=0, fanin_index=0, type='TensorSet', control={})::
        local expose_fanin_index = (fanin_index + 2)%2;
        local expose_fanout_index = (fanout_index + 2)%2;
        local fanout = pg.pnode({
            type: "SPNGFanout" + type + "s",
            name: name,
            data: {
                multiplicity: 2, // fixme: could be option
            } + wc.object_with(control, ["verbosity"]), 
        }, nin=1, nout=2);
        local fanin = pg.pnode({
            type: "SPNGFanin" + type + "s",
            name: name,
            data: {
                multiplicity: 2, // fixme: could be option
            } + wc.object_with(control, ["verbosity"]),
        }, nin=2, nout=1);
        pg.intern(innodes=[fanout], outnodes=[fanin],
                  iports=[fanout.iports[0], fanin.iports[expose_fanin_index]],
                  oports=[fanout.oports[expose_fanout_index], fanin.oports[0]],
                  edges=[pg.edge(fanout, fanin, fanout_index, fanin_index)]),

    /// Connect bypass between two sources and targets. 
    connect_bypass(bypass, sources, targets)::
        pg.intern(innodes=sources, outnodes=targets, centernodes=[bypass],
                  edges=[
                      pg.edge(sources[0], bypass, 0, 0),
                      pg.edge(sources[1], bypass, 0, 1),
                      pg.edge(bypass, targets[0], 0, 0),
                      pg.edge(bypass, targets[1], 1, 0)]),
    
    /// Return a subgraph that shunts across inputs of given multiplicity to
    /// outputs of same and that taps each into the given tensor set sink.  The
    /// sink is interned.
    ///
    sink_taps(name, multiplicity, sink, control={}, digitize=null)::
        local fanouts = [fans.fanout(name + std.toString(index), control=control)
                         for index in wc.iota(multiplicity)];
        local packer = $.tensorset_view_repacker(name, multiplicity, control=control);
        pg.intern(innodes=fanouts, centernodes=[packer, sink],
                  oports=[fo.oports[0] for fo in fanouts], // fanout port 0 of each is out
                  edges=[
                      pg.edge(fo.value, packer, 1, fo.index) // fanout port 1 to packer
                      for fo in wc.enumerate(fanouts)
                  ] + [
                      pg.edge(packer, sink) // packer to sink
                  ]),
        
}



