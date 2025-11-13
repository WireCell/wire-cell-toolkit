// Deal with things at frame level

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fans = import "fans.jsonnet";

/// Make a configuration for a subgraph that converts from IFrame to TDM tensor set.
{


    /// Make a FrameToTDM config object for a tpc with a single set of grouped
    /// views.  See the class documentation for details.  If more than one of
    /// these per tpc, use extra_name to distinguish.
    to_tdm(tpc, tag="", extra_name="", control={verbosity:0}):: pg.pnode({
        type: 'SPNGFrameToTdm',
        name: tpc.name + extra_name,
        data: {
            anode: wc.tn(tpc.anode),
            verbosity: control.verbosity,
            rules: [{
                tag: tag,
                groups: [{
                    wpids: g.wpids(tpc.ident),
                }, for g in tpc.view_groups]
            }],                 // just one rule
        }
    }, nin=1, nout=1, uses=[tpc.anode]),

    /// Map a tensors in a set by an iteration of a datapath pattern over a set
    /// of a count to a fan out of individual tensors.
    tensorset_unpacker(tpc,
                       extra_name="", // appended to tpc.name
                       // How to locate each tensor, %d is interpolated as an index
                       datapath_pattern="/frames/\\d+/tags/null/rules/0/groups/%d/traces")::
        local ngroups = std.length(tpc.view_groups);
        pg.pnode({
            type: 'SPNGTorchSetUnpacker',
            name: tpc.name + extra_name,
            data: {
                selections: [
                    {datapath:datapath_pattern % group}
                    for group in wc.iota(ngroups) ],
            },
        }, nin=1, nout=ngroups),

    /// Currently the pattern is to make a frame tensor set, pass that to the
    /// end of the graph and also split out the individual group tensors.  This
    /// returns that IFrame -> set + n tensor subgraph.
    set_plus_groups(tpc)::
        local tdm = $.to_tdm(tpc);
        local fanout = fans.tensor_fanout(tdm.name+"frame");
        local up = $.tensorset_unpacker(tpc); // 1->4
        // 1->1+4
        pg.intern(innodes=[tdm],
                  outnodes=[fanout, up],
                  // We add open fanout port zero
                  oports=[fanout.oports[0]]+up.oports,
                  edges=[pg.edge(tdm, fanout),
                         pg.edge(fanout, up, 1, 0)]),
        

}


