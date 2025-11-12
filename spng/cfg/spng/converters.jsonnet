local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";


/// Make a configuration for a subgraph that converts from IFrame to TDM tensor set.
{

    // Construct WPIDs based on faces and connections for one view.
    view_group_wpids(tpc, view_index)::
        // Start with each face separate.
        local split = [[wc.WirePlaneId(wc.wpid_index_to_layer(view_index), face, tpc.ident)]
                       for face in tpc.faces];
        if tpc.connections[view_index] > 0
        then [wc.flatten(split)]
        else split,

    // Combine across views into flat list-of-list-of-wcpid..
    tpc_group_wpids(tpc)::
        wc.flatten([$.view_group_wpids(tpc, index) for index in [0,1,2]]),

    view_group_indices(tpc, view_index)::
        local split = [[view_index] for face in tpc.faces];
        if tpc.connections[view_index] > 0
        then [view_index]
        else [view_index, view_index],

    // Get plane index for each group.
    tpc_group_planes(tpc):
        wc.flatten([$.view_group_indices(tpc, index) for index in [0,1,2]]),


    /// Make a FrameToTDM config object for a tpc with a single set of grouped
    /// views.  See the class documentation for details.  If more than one of
    /// these per tpc, use extra_name to distinguish.
    frame_to_tdm_views(tpc, tag="", extra_name="", control={verbosity:0}):: pg.pnode({
        type: 'SPNGFrameToTdm',
        name: tpc.name + extra_name,
        data: {
            anode: wc.tn(tpc.anode),
            verbosity: control.verbosity,
            rules: [{
                tag: tag,
                groups: [{
                    wpids: wpids,
                }, for wpids in $.tpc_group_wpids(tpc)],
            }],                 // just one rule
        }
    }, nin=1, nout=1, uses=[tpc.anode]),

    /// Map a tensors in a set by an iteration of a datapath pattern over a set
    /// of a count to a fan out of individual tensors.
    frame_tensorset_unpacker(tpc,
                             extra_name="", // appended to tpc.name
                             // How to locate each tensor, %d is interpolated as an index
                             datapath_pattern="/frames/\\d+/tags/null/rules/0/groups/%d/traces")::
        local ngroups = std.length($.tpc_group_planes(tpc));
        pg.pnode({
            type: 'SPNGTorchSetUnpacker',
            name: tpc.name + extra_name,
            data: {
                selections: [
                    {datapath:datapath_pattern % group}
                    for group in wc.iota(ngroups) ],
            },
        }, nin=1, nout=ngroups),
}


