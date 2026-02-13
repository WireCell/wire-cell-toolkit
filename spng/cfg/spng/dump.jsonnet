/// This provides m
/// dumped to file.
///
/// It accepts a pgraph object and returns one modified for the "hack" that
/// should then be used to construct the graph.  The selection is defined by the
/// pnode_wrappers object which maps a name a "wrapper" object.  This object
/// provides attributes "rule" that gives a rule for when the wrapper applies
/// and "wrap" which is a function that wraps a pnode into a subgraph to dump a
/// file - or really any subgraph.
///
/// The "hacked" pgraph module that is returned provides various dump_*
/// functions that are useful to call inside the wrap function.


/// information and a function

local real_pg = import 'pgraph.jsonnet';



{
    ///////////for debugging
    // Wrap a tensor producing pnode with a pickle dumper.
    local dump_tensors(tier, inode, pnode) =
        local wrap="tensors-%(tier)s-%(tpcid)s-%(itype)s-%(iname)s.pkl";
        local name = tier + '_' + inode.type + '_' + inode.name;
        local filename = wrap % {tier: tier, tpcid: std.toString(tpcid), itype:inode.type, iname:inode.name};
        local sink = tio.pickle_tensor_set(filename);
        local tap = frame.sink_taps(name, std.length(pnode.oports), sink);
        pg.shuntline(pnode, tap);

    // Wrap pnode in file dumper for tier if user has tier in dump
    local dump_tensors_maybe(tier, inode, pnode) =
        if std.findSubstr(tier,dump) == []
        then pnode
        else dump_tensors(tier, inode, pnode);

    local have_tier(tier, labels, inode) =
        std.findSubstr(tier, dump) != [] && std.findSubstr(labels[tier], inode.name) != []; 

    // Wrap pnode if a key of labels is in dump list and value matches inode.name
    local dump_tensors_matched(labels, inode, pnode) =
        local tiers = [tier
                       for tier in std.objectFields(labels)
                       if have_tier(tier, labels, inode)];
        if tiers == []
        then pnode
        else dump_tensors(tiers[0], inode, pnode);

    local pnode_wrappers = {
        SPNGRebinner: function(inode, pnode)
            dump_tensors_maybe("rebin", inode, pnode),
        SPNGResampler: function(inode, pnode)
            dump_tensors_matched({splat:"splat"}, inode, pnode),
        SPNGTransform: function(inode, pnode)
            dump_tensors_matched({dscale:"dnnroi_scale", sscale:"splat_scale"}, inode, pnode),
        SPNGKernelConvolve: function(inode, pnode)
            dump_tensors_matched({dnnroi:"dnnroi",wiener:"wiener",decon:"_group"}, inode, pnode),
        SPNGThreshold: function(inode, pnode)
            dump_tensors_matched({wthresh:"wthresh"}, inode, pnode),
        SPNGReduce: function(inode, pnode)
            dump_tensors_matched({stack:"stack"}, inode, pnode),
    };
    ///////////end for debugging
    local wpg = hacks.wrap_pnode(pnode_wrappers);
}    
