// wrap.jsonnet
// Provides functions to enable graph rewriting by wrapping selected nodes

local real_pg = import "pgraph.jsonnet";

{
    replace_pnode(labels, finder, make_wrapper, pg = real_pg)::
        real_pg + { pnode:: $.wrap_pnode(labels, finder, make_wrapper) },

    // Create a wrapped pnode function that selectively wraps nodes based on finder rules
    //
    // Args:
    //   labels: array of strings - abstract labels indicating which nodes to wrap
    //   finder: object mapping inode.type to object mapping label to name substring
    //           e.g., {"MyType": {"label1": "substring", "label2": "other"}}
    //   make_wrapper: function(label, inode, pnode) - called to wrap matching nodes
    //
    // Returns:
    //   A function with same signature as pg.pnode() that wraps matching nodes
    wrap_pnode(labels, finder, make_wrapper)::
        function(inode, nin=0, nout=0, uses=[], name=null)
            // First create the original pnode
            local original_pnode = real_pg.pnode(inode, nin, nout, uses, name);

            // Check if this inode type is in the finder
            local inode_type = inode.type;
            local has_type = std.objectHas(finder, inode_type);

            // If type not found, return original pnode
            if !has_type then
                original_pnode
            else
                // Type found, check each label for a match
                local type_finder = finder[inode_type];
                local inode_name = if name != null then name else
                                   if std.objectHas(inode, 'name') then inode.name else "";

                // Find first matching label
                local matching_label = std.foldl(
                    function(found, label)
                        // If already found a match, keep it
                        if found != null then found
                        // Check if this label exists in type_finder
                        else if std.objectHas(type_finder, label) then
                            local substring = type_finder[label];
                            // Check if substring matches the inode name
                            local matches = std.length(std.findSubstr(substring, inode_name)) > 0;
                            if matches then label else null
                        else null,
                    labels,
                    null  // initial value
                );

                // If a matching label was found, wrap the node
                if matching_label != null then
                    make_wrapper(matching_label, inode, original_pnode)
                else
                    original_pnode,
}
