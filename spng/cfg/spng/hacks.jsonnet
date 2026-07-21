// This file holds some "hacks" to do useful if questionable things.

local real_pg = import "pgraph.jsonnet";


{

    // This wraps pg.pnode() with a function provided by class_wrappers.
    //
    // The "pg" will be pgraph.jsonnet if not otherwise provided.
    //
    // Class wrappers maps a "type" name to a function(inode, pnode) that is
    // given the underlying inode and the pnode returned by the "real"
    // pg.pnode().  The object returned by this function is then the pnode.
    //
    // This can be used to provide a subgraph aroudn a "real" pnode that acts as
    // the pnode itself as far as neighbor nodes know.  For example, this may be
    // used to dump input/output data to/from all nodes of the target type.
    pnode_wrapper(type_wrappers, pg=real_pg) ::
        function(inode, nin=0, nout=0, uses=[], name=null) // pnode
            //local itype=std.trace('inode.type=%s'%inode.type, inode.type);
            local itype=inode.type;
            local w = std.get(type_wrappers, itype, function(inode, pnode) pnode);
            w(inode, pg.pnode(inode, nin=nin, nout=nout, uses=uses, name=name)),


    // Override pg.pnode with the wrapped version.
    wrap_pnode(type_wrappers, pg=real_pg, pnode_wrapper=$.pnode_wrapper)::
        pg + { pnode:: pnode_wrapper(type_wrappers)},
            

}
