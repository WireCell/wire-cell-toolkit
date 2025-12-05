local pg = import "pgraph.jsonnet";

{
    
    /// Apply a reduction operation
    reduce_one(op, dim=0, multiplicity=2, name="")::
        pg.pnode({
            type:'SPNGReduce',
            name: name,
            data: {
                operation: op,
                dim: dim,
                multiplicity: multiplicity,
            },
        }, nin=multiplicity, nout=1),

    /// A stackline of the same reduces
    reduce(op, dim=0, multiplicity=2, name="")::
        pg.stackline([
            $.reduce_one(op, dim, multiplicity, name="v"+std.toString(view_index)+name)
            for view_index in [0,1,2]]),


}
