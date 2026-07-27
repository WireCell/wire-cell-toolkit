// SPNG (torch) specific file I/O config

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local fans_js = import "spng/fans.jsonnet";
local fans = fans_js({});

// See and use wct/cfg/fileio.jsonnet for more

{
    // Simple torch tensor set sink.  
    pickle_tensor_set(filename, digitize=null):: pg.pnode({
        type: 'SPNGTensorSetPickleSink',
        name: filename, 
        data: {
            filename: filename,
        }
    }, nin=1, nout=0),


    // Return a "make_wrapper()" function used with wrap_pnode().
    // This assumes nodes emit individual tensors.
    pickle_wrapper(filename_pattern, is_set = function(itype) false)::
        function(label, inode, pnode)
            local filename = filename_pattern % {iname: inode.name, itype: inode.type };
            local name = label + '_' + inode.type + '_' + inode.name;
            local ttype = if is_set(inode.type)
                          then 'TensorSet'
                          else 'Tensor';
            local tap = fans.fanout(name, type=ttype);
            local pack =  pg.pnode({ // need an enclosing set for the sink
                type: 'SPNGTensorPacker',
                name: name,
                data: {
                    multiplicity: 1,
                }
            }, nin=1, nout=1);
            local sink = $.pickle_tensor_set(filename);
            local pipe =
                if is_set(inode.type)
                then pg.pipeline([pnode, tap, sink])
                else pg.pipeline([pnode, tap, pack, sink]);
            pg.intern(innodes=[pipe], oports=[tap.oports[1]]),

}

