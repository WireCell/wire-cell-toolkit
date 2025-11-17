local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

{
    pickle_tensor_set(filename): pg.pnode({
        type: 'TensorSetPickleSink',
        name: filename, 
        data: {
            filename: filename,
        }
    }, nin=1, nout=0),

}
