// SPNG (torch) specific file I/O config

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

// See and use wct/cfg/fileio.jsonnet for more

{
    // Simple torch tensor set sink.  
    pickle_tensor_set(filename, digitize=null): pg.pnode({
        type: 'SPNGTensorSetPickleSink',
        name: filename, 
        data: {
            filename: filename,
        }
    }, nin=1, nout=0),

}
