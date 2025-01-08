local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;


function(infile, outfile, detector="uboone", variant="nominal", tags=[], digitize="no", gzip="5", chunk="256")
    local params = high.params(detector, variant);
    local mid = high.api(detector, params);
    local anode = mid.anodes()[0];
    local source = high.fio.frame_file_source(infile, tags=tags);
    local tap = pg.pnode({
        type: "HDF5FrameTap",
        name: outfile,
        data: {
            filename: outfile,
            //anode: wc.tn(anode),
            trace_tags: tags,
            digitize: digitize == "yes",
            gzip: std.parseInt(gzip),
            // lazy: support only square chunk here, but could pass rectangular array [nrows,ncols]
            chunk: std.parseInt(chunk), 
        },
    }, nin=1, nout=1, /*uses=[anode]*/);
    local sink = pg.pnode({
        type: "DumpFrames",
        name: ""
    }, nin=1, nout=0);
    local graph = pg.pipeline([source, tap, sink]);
    high.main(graph, "Pgrapher", ["WireCellHio"])

          
