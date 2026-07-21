local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';

function() {
    local frame_source = g.pnode({
        type: 'SPNGTorchTensorFileSource',
        name: 'ttfsource',
        data: {
            inname: 'test.tar',
        },
    }, nin=0, nout=1),
    local frame_sink = g.pnode({
        type: 'SPNGTorchTensorFileSink',
        name: 'ttfsink',
        data: {
            outname: 'testout.npy',
            prefix: './'
        },
    }, nin=1, nout=0),
local sig_proc = g.pnode({
        type: 'SPNGSigProc',
        name: 'spngsigproc'
    }, nin=1, nout=1),

// local graph = g.pipeline([input.frame_source, filter.sig_proc]);
local graph = g.pipeline([frame_source, frame_sink]),

local app = {
    type: 'Pgrapher',
    data: {
        edges: g.edges(graph),
    },
},

local plugins = ["WireCellSpng", "WireCellPgraph"],

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: plugins,
        apps: ["Pgrapher"],
    }
},
seq: [cmdline] + g.uses(graph) + [app]
}.seq