local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
local util = import "util.jsonnet";

function(params, anode, tags=[], name=null)
    pg.pnode({
        type: 'Reframer',
        name: if std.type(name) == "null" then util.idents(anode) else name,
        data: {
            anode: wc.tn(anode),
            tags: tags,
            fill: 0.0,
            toffset: 0,
            tbin: params.ductor.tbin,
            nticks: params.ductor.binning.nticks - self.tbin,
        },
    }, nin=1, nout=1, uses=[anode])

