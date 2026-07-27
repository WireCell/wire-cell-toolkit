// drifter-fanpipe.jsonnet
local pg = import "pgraph.jsonnet";
local source = pg.pnode({
    type: 'DepoFileSource', name:"",
}, nin=0, nout=1);
local fanout = pg.pnode({
    type: 'DeposSetFanout', name:"",
}, nin=1, nout=2);
local drifters = [
    pg.pnode({
        type: 'DepoSetDrifter', name:"one",
    }, nin=1, nout=1),
    pg.pnode({
        type: 'DepoSetDrifter', name:"two",
    }, nin=1, nout=1)
];
local fanin = pg.pnode({
    type: 'DepoSetFanin', name:"",
}, nin=2, nout=1);
local sink = pg.pnode({
    type: 'DepoFileSink', name:"",
}, nin=1, nout=0);
local fanpipe = pg.fan.pipe('DepoSetFanout', drifters, 'DepoSetFanin');
local graph = pg.pipeline([source, fanpipe, sink]);
pg.main(graph)
