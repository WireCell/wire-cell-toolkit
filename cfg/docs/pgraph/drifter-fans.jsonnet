// drifter-fans.jsonnet
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
local graph = pg.intern(
    innodes=[source],
    centernodes=[fanin, fanout] + drifters,
    outnodes=[sink],
    edges = [
        pg.edge(source, fanout, 0, 0),
        pg.edge(fanout, drifters[0], 0, 0),
        pg.edge(fanout, drifters[1], 1, 0),
        pg.edge(drifters[0], fanin, 0, 0),
        pg.edge(drifters[1], fanin, 0, 1),
        pg.edge(fanin, sink),                            
    ]);
pg.main(graph)
