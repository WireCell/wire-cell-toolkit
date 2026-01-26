// drifter-intern.jsonnet
local pg = import "pgraph.jsonnet";
local source = pg.pnode({
    type: 'DepoFileSource', name:"",
}, nin=0, nout=1);
local drifter = pg.pnode({
    type: 'DepoSetDrifter', name:"",
}, nin=1, nout=1);
local sink = pg.pnode({
    type: 'DepoFileSink', name:"",
}, nin=1, nout=0);
local graph = pg.intern(innodes=[source], centernodes=[drifter], outnodes=[sink]);
pg.main(graph)
