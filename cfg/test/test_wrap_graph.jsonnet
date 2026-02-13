local real_pg = import "pgraph.jsonnet";
local wrap = import "wrap.jsonnet";

local labels = ["dump_source", "dump_body"];
local finder = {
    Source: { dump_source: "source" },
    Body: { dump_body: "body" }
};

local make_wrapper(label, inode, pnode) =
    local tap = real_pg.pnode({type:'Fanout',name:'fanout-'+label,data:{}}, nin=1, nout=2);
    local sink = real_pg.pnode({type:'Sink',name:'sink-'+label,data:{}}, nin=1, nout=0);
    real_pg.intern(innodes=[tap],
                   oports=[tap.oports[0]],
                   edges=[real_pg.edge(tap, sink, 1, 0)]);

local pg = real_pg + { pnode:: wrap.wrap_pnode(labels, finder, make_wrapper) };


local source = pg.pnode({type:'Source',name:"source", data:{}}, nin=0, nout=1);
local body = pg.pnode({type:'Body',name:"body", data:{}}, nin=1, nout=1);
local sink = pg.pnode({type:'Sink',name:"sink", data:{}}, nin=1, nout=0);

local graph = pg.pipeline([source, body, sink]);
pg.main(graph, 'Pgrapher')
