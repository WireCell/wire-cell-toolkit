// WARNING: node config objects here are all made up

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// Splicing


// Build nominal / original graph
local source = pg.pnode({ type:'Source', name:"origsource", data:{} }, nin=0, nout=1);
local fanout = pg.pnode({ type:'Fanout', name:"origfanout", data:{} }, nin=1, nout=2);
local pipes = [
    pg.pnode({ type:'Pipe', name:"orig0", data:{} }, nin=1, nout=1),
    pg.pnode({ type:'Pipe', name:"orig1", data:{} }, nin=1, nout=1),
];
local fanin = pg.pnode({ type:'Fanin', name:"origfanin", data:{} }, nin=2, nout=1);
local sink = pg.pnode({ type:'Sink', name:"origsink", data:{} }, nin=1, nout=0);

// contrive some abstract subgraphs
local origfe = pg.intern(centernodes=[source],
                         outnodes=[fanout],
                         edges=[
                             pg.edge(source,fanout),
                         ], name="origfe");
local origbe = pg.intern(centernodes=[sink],
                         innodes=[fanin],
                         edges=[
                             pg.edge(fanin, sink),
                         ], name="origbe");
local orig = pg.intern(innodes=[origfe], outnodes=[origbe], centernodes=pipes,
                       edges=[
                           pg.edge(fanout,pipes[0],0,0),
                           pg.edge(fanout,pipes[1],1,0),
                           pg.edge(pipes[0],fanin,0,0),
                           pg.edge(pipes[1],fanin,0,1),
                       ], name="orig");

local orig_mono = pg.intern(centernodes=[source,fanout,fanin]+pipes,
                       edges=[
                           pg.edge(source,fanout),
                           pg.edge(fanout,pipes[0],0,0),
                           pg.edge(fanout,pipes[1],1,0),
                           pg.edge(pipes[0],fanin,0,0),
                           pg.edge(pipes[1],fanin,0,1),
                           pg.edge(fanin,sink)
                       ], name="origmono");

// Build an incomplete subgraph ending in a second sink to be spliced.
local ofanin = pg.pnode({ type: 'Fanin', name:"outfanin", data:{} }, nin=2, nout=1);
local osink = pg.pnode({ type: 'Sink', name:"outsink", data:{} }, nin=1, nout=0);
local out = pg.intern(innodes=[ofanin,], centernodes = [osink],
                      edges=[ pg.edge(ofanin, osink) ], name="outsgr");

local edge_selector(e) = std.startsWith(e.tail.node, "Pipe:");
local fanout_factory(n,e) = { type:'Fanout', name:"splice%d"%n };

local spliced = pg.splice(orig, out, edge_selector, fanout_factory);

local main(graph, engine='TbbFlow') = pg.uses(graph) + [{type:engine, data: { edges: pg.edges(graph) }}];

local maybe_output = {
    orig:main(orig),
    out:main(out),
    spliced:main(spliced),
};

function(what) maybe_output[what]
