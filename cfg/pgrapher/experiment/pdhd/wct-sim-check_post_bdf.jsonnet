local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';


local tools = tools_maker(params);
local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);

local sim = sim_maker(params, tools);
local bagger = sim.make_bagger();

# Load depositions (same for all anodes - contains all depos)
local depo_loader = g.pnode({
  type: "NumpyDepoSetLoader",
  name: "depo_loader",
  data: { filename: "depositions.npz" }
}, nin=0, nout=1);

# Create fanout to distribute depos to all anodes
local depo_fanout = g.pnode({
  type: 'DepoSetFanout',
  name: 'depo_fanout_all_anodes',
  data: { multiplicity: nanodes }
}, nin=1, nout=nanodes);

# Load clusters for each anode
local cluster_loaders = [
  g.pnode({
    type: "ClusterFileSource",
    name: "cluster_loader_apa%d" % n,
    data: {
      inname: "clusters-apa-apa%d.tar.gz" % tools.anodes[n].data.ident,
      anodes: [wc.tn(tools.anodes[n])]
    }
  }, nin=0, nout=1, uses=[tools.anodes[n]]) for n in anode_iota
];


# Create BlobDepoFill nodes for each anode
local blob_depo_fills = [
  g.pnode({
    type: 'BlobDepoFill',
    name: 'blobdepofill_apa%d' % n,
    data: {
      speed: params.lar.drift_speed,
      time_offset: params.sim.ductor.start_time,
      nsigma: 3.0,
      pindex: 2,
    }
  }, nin=2, nout=1) for n in anode_iota
];

# Create output sinks for each anode
local truth_sinks = [
  g.pnode({
    type: "ClusterFileSink",
    name: "truth_sink_apa%d" % n,
    data: {
      outname: "truthfilled_clusters_apa%d.json" % tools.anodes[n].data.ident,
      format: "json"
    }
  }, nin=1, nout=0) for n in anode_iota
];

# Create processing pipeline for each anode
local per_anode_pipelines = [
  g.pipeline([
    blob_depo_fills[n],
    truth_sinks[n]
  ], name='blobfill_pipe_apa%d' % n) for n in anode_iota
];

# Create the complete graph
local complete_graph = g.intern(
  innodes=[depo_loader] + cluster_loaders,
  centernodes=[depo_fanout] + per_anode_pipelines,
  outnodes=[],
  edges=[
    # Connect depo loader to fanout
    // g.edge(depo_loader, bagger, 0, 0),          # IDepo → DepoBagger
    // g.edge(bagger, depo_fanout, 0, 0),          # IDepoSet → DepoSetFanout
    g.edge(depo_loader, depo_fanout, 0, 0),
  ] + [
    # Connect each cluster loader to corresponding BlobDepoFill (input 0)
    g.edge(cluster_loaders[n], per_anode_pipelines[n], 0, 0) for n in anode_iota
  ] + [
    # Connect depo fanout to each BlobDepoFill (input 1)
    g.edge(depo_fanout, per_anode_pipelines[n], n, 1) for n in anode_iota
  ]
);

local app = {
  type: 'Pgrapher',
  data: { edges: g.edges(complete_graph) }
};

local cmdline = {
  type: "wire-cell",
  data: {
    plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellImg"],
    apps: ["Pgrapher"]
  }
};

[cmdline] + g.uses(complete_graph) + [app]