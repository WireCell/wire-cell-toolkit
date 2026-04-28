# usage: wire-cell -l stdout wct-sim-check-with-img-clustering.jsonnet
# Complete simulation + signal processing + imaging + clustering + BlobDepoFill pipeline

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);

// Track definition
local track= {
  head: wc.point(10, 467.694, 35.7493, wc.cm),
  tail: wc.point(10, 328.196, 446.613, wc.cm),
};

local tracklist = [
  {
    time: 0 * wc.us,
    charge: -500, // negative means # electrons per step (see below configuration) 
    ray: track,
  },
];

local depos = sim.tracks(tracklist, step=0.1 * wc.mm); // MIP <=> 5000e/mm

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);
local anode_idents = [anode.data.ident for anode in tools.anodes];

local drifter = sim.drifter;
local bagger = sim.make_bagger();

// signal plus noise pipelines
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/pdhd/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = {
    sparse: false,
    use_roi_debug_mode: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2]
};

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

// IMAGING PIPELINE
local img = import 'pgrapher/experiment/pdhd/img_simple.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a,'single') for a in tools.anodes];

// =============================================================================
// BlobDepoFill Configuration
// =============================================================================

// Create BlobDepoFill nodes for each anode
local blob_depo_fill_nodes = [
  g.pnode({
    type: 'BlobDepoFill',
    name: 'blobdepofill%d' % n,
    data: {
      // Drift speed - must match simulation parameters
      speed: params.lar.drift_speed,
      
      // Time offset to account for processing delays 
      // This should match the start_time used in simulation
      time_offset: params.sim.ductor.start_time,
      
      // Number of sigma for Gaussian truncation (default: 3.0)
      nsigma: 3.0,
      
      // Primary plane index (default: 2 for collection plane)
      pindex: 2,
    }
  }, nin=2, nout=1) for n in anode_iota];

// Create cluster fanout nodes to split imaging output for both regular output and BlobDepoFill
local cluster_fanouts = [
  g.pnode({
    type: 'ClusterFanout', 
    name: 'cluster_fanout%d' % n,
    data: {
      multiplicity: 2,
      tags: [],
    }
  }, nin=1, nout=2) for n in anode_iota];

// Create cluster sinks to save the truth-filled clusters
local cluster_sinks = [
  g.pnode({
    type: "ClusterFileSink",
    name: "truthfilled_clusters%d" % n,
    data: {
      outname: "truthfilled_clusters_apa%d.json" % tools.anodes[n].data.ident,
      format: "json"
    }
  }, nin=1, nout=0) for n in anode_iota];

// =============================================================================
// Modified Pipeline Construction
// =============================================================================

// Standard reconstruction pipeline with cluster fanout for BlobDepoFill
local standard_pipes = [
  g.pipeline([
    sn_pipes[n],
    nf_pipes[n],
    sp_pipes[n],
    img_pipes[n],
    cluster_fanouts[n]  // Split clusters: one for output, one for BlobDepoFill
  ], name='standard_pipe%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)
];

// BlobDepoFill pipeline that takes clusters + depos and outputs truth-filled clusters
local blobfill_pipes = [
  g.pipeline([
    blob_depo_fill_nodes[n],
    cluster_sinks[n]  // Save the truth-filled clusters
  ], name='blobfill_pipe%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)
];

// =============================================================================
// Graph Assembly with BlobDepoFill Integration
// =============================================================================

local magoutput = 'protodunehd-sim-check-full-blobfill.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local magnifyio = magnify(tools, magoutput);

// Create a more complex graph that handles the depo splitting and routing
local main_graph = g.intern(
  innodes=[],
  centernodes=standard_pipes + blobfill_pipes,
  outnodes=[],
  edges=[
    // Connect cluster fanout first output to regular processing (if needed)
    // Connect cluster fanout second output to BlobDepoFill first input (clusters)
    g.edge(standard_pipes[n], blobfill_pipes[n], 1, 0) for n in std.range(0, std.length(tools.anodes) - 1)
  ]
);

// We need to create a custom fanout structure to handle both the standard pipeline
// and routing depos to BlobDepoFill
local depo_splitter = g.pnode({
  type: 'DepoSetFanout',
  name: 'main_depo_fanout',
  data: {
    multiplicity: 2,  // One copy for standard pipeline, one copy for BlobDepoFill
    tags: []
  }
}, nin=1, nout=2);

// Route the first depo copy through the standard multifanout for signal processing
local standard_multifanout = f.multifanout('DepoSetFanout', standard_pipes, [1,nanodes], [nanodes,1], 'sn_mag_nf_img');

// Route the second depo copy to BlobDepoFill nodes (this needs custom routing)
local depo_to_blobfill = g.pnode({
  type: 'DepoSetFanout',
  name: 'depo_to_blobfill_fanout', 
  data: {
    multiplicity: nanodes,
    tags: []
  }
}, nin=1, nout=nanodes);

// Final integrated graph with proper depo routing
local integrated_graph = g.intern(
  innodes=[depo_splitter],
  centernodes=[standard_multifanout, depo_to_blobfill] + blobfill_pipes,
  outnodes=[],
  edges=[
    // Route first depo copy to standard processing
    g.edge(depo_splitter, standard_multifanout, 0, 0),
    
    // Route second depo copy to BlobDepoFill fanout
    g.edge(depo_splitter, depo_to_blobfill, 1, 0),
  ] + [
    // Route individual depo copies to each BlobDepoFill second input
    g.edge(depo_to_blobfill, blobfill_pipes[n], n, 1) for n in std.range(0, std.length(tools.anodes) - 1)
  ]
);

local sink = sim.frame_sink;

// Build final graph: depos -> drift -> bagger -> integrated processing
local graph = g.pipeline([depos, drifter, bagger, integrated_graph]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        // NOTE: Added WireCellImg plugin which is required for BlobDepoFill
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", 
                  "WireCellRoot", "WireCellTbb", "WireCellImg"],
        apps: ["Pgrapher"]
    }
};

// Finally, the configuration sequence which is emitted.
[cmdline] + g.uses(graph) + [app]