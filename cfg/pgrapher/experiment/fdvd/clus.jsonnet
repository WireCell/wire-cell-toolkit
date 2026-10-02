// FD-VD low-energy clustering nodes (fdvd_sim doc 16).  Counterpart of fdvd_sim/cfg/wct-clus-fdvd.jsonnet
// (wcp WireCell/wcp-porting-validation; doc 04 sec 4), split into builder functions so the in-toolkit
// low-energy reconstruction (lowe-reco.jsonnet) can tap the per-CRM imaging stream for the Q-L blob table.
// Node types, names and data are those of the doc 04 job, so the compiled clustering is the same:
//   per CRM  ClusterFileSource(clusters-apa-crm<N>-<pipe>.tar.gz) -> PointTreeBuilding (BlobSampler
//   'stepped' with center_fallback), all -> PointTreeMerging -> MultiAlgBlobClustering
//   (pointed, close 1.2 cm, isolated at all-CRM scope with small_small = small_big = iso_cm, bbox prefilter)
local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local clus = import 'pgrapher/common/clus.jsonnet';

function(params, anodes) {
  local v = params.lar.drift_speed,
  // one metadata block per CRM wire-plane id; every CRM is the same (single face 0, drift to -x)
  local crm_md = {
    drift_speed: v,
    tick: 0.5 * wc.us,
    tick_drift: self.drift_speed * self.tick,
    time_offset: -8.7 * wc.us,      // doc 02 sec 6.3: reco time = arrival + 8.7 us
    trigger_offset: 0,
    nticks_live_slice: 4,
    FV_xmin: -3250.0 * wc.mm,       // cathode
    FV_xmax: 3234.3 * wc.mm,        // shield plane (the anode slab is dead, doc 02 sec 6.1)
    FV_xmin_margin: 2 * wc.cm,
    FV_xmax_margin: 2 * wc.cm,
  },
  dv:: {
    type: 'DetectorVolumes',
    name: 'dv-fdvd',
    data: {
      anodes: [wc.tn(a) for a in anodes],
      metadata: {
        overall: {
          FV_xmin: -3250.0 * wc.mm, FV_xmax: 3234.3 * wc.mm,
          FV_ymin: -6751.0 * wc.mm, FV_ymax: 6751.0 * wc.mm,
          FV_zmin: 8.5 * wc.mm, FV_zmax: 21003.3 * wc.mm,
          FV_xmin_margin: 2 * wc.cm, FV_xmax_margin: 2 * wc.cm,
          FV_ymin_margin: 2.5 * wc.cm, FV_ymax_margin: 2.5 * wc.cm,
          FV_zmin_margin: 3 * wc.cm, FV_zmax_margin: 3 * wc.cm,
          vertical_dir: [0, 1, 0],
          beam_dir: [0, 0, 1],
        },
      } + { ['a%df0pA' % a.data.ident]: crm_md for a in anodes },
    },
    uses: anodes,
  },
  local dv = self.dv,
  local pcts = {
    type: 'PCTransformSet',
    name: dv.name,
    data: { detector_volumes: wc.tn(dv) },
    uses: [dv],
  },

  // imaging archive of one CRM
  source(a, indir, pipe):: g.pnode({
    type: 'ClusterFileSource',
    name: 'crm%d' % a.data.ident,
    // restore_corners: carry the imaging-time blob corners onto the loaded blobs for FdvdBlobTable (the light
    // prediction's blob centre).  It feeds only the dead-blob "corner" PC otherwise; FD-VD has no dead blobs,
    // so the clustering is unchanged (doc 16 gate B).
    data: { inname: '%s/clusters-apa-%s-%s.tar.gz' % [indir, a.name, pipe], anodes: [wc.tn(a)],
            restore_corners: true },
  }, nin=0, nout=1, uses=[a]),

  point_tree(a):: g.pnode({
    type: 'PointTreeBuilding',
    name: 'crm%d' % a.data.ident,
    data: {
      samplers: { '3d': wc.tn($.sampler(a)) },
      multiplicity: 1,
      tags: ['live'],
      anode: wc.tn(a),
      face: 0,
      detector_volumes: wc.tn(dv),
    },
  }, nin=1, nout=1, uses=[$.sampler(a), dv]),

  sampler(a):: {
    type: 'BlobSampler',
    name: 'live-crm%d' % a.data.ident,
    data: {
      drift_speed: v,
      time_offset: -8.7 * wc.us,
      strategy: [{ name: 'stepped', center_fallback: true }],
      extra: ['.*wire_index', '.*charge_val', '.*charge_unc', 'wpid'],
    },
  },

  merge(n):: g.pnode({
    type: 'PointTreeMerging',
    name: 'fdvd-all',
    // tolerate_missing: the FD-VD point trees carry no 'dead' subtree (no dead-channel archive)
    data: { multiplicity: n, inpath: 'pointtrees/%d', outpath: 'pointtrees/%d', tolerate_missing: true },
  }, nin=n, nout=1),

  // bee_zip '' = no Bee output
  mabc(iso=true, iso_cm=30, bbox=true, run=0, subrun=0, event=0, bee_zip=''):
    local cm = clus.clustering_methods(prefix='fdvd', detector_volumes=dv, pc_transforms=pcts,
                                       coords=['x', 'y', 'z']);
    local pipeline = [cm.pointed(), cm.close(length_cut=1.2 * wc.cm)] +
      (if iso then [cm.isolated(small_small_dis_cut=if iso_cm == null then null else iso_cm * wc.cm,
                                small_big_dis_cut=if iso_cm == null then null else iso_cm * wc.cm,
                                bbox_prefilter=bbox)] else []);
    g.pnode({
      type: 'MultiAlgBlobClustering',
      name: 'fdvd-all',
      data: {
        inpath: 'pointtrees/%d',
        outpath: 'pointtrees/%d',
        perf: true,
        initial_index: 0,
        use_config_rse: true,
        runNo: run,
        subRunNo: subrun,
        eventNo: event,
        anodes: [wc.tn(a) for a in anodes],
        detector_volumes: wc.tn(dv),
        pipeline: wc.tns(pipeline),
      } + (if bee_zip != '' then {
        bee_zip: bee_zip,
        bee_detector: 'dunefdvd-1x8x14',
        bee_points_sets: [{
          name: 'clustering',
          detector: 'dunefdvd-1x8x14',
          algorithm: 'clustering',
          pcname: '3d',
          coords: ['x', 'y', 'z'],
          individual: false,
        }],
      } else {}),
    }, nin=1, nout=1, uses=anodes + [dv, pcts] + pipeline),
}
