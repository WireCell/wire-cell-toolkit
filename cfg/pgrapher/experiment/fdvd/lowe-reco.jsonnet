// FD-VD low-energy (solar) reconstruction in one wire-cell job (fdvd_sim doc 16).
//
//   wcsonnet -A indir=<per-CRM imaging archives> -A opwf=<opwf tensor archive> \
//            -A library=<fdvd-photlib-vis-comb-10cm.json> -A calibration=<fdvd-lowe-qlcal-R20pe.json> \
//            -A rootfile=lowe.root [-A outfile=lowe.tar.gz] [-A drift_model=<m3 .ts> -A frames='<dir>/sp-anode%d.tar.gz' \
//            -A drift_device=gpu] [-A drift_in=<drift tensor archive>] [-A bee_zip=mabc.zip] \
//            [--tla-code event=N] cfg/pgrapher/experiment/fdvd/lowe-reco.jsonnet
//
// WIRECELL_PATH must also hold the FD-VD wires file dunevd10kt_3view_30deg_v7_refactored_1x8x14.json.bz2
// (fdvd_sim/data; not in wire-cell-data).
//
// Graph:
//   charge  per CRM  ClusterFileSource -> ClusterFanout(2) -> PointTreeBuilding -> PointTreeMerging(N)
//                                                                         -> MultiAlgBlobClustering ---> QL port 0
//                                                         \-> FdvdBlobTable(N) ----------------------> QL port 1
//   light   opwf:    TensorFileSource(opwf_) -> FdvdOpHitFinder -> OpFlashFinder(frozen doc 04 point) -> QL port 2
//           opflash: TensorFileSource(opflash_) (an OpFlashFinder archive)                            -> QL port 2
//   drift   drift_model + frames: the matcher runs FdvdDriftRegressor itself (TorchTensorSetService with the
//           TorchScript M3, fdvd_sim/stageB/export_drift_ts.py; SP gauss frames read from the per-CRM archives
//           named by frames = '<dir>/sp-anode%d.tar.gz'), or
//           drift_in: TensorFileSource(drift_) ("drift" [n, 3] = rep blob, mu, sigma)                 -> QL port 3
//   FdvdLowEQLMatching -> [FdvdLowERootWriter(rootfile), pass-through] -> TensorFileSink(lowe_; dump_mode if outfile '')
// Inputs are the per-CRM imaging archives (the "assembled" mode of doc 16: exactly what the python chain
// clustered).
local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local params = import 'pgrapher/experiment/fdvd/params.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);
local clus_maker = import 'pgrapher/experiment/fdvd/clus.jsonnet';
local flash = import 'pgrapher/experiment/fdvd/flash.jsonnet';

function(indir, library, calibration, outfile='', rootfile='', opwf='', opflash='', drift_in='', drift_model='', frames='',
         drift_device='cpu', dump_crops=false, pipe='simple', iso_cm=30, event=0, bee_zip='', dump_tables=false,
         anode_indices=std.range(0, std.length(tools.anodes) - 1))
  assert (opwf != '') != (opflash != '') : 'give exactly one of opwf / opflash';
  local anodes = [tools.anodes[i] for i in anode_indices];
  local n = std.length(anodes);
  local C = clus_maker(params, anodes);

  local src = [C.source(a, indir, pipe) for a in anodes];
  local fan = [g.pnode({ type: 'ClusterFanout', name: 'crm%d' % a.data.ident, data: { multiplicity: 2 } },
                       nin=1, nout=2) for a in anodes];
  local ptb = [C.point_tree(a) for a in anodes];
  local merge = C.merge(n);
  local mabc = C.mabc(iso_cm=iso_cm, event=event, bee_zip=bee_zip);
  local btab = g.pnode({
    type: 'FdvdBlobTable',
    name: 'fdvd-all',
    // python blob order: archives sorted by name (img_eval_all / sigbkg_eval.event_blobs glob)
    data: { multiplicity: n, labels: [a.name for a in anodes], tick: 0.5 * wc.us },
  }, nin=n, nout=1);

  local light = if opwf != '' then g.pipeline([
    g.pnode({ type: 'TensorFileSource', name: 'opwf', data: { inname: opwf, prefix: 'opwf_' } }, nin=0, nout=1),
    g.pnode({ type: 'FdvdOpHitFinder', name: '', data: {} }, nin=1, nout=1),
    g.pnode({ type: 'OpFlashFinder', name: 'fdvd', data: flash }, nin=1, nout=1),
  ]) else g.pnode({ type: 'TensorFileSource', name: 'opflash', data: { inname: opflash, prefix: 'opflash_' } },
                  nin=0, nout=1);

  assert !(drift_in != '' && drift_model != '') : 'give drift_in or drift_model, not both';
  assert (drift_model == '') == (frames == '') : 'drift_model needs frames';
  local with_drift = drift_in != '';
  local torch = {
    type: 'TorchTensorSetService',
    name: 'drift',
    data: { model: drift_model, device: drift_device },
  };
  local ql = g.pnode({
    type: 'FdvdLowEQLMatching',
    name: 'fdvd',
    data: {
      multiplicity: if with_drift then 4 else 3,
      inpath: 'pointtrees/%d',
      library: library,
      calibration: calibration,
      geom_file: 'pgrapher/experiment/fdvd/fdvd-opdet-geom.json',
      dump_tables: dump_tables,
    } + (if drift_model != '' then {
      drift: { forward: wc.tn(torch), frames: frames, frame_tag: 'gauss', dump_crops: dump_crops },
    } else {}),
  }, nin=if with_drift then 4 else 3, nout=1, uses=if drift_model != '' then [torch] else []);
  local drift = g.pnode({ type: 'TensorFileSource', name: 'drift', data: { inname: drift_in, prefix: 'drift_' } },
                        nin=0, nout=1);
  assert outfile != '' || rootfile != '' : 'give outfile and/or rootfile';
  local sink = g.pnode({ type: 'TensorFileSink', name: 'lowe',
                         data: { outname: if outfile != '' then outfile else 'unused.tar.gz', prefix: 'lowe_',
                                 dump_mode: outfile == '' } }, nin=1, nout=0);
  local writer = g.pnode({ type: 'FdvdLowERootWriter', name: 'lowe', data: { output_filename: rootfile, event: event } },
                         nin=1, nout=1);
  local tail = if rootfile != '' then g.pipeline([writer, sink]) else sink;

  local edges =
    [g.edge(src[k], fan[k]) for k in std.range(0, n - 1)] +
    [g.edge(fan[k], ptb[k], 0, 0) for k in std.range(0, n - 1)] +
    [g.edge(fan[k], btab, 1, k) for k in std.range(0, n - 1)] +
    [g.edge(ptb[k], merge, 0, k) for k in std.range(0, n - 1)] +
    [g.edge(merge, mabc), g.edge(mabc, ql, 0, 0), g.edge(btab, ql, 0, 1), g.edge(light, ql, 0, 2), g.edge(ql, tail)] +
    (if with_drift then [g.edge(drift, ql, 0, 3)] else []);
  local graph = g.intern(innodes=src + [light] + (if with_drift then [drift] else []), centernodes=fan + ptb + [merge, mabc, btab, ql],
                         outnodes=[tail], edges=edges);
  local app = { type: 'Pgrapher', data: { edges: g.edges(graph) } };
  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellSio', 'WireCellAux', 'WireCellImg', 'WireCellClus',
                'WireCellFlash', 'WireCellMatch'] + (if drift_model != '' then ['WireCellPytorch'] else [])
               + (if rootfile != '' then ['WireCellRoot'] else []),
      apps: ['Pgrapher'],
    },
  };
  [cmdline] + g.uses(graph) + [app]
