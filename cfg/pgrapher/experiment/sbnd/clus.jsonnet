// Canonical SBND per-APA and all-APA clustering using MultiAlgBlobClustering.
//
// This is the single source of truth for the SBND clustering graph and fiducial
// volume.  Both callers import it:
//   - cfg/.../sbnd/wcls-img-clus.jsonnet  (LArSoft production) via per_volume()
//   - sbnd_xin/wct-clustering.jsonnet     (standalone dev chain) via per_apa()/all_apa()
// The standalone sbnd_xin/clus.jsonnet is a thin re-export of this file.
//
// Schema: DetectorVolumes + PCTransformSet + pgrapher/common/clus.jsonnet
// clustering_methods() components wired into MABC's "pipeline".  (The older
// SimpleClusGeomHelper + func_cfgs style is dead against current MABC, which
// reads cfg["pipeline"] as a list of component type:name strings and requires
// cfg["detector_volumes"].)

local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local clus = import 'pgrapher/common/clus.jsonnet';
local dead_regions = import 'pgrapher/experiment/sbnd/dead_regions.jsonnet';

local time_offset = -205 * wc.us;  // = -tick0_time (cfg/.../sbnd/params.jsonnet sim.tick0_time)
local drift_speed = 1.563 * wc.mm / wc.us;
local bee_dir = 'data';

local common_coords = ['x', 'y', 'z'];

// Per-TPC transverse (Y,Z) position offset, materialized in the post-QLMatching
// scope by T0Correction as y_cor/z_cor (see match/docs/cathode-offset-correction.md).
// One flag drives BOTH the metadata injection (which the C++ keys on for the
// y_cor/z_cor scope) and the corrected-coords names below, so jsonnet and C++ stay
// in lockstep.  pos_offset is a DATA-ONLY calibration: it was measured from data
// cathode-crossers (data transverse ~1.4 cm vs MC ~0, see
// match/docs/cathode-offset-correction.md / project cathode-crossing diagnosis).
// MC has no such misalignment, so applying it to MC would inject a spurious shift.
// Hence pos_offset_on is gated on reality at the function entry below
// (reality='data' -> on, 'sim' -> off).  x component is 0 (drift stays with the
// t0/flash_x_offset term).  Values = symmetric split of the measured
// T_yz=(-0.22,+1.34) cm cathode gap.
local pos_offset_a0 = [0, -0.11 * wc.cm, 0.67 * wc.cm];   // TPC0 (East, x<0)
local pos_offset_a1 = [0, 0.11 * wc.cm, -0.67 * wc.cm];   // TPC1 (West, x>=0)

local common_corr_coords(pos_offset_on) =
    if pos_offset_on then ['x_t0cor', 'y_cor', 'z_cor'] else ['x_t0cor', 'y', 'z'];

// SBND cathode-crossing connector: connect the two halves of a cathode-crossing
// cosmic that the generic all-APA merge passes leave unmerged (their closest-point
// distance lands just over the 3 cm lenient-merge cap; see
// clus/docs/cathode-crossing-clustering.md).  Narrow cathode-specific cut set
// (collinear + close + opposite TPCs + both ends at the cathode), so it cannot fire
// within a single TPC.  SBND committed ON; set false to recover the pre-connector
// all-APA pipeline (bit-identical).  Retire (flip false) when the pos_offset / SCE
// transverse calibration tightens enough that the generic 3 cm path catches these.
local cathode_connect_on = true;

// FV = sbnd-wires-geometry-v0206 bbox - 1 cm inset on every face.
// X anode = W (collection) plane; X inner = data CPA face (DENT-gap geometry, ±1.5 cm).
// See wire-cell-bee3/docs/sbnd_geometry.md §8.1 for the source TPC bounding boxes.
local dvm(pos_offset_on) = {
    overall: {
        FV_xmin: -201.05  * wc.cm,  // W plane (-202.05) + 1 cm inset
        FV_xmax:  201.05  * wc.cm,  // W plane (+202.05) - 1 cm
        FV_ymin: -199.312 * wc.cm,  // wires Y bbox (-200.312) + 1 cm
        FV_ymax:  199.312 * wc.cm,  // wires Y bbox (+200.312) - 1 cm
        FV_zmin:    0.85  * wc.cm,  // wires Z bbox (-0.15) + 1 cm; legacy 4.05 was stale v02_02 z_J frame
        FV_zmax:  500.15  * wc.cm,  // wires Z bbox (+501.15) - 1 cm
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
        FV_ymin_margin: 2.5 * wc.cm,
        FV_ymax_margin: 2.5 * wc.cm,
        FV_zmin_margin: 3 * wc.cm,
        FV_zmax_margin: 3 * wc.cm,
        vertical_dir: [0, 1, 0],
        beam_dir: [0, 0, 1],
    },
    a0f0pA: {
        drift_speed: drift_speed,
        tick: 0.5 * wc.us,
        tick_drift: self.drift_speed * self.tick,
        time_offset: time_offset,
        nticks_live_slice: 4,
        FV_xmin: -201.05 * wc.cm,  // W plane (-202.05) + 1 cm
        FV_xmax:   -2.5  * wc.cm,  // data CPA face (-1.5) - 1 cm toward TPC0 interior
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
        // y/z fiducial bounds + margins are detector-wide (SBND TPCs span the full
        // height/length), so the per-(APA,face) FV reuses the overall values; only
        // the x-bounds above are genuinely per-TPC.  These complete the per-face FV
        // so the scope-aware select_scope_fv (clustering_separate / clustering_neutrino)
        // has y/z + dirs for a single-APA pass without falling back to "overall".
        FV_ymin: $.overall.FV_ymin,
        FV_ymax: $.overall.FV_ymax,
        FV_zmin: $.overall.FV_zmin,
        FV_zmax: $.overall.FV_zmax,
        FV_ymin_margin: $.overall.FV_ymin_margin,
        FV_ymax_margin: $.overall.FV_ymax_margin,
        FV_zmin_margin: $.overall.FV_zmin_margin,
        FV_zmax_margin: $.overall.FV_zmax_margin,
    } + (if pos_offset_on then { pos_offset: pos_offset_a0 } else {}),
    a1f0pA: $.a0f0pA + {
        FV_xmin:    2.5  * wc.cm,  // data CPA face (+1.5) + 1 cm toward TPC1 interior
        FV_xmax:  201.05 * wc.cm,  // W plane (+202.05) - 1 cm
    } + (if pos_offset_on then { pos_offset: pos_offset_a1 } else {}),  // override a0's
};

local anodes_name(anodes, face='') =
    std.join('-', [std.toString(a.data.ident) for a in anodes])
    + if face == '' then '' else '-' + std.toString(face);

local detector_volumes(anodes, face='', pos_offset_on=true) = {
    local m = dvm(pos_offset_on),
    type: 'DetectorVolumes',
    name: 'dv-apa' + anodes_name(anodes, face),
    data: {
        anodes: [wc.tn(anode) for anode in anodes],
        metadata:
            { overall: m['overall'] } +
            { a0f0pA: m['a0f0pA'] } +
            { a1f0pA: m['a1f0pA'] },
    },
    uses: anodes,
};

local pctransforms(dv) = {
    type: 'PCTransformSet',
    name: dv.name,
    data: { detector_volumes: wc.tn(dv) },
    uses: [dv],
};

local bs_live_face(apa, face) = {
    type: 'BlobSampler',
    name: 'live-%s-%d' % [apa, face],
    data: {
        drift_speed: drift_speed,
        time_offset: time_offset,
        strategy: ['stepped'],
        extra: ['.*wire_index', '.*charge_val', '.*charge_unc', 'wpid'],
    },
};
local bs_dead_face(apa, face) = {
    type: 'BlobSampler',
    name: 'dead-%s-%d' % [apa, face],
    data: {
        strategy: ['center'],
        extra: ['.*'],
    },
};

// Per-APA / per-face clustering.  The active pipeline matches the sbnd_xin
// standalone chain (pointed .. connect1).  The original cfg func_cfgs tail
// (deghost -> examine_x_boundary -> isolated) is retained below as commented
// lines so it can be re-enabled without re-deriving it.
local clus_per_face(anode, face, dump, output_dir, runNo, subRunNo, eventNo, bee_sink=null, rse_from_ident=false, pos_offset_on=true) = {
    local dv = detector_volumes([anode], face, pos_offset_on),
    local pcts = pctransforms(dv),
    local bsl = bs_live_face(anode.name, face),
    local bsd = bs_dead_face(anode.name, face),
    local ptb = g.pnode({
        type: 'PointTreeBuilding',
        name: '%s-%d' % [anode.name, face],
        data: {
            samplers: { '3d': wc.tn(bsl), dead: wc.tn(bsd) },
            multiplicity: 2,
            tags: ['live', 'dead'],
            anode: wc.tn(anode),
            face: face,
            detector_volumes: wc.tn(dv),
            // Hand-declared dead winds at the known-bad Y-Z region (W dead + U/V
            // distorted) so examine_bundles' relaxed-graph bridge crosses it
            // instead of fragmenting a single track.  Per-anode (TPC0/TPC1).
            inject_dead_winds: [dead_regions.region(anode.data.ident)],
        },
    }, nin=2, nout=1, uses=[bsl, bsd, dv]),
    local cluster2pct = ptb,
    local face_name = '%s-%d' % [anode.name, face],
    local cm = clus.clustering_methods(
        prefix=face_name,
        detector_volumes=dv,
        pc_transforms=pcts,
        coords=common_coords),
    local cm_pipeline = [
        cm.pointed(),
        cm.live_dead(dead_live_overlap_offset=2),
        cm.extend(flag=4, length_cut=60 * wc.cm, num_try=0, length_2_cut=15 * wc.cm, num_dead_try=1),
        cm.regular(name='-one', length_cut=60 * wc.cm, flag_enable_extend=false),
        cm.regular(name='_two', length_cut=30 * wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35 * wc.cm),
        cm.close(length_cut=1.2 * wc.cm),
        cm.extend_loop(num_try=3),
        // Raise the convex-hull point cap (default 10000) so full-detector
        // multi-track overclusters (>10k points) are still considered for
        // separation; otherwise get_hull returns empty and separation is skipped.
        cm.separate(use_ctpc=true, max_hull_points=100000, sbnd_boundary_tag=true),
        // SBND: cap the isochronous-relaxed connection on the real closest-point
        // distance.  Without it, connect1 merges two genuinely-separate isochronous
        // cosmics (e.g. evt 183888, ~7.3 cm apart in drift) on the misleadingly small
        // infinite-line distance.  5 cm < the 7.3 cm real gap, above SBND broken-track
        // gaps.  Default OFF (-1) elsewhere keeps production bit-identical.
        cm.connect1(iso_max_dis=5 * wc.cm),
        // MicroBooNE-style clustering tail: produce cluster groups (one main +
        // associated small clusters) carried as the "isolated"/"perblob" per-blob
        // array (main blobs tagged -1). examine_bundles MUST follow neutrino/isolated,
        // which produce the perblob it consumes. Group-aware QLMatching downstream
        // decomposes these groups into main+others for prototype-style bundle building.
        cm.deghost(),
        cm.examine_x_boundary(),
        cm.protect_overclustering(),
        cm.neutrino(),
        // SBND: tighten the isolated small/big length_cut from the 20 cm default
        // to 15 cm so a ~16 cm EM (gamma) blob is no longer auto-classified
        // "small" and absorbed into a nearby long cosmic track by the
        // angle-less 80 cm small->big merge.  See sbnd_xin/docs/
        // overclustering-evt11-gamma.md.  range_cut left at its 150 default.
        cm.isolated(length_cut=15 * wc.cm),
        cm.examine_bundles(),
    ],
    local bee_zip_path = (if output_dir == '' then '' else output_dir + '/')
                         + 'mabc-%s-face%d.zip' % [anode.name, face],
    local mabc = g.pnode({
        local name = '%s-%d' % [anode.name, face],
        type: 'MultiAlgBlobClustering',
        name: name,
        data: {
            inpath: 'pointtrees/%d',
            outpath: 'pointtrees/%d',
            perf: true,
            bee_dir: bee_dir,
            bee_zip: bee_zip_path,
            // When a shared Bee sink is supplied, all Bee writes go to that one
            // zip (bee_zip is then ignored) and the per-APA dead-area is dropped
            // to avoid a duplicate channel-deadarea-* entry (the all-APA node
            // writes byte-identical dead-area for both APAs into the same zip).
            [if bee_sink != null then 'bee_sink']: wc.tn(bee_sink),
            bee_detector: 'sbnd',
            initial_index: 0,
            use_config_rse: true,
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            // Take the Bee event number from each event's tensor ident (run/subrun
            // = 0).  Conditional key: omitted when off, so production stays
            // byte-identical (mirrors the bee_sink conditional above).
            [if rse_from_ident then 'rse_from_ident']: true,
            save_deadarea: bee_sink == null,
            dead_area_version: 2,
            anodes: [wc.tn(anode)],
            face: face,
            detector_volumes: wc.tn(dv),
            // Renumber cluster idents (insertion order, 1..N) after every
            // clustering step.  Without this, clusters created mid-pipeline keep
            // the unset-ident sentinel (-1) and collapse together in the Bee
            // display.  See clustering.md "Cluster id numbering".
            cluster_id_order: 'tree',
            bee_points_sets: [{
                name: 'clustering',
                detector: 'sbnd',
                algorithm: 'clustering',
                pcname: '3d',
                coords: ['x', 'y', 'z'],
                individual: true,
            }],
            pipeline: wc.tns(cm_pipeline),
        },
    }, nin=1, nout=1, uses=[dv, anode, pcts] + cm_pipeline
              + (if bee_sink != null then [bee_sink] else [])),
    local sink = g.pnode({
        type: 'TensorFileSink',
        name: 'clus_per_face-%s-%d' % [anode.name, face],
        data: {
            outname: 'trash-%s-face%d.tar.gz' % [anode.name, face],
            prefix: 'clustering_',
            dump_mode: true,
        },
    }, nin=1, nout=0),
    local end = if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),
    ret:: g.pipeline([cluster2pct, end], 'clus_per_face-%s-%d' % [anode.name, face]),
}.ret;

// premerged=true: the upstream node (joint QLMatching) has already merged the
// per-APA cluster trees into one, so skip the PointTreeMerging fanin and feed the
// single pre-merged input straight to the all-APA MABC.  Default false = the
// historical per-APA path (two QLMatching nodes -> PointTreeMerging -> MABC).
// tensor_outname (default ''): when set, the terminal TensorFileSink becomes a
// REAL sink writing the post-QL point-cloud tree tensors (live+dead, prefix
// 'clustering_') to that file -- the persistent intermediate format consumed by
// the downstream pattern-recognition job (see sbnd/docs/sbnd-pattern-recognition.md).
// Default '' keeps the historical dump_mode no-op sink (byte-identical).
local clus_all_apa(anodes, dump, output_dir, runNo, subRunNo, eventNo, bee_sink=null, premerged=false, rse_from_ident=false, pos_offset_on=true, tensor_outname='') = {
    local nanodes = std.length(anodes),
    local pcmerging = g.pnode({
        type: 'PointTreeMerging',
        name: 'clus_all_apa',
        data: {
            multiplicity: nanodes,
            inpath: 'pointtrees/%d',
            outpath: 'pointtrees/%d',
            // Carry the per-APA self-contained optical flash display PC
            // (written by QLMatching, keyed by global flash id) into the merged
            // all-APA grouping so the all-APA MABC can dump the op/flash Bee
            // display. Other per-anode root PCs are intentionally not merged.
            root_pcs_to_merge: ['opflash'],
        },
    }, nin=nanodes, nout=1),
    local dv = detector_volumes(anodes, '', pos_offset_on),
    local pcts = pctransforms(dv),
    local cm_old = clus.clustering_methods(
        prefix='all', detector_volumes=dv, pc_transforms=pcts, coords=common_coords),
    local cm = clus.clustering_methods(
        prefix='all', detector_volumes=dv, pc_transforms=pcts, coords=common_corr_coords(pos_offset_on)),
    // Combined (all-APA) clustering runs AFTER QL charge-light matching, so every
    // cluster carries a matched flash time (cluster_t0).  switch_scope applies the
    // per-cluster T0 correction (x_t0cor scope) and drops any stale per-APA
    // "isolated"/"perblob" array (it destroys+recreates every cluster).  The
    // merging steps below set use_flash_t0=true so they only merge clusters
    // coincident in flash time, and examine_bundles collapses each flash group into
    // one cluster carrying a fresh "isolated"/"perblob" array (main sub-component
    // = -1), like clustering_isolated but grouped by flash time instead of geometry.
    local cm_pipeline = [
        cm_old.switch_scope(),
        cm.extend(flag=4, length_cut=60 * wc.cm, num_try=0, length_2_cut=15 * wc.cm, num_dead_try=1, use_flash_t0=true),
        cm.regular(name='1', length_cut=60 * wc.cm, flag_enable_extend=false, use_flash_t0=true),
        cm.regular(name='2', length_cut=30 * wc.cm, flag_enable_extend=true, use_flash_t0=true),
        cm.parallel_prolong(length_cut=35 * wc.cm, use_flash_t0=true),
        cm.close(length_cut=1.2 * wc.cm, use_flash_t0=true),
        cm.extend_loop(num_try=3, use_flash_t0=true),
    ]
    // Cathode-crossing connector: after the generic merge passes (so it only ADDS
    // merges they missed) and before examine_bundles (so a connected crosser is one
    // cluster before the flash-bundle collapse).  SBND-on; off => list unchanged.
    + (if cathode_connect_on then [cm.cathode_connect(cathode_x_cut=5*wc.cm, drift_cut=8*wc.cm, min_length_short=2*wc.cm, short_dir_len=25*wc.cm, conn_short_cut=30.0, flash_t0_window=800*wc.ns)] else [])
    + [
        // flags_from_longest: the flash-time merge here collapses a bundle's
        // clusters into one; without this the merged cluster inherits its flags
        // from an arbitrary member, so a matched main that absorbs a tiny
        // co-merged fragment loses flag_main_cluster to it (SBND evt284349:
        // the 2173-pt beam track lost it to a 3-pt TPC1 speck, leaving the flag
        // only on its own out-of-volume shard).  The taggers key on that flag.
        cm.examine_bundles(use_flash_t0=true, flags_from_longest=true),
    ],
    local bee_zip_path = (if output_dir == '' then '' else output_dir + '/') + 'mabc-all-apa.zip',
    local mabc = g.pnode({
        type: 'MultiAlgBlobClustering',
        name: 'clus_all_apa',
        data: {
            inpath: 'pointtrees/%d',
            outpath: 'pointtrees/%d',
            perf: true,
            bee_dir: bee_dir,
            bee_zip: bee_zip_path,
            // When a shared Bee sink is supplied, the all-APA views (img/
            // clustering/op + dead-area for both APAs) are written into that one
            // shared zip together with the per-APA views; bee_zip is ignored.
            [if bee_sink != null then 'bee_sink']: wc.tn(bee_sink),
            bee_detector: 'sbnd',
            initial_index: 0,
            use_config_rse: true,
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            // Take the Bee event number from each event's tensor ident (run/subrun
            // = 0).  Conditional key: omitted when off, so production stays
            // byte-identical (mirrors the bee_sink conditional above).
            [if rse_from_ident then 'rse_from_ident']: true,
            save_deadarea: true,
            dead_area_version: 2,
            // Dump the optical flash / charge-light "op" display into this same
            // mabc-all-apa.zip (reads the merged-root "opflash" PC + per-cluster
            // matched-flash association written by QLMatching). bee_detector
            // above supplies the Bee geom.
            save_opflash: true,
            // Emit one "op" row per flash carrying ALL matched cluster ids
            // (op_cluster_ids array) with summed predicted PE, so a flash matched
            // to several clusters shows them together (MicroBooNE-style).
            bee_flash_per_flash: true,
            // Group flashes from the two TPC sides by this ±time window and stash
            // a per-flash "group" array on the root opflash PC (pre-pipeline, so
            // the op dump and every later step can read it).  0 = off.
            flash_group_window: 80 * wc.ns,
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            // Renumber cluster idents (insertion order, 1..N) after every step;
            // see clus_per_face above and clustering.md "Cluster id numbering".
            cluster_id_order: 'tree',
            bee_points_sets: [
                {
                    name: 'img',
                    detector: 'sbnd',
                    algorithm: 'img',
                    pcname: '3d',
                    coords: ['x', 'y', 'z'],
                    individual: false,
                },
                {
                    name: 'clustering',
                    detector: 'sbnd',
                    algorithm: 'clustering',
                    pcname: '3d',
                    // Same corrected coords as the clustering scope, so the Bee
                    // display reflects the transverse shift when it is on (makes the
                    // separate Bee-zip transverse shift redundant -- pick one).
                    coords: common_corr_coords(pos_offset_on),
                    individual: false,
                },
            ],
            pipeline: wc.tns(cm_pipeline),
        },
    }, nin=1, nout=1, uses=anodes + [dv, pcts] + cm_pipeline
              + (if bee_sink != null then [bee_sink] else [])),
    local sink = g.pnode({
        type: 'TensorFileSink',
        name: 'clus_all_apa',
        data: {
            outname: if tensor_outname == '' then 'trash-all-apa.tar.gz' else tensor_outname,
            prefix: 'clustering_',
            dump_mode: tensor_outname == '',
        },
    }, nin=1, nout=0),
    local end = if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),
    // premerged: input is already one merged tree (joint QLMatching) -> feed MABC
    // directly, no PointTreeMerging.  Else: fan the per-APA inputs into pcmerging.
    ret:: if premerged then end else g.intern(
        innodes=[pcmerging],
        centernodes=[],
        outnodes=[end],
        edges=[g.edge(pcmerging, end, 0, 0)]
    ),
}.ret;

// Pattern-recognition (PR) stage: consume the persisted post-QL point-cloud tree
// (the tensor_outname tarball written by clus_all_apa, reloaded through a
// TensorFileSource) and run the PR-tail visitors on it.  See
// sbnd/docs/sbnd-pattern-recognition.md.  pipeline_names selects the visitors by
// name from the map below (empty = pass-through, used by the round-trip identity
// gate).  The Bee config mirrors clus_all_apa so the 'clustering' layer of the
// PR job's zip is directly comparable to mabc-all-apa.zip, except save_opflash
// is off: the op display needs the per-cluster flashpred pcarray, which is
// consumed by the Q/L job's pre-pipeline op dump and is not in the tarball.
local clus_pr(anodes, dump, output_dir, runNo, subRunNo, eventNo, rse_from_ident=false, pos_offset_on=true, pipeline_names=[], tensor_outname='',
              trackfitting_config_file='', particle_dataset=null, extra_uses=[], dl_weights='', beam_window=[0, 0],
              tgm_neutrino_candidate=false, tgm_chord_charge=false,
              tgm_chord_mode='chord', tgm_component_extremes=false) = {
    local dv = detector_volumes(anodes, '', pos_offset_on),
    local pcts = pctransforms(dv),
    // DetectorVolumes implements IFiducial (box FV from its metadata) -- used by
    // MakeFiducialUtils / the taggers' inside_fiducial_volume().
    local cm_old = clus.clustering_methods(
        prefix='pr', detector_volumes=dv, pc_transforms=pcts, fiducial=dv, coords=common_coords),
    local cm = clus.clustering_methods(
        prefix='pr', detector_volumes=dv, pc_transforms=pcts, fiducial=dv, coords=common_corr_coords(pos_offset_on)),
    // Box-model recombination at the SBND drift field (uBooNE used 0.273 kV/cm).
    local sbnd_box_recomb = {
        type: 'BoxRecombination',
        name: 'sbnd_box_recomb',
        data: { A: 1.0, B: 0.255, Efield: 0.5, rho: 1.38, Wi: 23.6e-6 },
    },
    // TGM fiducial: ONE box spanning BOTH TPCs (the overall FV bounds of
    // dvm above), so a cathode-crossing track is not an "exiter" at x=0.
    // The default fiducial=dv cannot serve here: DetectorVolumes::contained()
    // is the union of per-face sensitive volumes, which excludes the CPA slab
    // (|x| < 0.45 cm).  Margins go in via the tagger's fv_tolerance instead
    // of the box, mirroring the metadata *_margin values.
    local sbnd_pr_fv = {
        type: 'BoxFiducial',
        name: 'sbnd_pr_fv',
        data: {
            bounds: {
                tail: { x: -201.05 * wc.cm, y: -199.312 * wc.cm, z: 0.85 * wc.cm },
                head: { x: 201.05 * wc.cm, y: 199.312 * wc.cm, z: 500.15 * wc.cm },
            },
        },
    },
    local sbnd_pr_fv_margins = [-2 * wc.cm, -2 * wc.cm, -2.5 * wc.cm, -2.5 * wc.cm, -3 * wc.cm, -3 * wc.cm],
    // Retiler for the steiner stage: same 'stepped' samplers that built the 3d
    // PC (PointTreeBuilding), one per (APA, face 0).
    local improve2 = cm.improve_cluster_2(
        anodes=anodes,
        samplers=[clus.sampler(bs_live_face(a.name, 0), apa=a.data.ident, face=0) for a in anodes]),
    // Visitors available to the PR pipeline, by name.  switch_scope re-applies
    // the per-cluster T0 correction on the loaded tree (the corrected scope is
    // runtime state and does not persist through the tarball); it recomputes
    // deterministically from cluster_t0.  fiducialutils MUST precede any
    // tagger (they silently no-op without it).
    local cm_by_name = {
        switch_scope: cm_old.switch_scope(),
        // SBND has no beam_flash flag (QLMatching sets main/associated_cluster
        // instead) -- process every scope-passing cluster.
        steiner: cm.steiner(retiler=improve2, perf=true, require_beam_flash=false),
        fiducialutils: cm.fiducialutils(),
        tagger_check_stm: cm.tagger_check_stm(
            trackfitting_config_file=trackfitting_config_file,
            particle_dataset=wc.tn(particle_dataset),
            recombination_model=wc.tn(sbnd_box_recomb),
            require_in_scope=true),
        // Through-going-muon tagger (prototype check_tgm port).  Runs on every
        // matched main cluster.  tgm_neutrino_candidate (C++ default false;
        // key omitted when off => byte-identical): enable the ported
        // check_neutrino_candidate veto so in-beam-window bundles may be
        // tagged; when off in-beam bundles are never tagged (conservative v1
        // behavior).  Must run after fiducialutils (dead-region /
        // signal-processing checks) and before tagger_check_stm (which skips
        // TGM-flagged mains).
        // require_in_scope: evaluate only clusters that pass switch_scope's
        // active-volume filter.  Without it the tagger also sees the
        // out-of-volume shards switch_scope splits off (which keep an inherited
        // flag_main_cluster) -- on the SBND 10-event reco1 sample those were 43%
        // of all evaluated mains and tagged TGM at 85% vs 44% for real clusters.
        tagger_check_tgm: cm.tagger_check_tgm(
            fiducial=wc.tn(sbnd_pr_fv),
            fv_tolerance=sbnd_pr_fv_margins,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1],
            require_in_scope=true,
            check_neutrino_candidate=tgm_neutrino_candidate,
            require_chord_charge=tgm_chord_charge,
            // C++ default "chord". Key omitted then => byte-identical.
            chord_charge_mode=tgm_chord_mode,
            component_extremes=tgm_component_extremes),
        // Fully-contained tagger.  Independent of TGM/STM: it evaluates every
        // in-scope main cluster and only records a containment verdict (flag
        // "FC"), so it neither vetoes nor is vetoed by them.  Placed LAST in
        // the pipeline on purpose -- cluster_fc_check lazily populates
        // PCA/hough/steiner-boundary caches, so running it after the other
        // taggers leaves their inputs pristine rather than relying on those
        // caches being order-independent.  Coverage is unaffected by position.
        // fiducial/fv_tolerance are the SAME objects tagger_check_tgm gets, so
        // "contained" means one thing across both verdicts.  Without them FC
        // fell back to FiducialUtils (per-face sensitive volumes, no margin),
        // which is both more permissive at every wall than TGM's inset box and
        // holed at the CPA slab -- see docs/27_fc-tgm-consistent-fv.md.
        tagger_check_fc: cm.tagger_check_fc(
            fiducial=wc.tn(sbnd_pr_fv),
            fv_tolerance=sbnd_pr_fv_margins,
            require_in_scope=true),
        // Neutrino pattern recognition on the beam-coincident bundle.  The
        // beam_window gate (on cluster_t0 = matched flash time) replaces
        // uBooNE's single-main + beam_flash selection; companions are the
        // associated clusters sharing the main's matched_flash_gid.
        // dl_weights='' = geometric vertex (the SCN net is uBooNE-trained).
        tagger_check_neutrino: cm.tagger_check_neutrino(
            trackfitting_config_file=trackfitting_config_file,
            particle_dataset=wc.tn(particle_dataset),
            recombination_model=wc.tn(sbnd_box_recomb),
            perf=true,
            dl_weights=dl_weights,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1]),
    },
    local cm_pipeline = [cm_by_name[n] for n in pipeline_names],
    // The taggers' configs only name the recombination/particle-dataset
    // components; emit them (and the caller's LinterpFunctions etc. via
    // extra_uses) when a tagger is in the pipeline.
    local tagger_uses = (if std.member(pipeline_names, 'tagger_check_stm')
                         || std.member(pipeline_names, 'tagger_check_neutrino')
                         then [sbnd_box_recomb] + extra_uses else [])
                        + (if std.member(pipeline_names, 'tagger_check_tgm')
                           || std.member(pipeline_names, 'tagger_check_fc')
                           then [sbnd_pr_fv] else []),
    local bee_zip_path = (if output_dir == '' then '' else output_dir + '/') + 'mabc-pr.zip',
    local mabc = g.pnode({
        type: 'MultiAlgBlobClustering',
        name: 'clus_pr',
        data: {
            inpath: 'pointtrees/%d',
            outpath: 'pointtrees/%d',
            perf: true,
            bee_dir: bee_dir,
            bee_zip: bee_zip_path,
            bee_detector: 'sbnd',
            initial_index: 0,
            use_config_rse: true,
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            [if rse_from_ident then 'rse_from_ident']: true,
            save_deadarea: true,
            dead_area_version: 2,
            save_opflash: false,
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            cluster_id_order: 'tree',
            bee_points_sets: [
                {
                    name: 'clustering',
                    detector: 'sbnd',
                    algorithm: 'clustering',
                    pcname: '3d',
                    coords: common_corr_coords(pos_offset_on),
                    individual: false,
                    // Same filter semantics as clus_all_apa (default 1): the dump
                    // keys on the per-cluster runtime scope-filter flag, which is
                    // NOT persisted through the tarball -- it is re-established by
                    // running switch_scope at the head of the PR pipeline.  A
                    // pass-through job (pipeline_names=[]) therefore dumps
                    // nothing; the round-trip identity gate runs
                    // pipeline_names=['switch_scope'].
                },
                // Neutrino-PR layers, dumped after TaggerCheckNeutrino runs (the
                // visitor key is the full type:name; prefix='pr').  Inert unless
                // tagger_check_neutrino is in pipeline_names.  PRGraph fit points
                // are already T0-corrected, hence plain x/y/z.
                {
                    name: 'track_fit',
                    visitor: 'TaggerCheckNeutrino:pr',
                    grouping: 'live',
                    detector: 'sbnd',
                    algorithm: 'track_fit',
                    pcname: '3d',            // not used for PRGraph dumps, but required
                    coords: ['x', 'y', 'z'],
                    individual: false,
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                },
                {
                    name: 'shower_track',    // associated points: q=15000 shower, q=0 track
                    visitor: 'TaggerCheckNeutrino:pr',
                    grouping: 'live',
                    detector: 'sbnd',
                    algorithm: 'shower_track',
                    pcname: '3d',
                    coords: ['x', 'y', 'z'],
                    individual: false,
                    use_associate_points: true,
                },
                {
                    name: 'vertices',        // PR graph vertices; main vertex q=15000
                    visitor: 'TaggerCheckNeutrino:pr',
                    grouping: 'live',
                    detector: 'sbnd',
                    algorithm: 'vertices',
                    pcname: '3d',
                    coords: ['x', 'y', 'z'],
                    individual: false,
                    use_graph_vertices: true,
                },
            ],
            // Particle-flow Bee output ("mc" jsTree JSON), emitted once after
            // TaggerCheckNeutrino runs; inert when the visitor is not in the pipeline.
            bee_pf: [
                {
                    name: 'mc',
                    visitor: 'TaggerCheckNeutrino:pr',
                    grouping: 'live',
                },
            ],
            pipeline: wc.tns(cm_pipeline),
        },
    }, nin=1, nout=1, uses=anodes + [dv, pcts] + cm_pipeline + tagger_uses),
    local sink = g.pnode({
        type: 'TensorFileSink',
        name: 'clus_pr',
        data: {
            outname: if tensor_outname == '' then 'trash-pr.tar.gz' else tensor_outname,
            prefix: 'clustering_',
            dump_mode: tensor_outname == '',
        },
    }, nin=1, nout=0),
    ret:: if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),
}.ret;

// rse_from_ident (default false): when true, the MultiAlgBlobClustering nodes take
// the Bee event number from each event's tensor ident (run/subrun = 0) instead of
// the configured runNo/eventNo auto-increment.  Used by the bundled standalone chain
// (one wire-cell call over many events) whose ident already carries the real event
// id.  Default false keeps production byte-identical (the key is omitted).
function(output_dir='.', runNo=0, subRunNo=0, eventNo=0, rse_from_ident=false, reality='data') {
    // pos_offset (per-TPC transverse y,z calibration) is data-only -- see the
    // pos_offset_a0/a1 comment above.  reality='data' (default; keeps the data
    // chain byte-identical to the previous always-on state) -> on; reality='sim'
    // (MC) -> off, so the MC chain carries no transverse shift.
    local pos_offset_on = reality == 'data',
    // bee_sink (default null): when set to a shared IBeeSink node, all Bee
    // output for this node goes into that single shared zip instead of this
    // node's own bee_zip.  Default null -> own zip (production byte-identical).
    per_face(anode, face=0, dump=true, bee_sink=null)::
        clus_per_face(anode, face=face, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on),
    per_apa(anode, dump=true, bee_sink=null)::
        clus_per_face(anode, face=0, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on),
    // Production (LArSoft) entry point used by wcls-img-clus.jsonnet.
    per_volume(anode, face=0, dump=true, bee_sink=null)::
        clus_per_face(anode, face=face, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on),
    all_apa(anodes, dump=true, bee_sink=null, premerged=false, tensor_outname='')::
        clus_all_apa(anodes, dump=dump,
                     output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                     bee_sink=bee_sink, premerged=premerged, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                     tensor_outname=tensor_outname),
    // PR job: input is the reloaded post-QL tarball (see clus_pr above).
    pr(anodes, dump=true, pipeline_names=[], tensor_outname='',
       trackfitting_config_file='', particle_dataset=null, extra_uses=[],
       dl_weights='', beam_window=[0, 0], tgm_neutrino_candidate=false,
       tgm_chord_charge=false, tgm_chord_mode='chord',
       tgm_component_extremes=false)::
        clus_pr(anodes, dump=dump,
                output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                pipeline_names=pipeline_names, tensor_outname=tensor_outname,
                trackfitting_config_file=trackfitting_config_file,
                particle_dataset=particle_dataset, extra_uses=extra_uses,
                dl_weights=dl_weights, beam_window=beam_window,
                tgm_neutrino_candidate=tgm_neutrino_candidate,
                tgm_chord_charge=tgm_chord_charge,
                tgm_chord_mode=tgm_chord_mode,
                tgm_component_extremes=tgm_component_extremes),
    detector_volumes(anodes, face=0):: detector_volumes(anodes=anodes, face=face, pos_offset_on=pos_offset_on),
}
