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

// SBND top face of the active volume, cm.  Used to re-anchor the PR chain's
// uBooNE-calibrated cosmic_tagger y cuts (docs/pr/2 sec 2e(iv)): the uBooNE
// prototype was written against y in [-116,+117] cm, z in [0,1037] cm
// (prototype_base/pid/apps/wire-cell-prod-nue.cxx:417), while SBND's sensitive
// box is |y| <= 199.965 cm, z in [0,501.0] cm (see the DetectorVolumes note in
// clus_pr below).  Change THIS number, not the per-cut arithmetic in the
// clus_pr defaults below.
local sbnd_y_top = 200.0;

// DIAGNOSTIC (default OFF): dump the Bee "clustering" layer once per pipeline
// step, so the step that merged two pieces into one cluster can be named instead
// of guessed.  MABC's bee_points_sets entries accept a "visitor" = the step's
// type:name, and MultiAlgBlobClustering.cxx:2270 dumps that set right after that
// step runs (AFTER the per-step enumerate_idents renumbering, so match pieces by
// point coordinates, never by cluster id -- ids shift between layers).
// Layer names are tr<NN>_<Type>; the index disambiguates a type used twice
// (e.g. the two ClusteringRegular passes).
// Off => bee_points_sets list unchanged => compiled config byte-identical.
local trace_sets(pipeline, coords) = [
    {
        name: 'tr%02d_%s' % [i, std.split(wc.tn(pipeline[i]), ':')[0]],
        detector: 'sbnd',
        algorithm: 'clustering',
        pcname: '3d',
        coords: coords,
        individual: false,
        visitor: wc.tn(pipeline[i]),
    }
    for i in std.range(0, std.length(pipeline) - 1)
];

// SBND Space-Charge-Effect displacement field (per-TPC TH3 maps).  Wired into
// DetectorVolumes per-APA metadata (key "sce_field"); PCTransformSet picks it up
// by TypeName and builds a real (non-no-op) SCECorrection.  Used for sim
// (use_sce=true) so the all-APA clustering + Bee run in SCE true space (x_sce).
local sce_field = {
    type: 'SCEFieldTH3',
    name: 'sbnd_dualmap',
    data: {
        sce_map_file: '/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_42_00/SCEoffsets/SCEoffsets_SBND_E500_dualmap_CV_voxelTH3.root',
        // sign=+1: the map holds the TrueBkwd (reco->true) offset, so
        // SCECorrection::forward (x_t0 + displacement) gives the true position.
        sign: 1,
    },
};
// FORWARD (true->reco) SCE displacement from the same dualmap file: used by the
// truth labeler to shift the priorSCE (true position) SimEnergyDeposits onto the
// spatially-distorted charge the blobs were reconstructed from.
local sce_field_fwd = {
    type: 'SCEFieldTH3',
    name: 'sbnd_dualmap_fwd',
    data: sce_field.data {
        th3_name_E:   'TrueFwd_Displacement_X_E',
        th3_name_W:   'TrueFwd_Displacement_X_W',
        th3_name_E_y: 'TrueFwd_Displacement_Y_E',
        th3_name_W_y: 'TrueFwd_Displacement_Y_W',
        th3_name_E_z: 'TrueFwd_Displacement_Z_E',
        th3_name_W_z: 'TrueFwd_Displacement_Z_W',
        sign: 1,
    },
};

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

local common_t0cor_coords(pos_offset_on) =
    if pos_offset_on then ['x_t0cor', 'y_cor', 'z_cor'] else ['x_t0cor', 'y', 'z'];
// SCE-corrected (true) scope, produced by a backward SCECorrection step.
local common_sce_coords = ['x_sce', 'y_sce', 'z_sce'];
// use_sce=true -> SCE true space; false -> the T0-corrected reco scope above.
local common_corr_coords(pos_offset_on, use_sce=false) =
    if use_sce then common_sce_coords else common_t0cor_coords(pos_offset_on);

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
        // SCE displacement field for this APA (PCTransformSet reads this key to
        // build a real SCECorrection; a1f0pA inherits it via $.a0f0pA below).
        sce_field: wc.tn(sce_field),
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
    uses: [dv, sce_field],
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
// sep_vertex_veto (SBND default TRUE since doc pr/15, owner decision 2026-08-01):
// separate() un-splits a neutrino-vertex "V" whose two dominant pieces both END
// at their mutual closest approach (run 18255 evt 56463: the nu was cut in two
// at its vertex by the top-cosmic angle ladder).  false omits the key => the
// compiled config is byte-identical to before the knob existed.
// nu_iso_band_guard (SBND default TRUE since doc pr/18, owner decision 2026-08-01):
// the neutrino stage may not merge an isochronous band (narrow drift slab,
// large y-z footprint) with a non-band cluster that spans > 20 cm of drift,
// even when they touch (run 18255 evt 10550: separate correctly splits the nu
// candidate off the cosmic band, then neutrino re-merged them at 0.31 cm
// touch).  false omits both keys => byte-identical pre-knob config.
// iso_cathode_guard (SBND default FALSE, doc pr/19 campaign): clustering_isolated
// declines the angle-less 80 cm small->big absorb for a small cluster that
// reaches within 30 cm of the cathode and is farther from every big cluster
// than from the cathode (run 18253 evt 444187: ~560 pts of near-cathode nueCC
// shower fragments absorbed by a cosmic 46-76 cm away; the true parent sat
// 1.9 cm across the cathode, invisible per-APA).  false omits the key =>
// compiled config byte-identical to before the knob existed.
// nu_band_veto (SBND PRODUCTION ON, owner flip 2026-08-12, doc pr/66):
// nu_iso_band_guard above stops the per-APA chain from merging a band with a
// drift-spanning partner, but the SEPARATE all-APA clustering chain
// (clus_all_apa below) has no iso-band guard of its own and can re-merge the
// exact pair per-APA just refused (run 18255 evt 10550: the 1e1p nu candidate
// and a TGM cosmic band are correctly split at img-global, then fused again
// by the time Q/L runs). When true, each per-APA refusal is recorded as a
// per-blob "nu_band_veto_role" provenance array that merge_clusters() (and
// the cathode bundle rescue's candidate selection) honor everywhere
// downstream, including the all-APA chain -- see
// clus/src/ClusteringFuncs.cxx band_veto_forbids(). false (legacy escape via
// SBND_NU_BAND_VETO=0) omits the key => compiled config byte-identical to
// before the knob existed; only meaningful with nu_iso_band_guard on.
local clus_per_face(anode, face, dump, output_dir, runNo, subRunNo, eventNo, bee_sink=null, rse_from_ident=false, pos_offset_on=true, trace_bee=false, save_assoc_id=false, sep_vertex_veto=true, nu_iso_band_guard=true, iso_cathode_guard=false, nu_band_veto=true) = {
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
        // vertex_veto: un-split a neutrino-vertex "V" (both dominant pieces END
        // at their mutual closest approach) that the top-cosmic angle ladders
        // mistook for a cosmic -- doc pr/15, run 18255 evt 56463 nu cut in two
        // at its vertex.  C++ default false; SBND ON via sep_vertex_veto.
        cm.separate(use_ctpc=true, max_hull_points=100000, sbnd_boundary_tag=true,
                    vertex_veto=sep_vertex_veto),
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
        // nu_iso_band_guard: see the clus_per_face header comment (doc pr/18).
        // nu_band_veto: see the clus_per_face header comment (doc pr/66).
        cm.neutrino(protect_iso_band=nu_iso_band_guard,
                    protect_iso_band_xext=(if nu_iso_band_guard then 20 * wc.cm else null),
                    record_band_veto=nu_band_veto),
        // SBND: tighten the isolated small/big length_cut from the 20 cm default
        // to 15 cm so a ~16 cm EM (gamma) blob is no longer auto-classified
        // "small" and absorbed into a nearby long cosmic track by the
        // angle-less 80 cm small->big merge.  See sbnd_xin/docs/
        // overclustering-evt11-gamma.md.  range_cut left at its 150 default.
        // save_assoc_id: record the main + associated partition this pass creates
        // into "assoc_cluster_id"/"assoc_cluster_main" so it can be undone before
        // the taggers (doc 52).  C++ default false; key omitted when off =>
        // byte-identical compiled config.
        // iso_cathode_guard: see the clus_per_face header comment (doc pr/19).
        // C++ default 0 (guard off); 30 cm covers the 444187 fragment family
        // (cathode reach 0.3-24 cm) while leaving delta rays untouched (they
        // hug their track: big-gap < cathode-distance fails the guard test).
        // dis_floor 40 cm: a small cluster whose big absorber is within 40 cm
        // keeps the legacy absorb (nu-vertex fragments, nearby debris); only
        // DISTANT angle-less absorbs are declined (444187 family: 46-76 cm).
        cm.isolated(length_cut=15 * wc.cm, save_assoc_id=save_assoc_id,
                    cathode_guard_xcut=(if iso_cathode_guard then 30 * wc.cm else null),
                    cathode_guard_dis_floor=(if iso_cathode_guard then 40 * wc.cm else null)),
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
            // The isolated grouping's provenance pair is written HERE, at per-APA
            // scope, so it must be homogenized on THIS node's tensor output too:
            // TensorDM concatenates same-named local PCs across cluster nodes and
            // silently drops a key that is absent from the first-seen node (doc
            // 38's gotcha).  clustering_isolated marks only the clusters it
            // merged, so without this the arrays never reach the all-APA stage and
            // the all-APA save-time fill-in makes every blob look like a main --
            // measured: assoc_cluster_main all 1, the un-merge splits nothing.
            // Same flag as ClusteringIsolated's save_assoc_id: the two are only
            // ever wanted together.  C++ default false; key omitted when off.
            [if save_assoc_id then 'save_assoc_cluster_id']: true,
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
            }] + (if trace_bee then trace_sets(cm_pipeline, ['x', 'y', 'z']) else []),
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
// cathode_rescue_on (SBND default TRUE since doc pr/14 §7.4 validation, owner
// decision 2026-08-01): the cathode BUNDLE rescue -- the isolated patch for
// the flash-reco absorbing-window defect (sbnd_xin/docs/pr/14).  false restores
// the pre-pr/14 pipeline (compiled config byte-identical to before the knob
// existed).  Retire the knob (and its pipeline entry below) when the light
// reconstruction is fixed.
// cathode_rescue_unmatched (SBND default TRUE since doc pr/17, needs
// cathode_rescue_on): second rescue pass adopting NON-MATCHED clusters into
// the beam bundle when they geometrically continue a beam-window cluster
// across the cathode (the 56463 veto-ON mode: the whole far half is
// flashless; sbnd_xin/docs/pr/17; fires 1/1000 mcp1k, 0/48 nueCC48).
// false => rescue_unmatched key suppressed => compiled config byte-identical
// to pre-knob.
// adopt_nu_fragments (SBND default FALSE, doc pr/19 campaign): after the
// cathode rescue passes, adopt small flashless near-cathode fragments (the
// population iso_cathode_guard leaves unabsorbed) into a beam-window cluster
// on raw proximity under the beam-T0 hypothesis.  false omits the keys =>
// compiled config byte-identical to before the knob existed.
// use_sce / reality (from master, merged 2026-08-03): use_sce=true runs the
// all-APA clustering + Bee in SCE true space (x_sce) instead of the T0-corrected
// reco scope (x_t0cor).  Both SBND realities currently set use_sce=false (see
// the reco table in the tail function), so this is a no-op for our chain.
local clus_all_apa(anodes, dump, output_dir, runNo, subRunNo, eventNo, bee_sink=null, premerged=false, rse_from_ident=false, pos_offset_on=true, tensor_outname='', save_real_cluster_id=false, trace_bee=false, save_assoc_cluster_id=false, real_cluster_id_global=null, cathode_rescue_on=true, cathode_rescue_unmatched=true, adopt_nu_fragments=false, save_bundle_main_provenance=false, rescue_allow_in_beam_far=true, rescue_geom_first=true, rescue_pierce_test=true, rescue_pierce_cut=null, rescue_dest_beam_for_new=true, rescue_beam_main_only=true, bee_flash_pred_min=null, use_sce=false, reality='data') = {
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
        prefix='all', detector_volumes=dv, pc_transforms=pcts, coords=common_corr_coords(pos_offset_on, use_sce)),
    // Combined (all-APA) clustering runs AFTER QL charge-light matching, so every
    // cluster carries a matched flash time (cluster_t0).  switch_scope applies the
    // per-cluster T0 correction (x_t0cor scope) and drops any stale per-APA
    // "isolated"/"perblob" array (it destroys+recreates every cluster).  The
    // merging steps below set use_flash_t0=true so they only merge clusters
    // coincident in flash time, and examine_bundles collapses each flash group into
    // one cluster carrying a fresh "isolated"/"perblob" array (main sub-component
    // = -1), like clustering_isolated but grouped by flash time instead of geometry.
    local cm_pipeline = [
        // Scope step: use_sce=true -> SCECorrection (x_sce = reco->true, so the whole
        // downstream pipeline + Bee run in SCE true space); false -> T0Correction (x_t0cor).
        if use_sce then cm_old.switch_scope(correction_name='SCECorrection') else cm_old.switch_scope(),
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
    //
    // tip_touch_cut / crosser_pca_angle (doc pr/20 Part II, A1 + A2).  Both are
    // C++ default 0 = OFF; with the keys omitted the compiled config is what
    // SBND shipped before pr/20.  They relax the two direction terms that were
    // rejecting genuine cathode crossers:
    //   tip_touch_cut=4cm -- when the two tips approach within 4 cm, drop the
    //     cc_pca (connection-vector) term and accept on cluster-PCA collinearity
    //     alone.  cc_pca measures the cathode gap, not the track, at that
    //     separation: of the 183 pairs already accepted by the primary Hough
    //     test 123 (67%) violate the 30 deg conn_far_cut bound.  PDHD ships the
    //     same knob at 3 cm.  +10 edges / 500 SBND data events (evt 315497).
    //   crosser_pca_angle=20 -- raise the tt_pca (cluster-PCA agreement) bound
    //     from angle_cut=10 to 20 deg inside the near-cathode branch, for a
    //     gently curving crosser whose two halves' global PCA axes differ.
    //     +3 edges on the same sample (evt 406796 at 18.3 deg).
    // NOT bit-identical for SBND -- this is a behaviour change delivered as
    // config; gates in sbnd_xin/docs/pr/20.
    + (if cathode_connect_on then [cm.cathode_connect(cathode_x_cut=5*wc.cm, drift_cut=8*wc.cm, min_length_short=2*wc.cm, short_dir_len=25*wc.cm, conn_short_cut=30.0, tip_touch_cut=4*wc.cm, crosser_pca_angle=20.0, flash_t0_window=800*wc.ns)] else [])
    // Cathode BUNDLE rescue (default OFF => list unchanged => byte-identical):
    // joins a crosser whose two halves are in DIFFERENT flash bundles because
    // the flash reco's absorbing window hid the true flash on one side --
    // exactly the pairs cathode_connect's flash gate refuses (doc pr/14).
    // After cathode_connect (only acts on cross-bundle leftovers, with each
    // half maximally assembled) and before examine_bundles (the merged crosser
    // must be ONE flash-collapse member so the PR unmerge keeps it whole).
    // Geometry mirrors the cathode_connect operating point above; the beam
    // window mirrors clus_pr's beam_window / qlmatching's beam_pref_tlow/thigh
    // (keep the three in sync).  Asymmetric search window = owner decision,
    // doc pr/14 sec 3 (hand-scan dt0 spans -7.2 .. +12.1 us).
    + (if cathode_rescue_on then [cm.cathode_bundle_rescue(
        beam_window_low=0.2*wc.us, beam_window_high=2.2*wc.us,
        rescue_t0_early=8*wc.us, rescue_t0_late=13*wc.us,
        cathode_x_cut=5*wc.cm, drift_cut=8*wc.cm,
        min_length_short=2*wc.cm, short_dir_len=25*wc.cm, conn_short_cut=30.0,
        // false => key suppressed => byte-identical (doc pr/17); floors stay
        // at the C++ defaults (30 cm / 200 pts).
        rescue_unmatched=cathode_rescue_unmatched,
        // Fragment adoption (doc pr/19).  false => keys suppressed =>
        // byte-identical.  adopt_dis 13 cm covers the observed 444187 family
        // hop spacing (up to 12.1 cm); other floors stay at the C++ defaults
        // (adopt_xcut 30 cm, frag_max_length 60 cm, min 5 pts, beam >= 10 cm).
        adopt_nu_fragments=adopt_nu_fragments,
        adopt_dis=(if adopt_nu_fragments then 13*wc.cm else null),
        // Rounds 2+3 (sbnd_xin/docs/73).  ALL FIVE SBND PRODUCTION ON since
        // 2026-08-17 (owner flip on the sec-12 round-3 validation: 3-sample
        // valfast clean, 3000-event census PR-examined, hand scan "clearly
        // improvements"; the round-2 knobs alone were adverse -- sec 11 --
        // and ship only together with rescue_beam_main_only + the PR-side
        // round-3 fixes).  NOT bit-identical -- a behaviour change delivered
        // as config; false on all knobs omits every key and restores the
        // byte-identical pre-round-2 compiled config.  C++ defaults for the
        // companion cuts are 8 cm (geom_first_dis) / 8 cm (pierce_cut) / 0.8
        // / 8 cm; only pierce_cut is exposed here, because it is the one the
        // docs/73 sec 6 sweep tunes.
        rescue_allow_in_beam_far=rescue_allow_in_beam_far,
        rescue_geom_first=rescue_geom_first,
        rescue_pierce_test=rescue_pierce_test,
        pierce_cut=rescue_pierce_cut,
        rescue_dest_beam_for_new=rescue_dest_beam_for_new,
        // Round 3 (docs/73 sec 12): the beam-side donor must be its bundle's
        // matched main (evt 51128: a 3.8 cm associated fragment displaced the
        // real 57.7 cm main).  false => key suppressed => byte-identical.
        rescue_beam_main_only=rescue_beam_main_only)] else [])
    + [
        // flags_from_longest: the flash-time merge here collapses a bundle's
        // clusters into one; without this the merged cluster inherits its flags
        // from an arbitrary member, so a matched main that absorbs a tiny
        // co-merged fragment loses flag_main_cluster to it (SBND evt284349:
        // the 2173-pt beam track lost it to a 3-pt TPC1 speck, leaving the flag
        // only on its own out-of-volume shard).  The taggers key on that flag.
        // save_bundle_main_provenance (doc pr/20 Part I P1; C++ default false,
        // key omitted when off => byte-identical pre-knob config): also record,
        // per blob, which members were matched bundle MAINS before this merge
        // demoted them.  Only this merge destroys that fact, and only the
        // all-APA instance runs it, so per-APA output cannot move.
        cm.examine_bundles(use_flash_t0=true, flags_from_longest=true,
                           save_bundle_main_provenance=save_bundle_main_provenance),
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
            // Homogenize the perblob PC key set at tensor-save time so the
            // flash-merge per-blob provenance (real_cluster_id /
            // real_cluster_main) survives into the saved pctree tarball --
            // the serializer silently drops heterogeneous keys.  C++ default
            // false.  Key omitted when off => byte-identical pre-fix tarball.
            [if save_real_cluster_id then 'save_real_cluster_id']: true,
            // Same homogenization for the isolated grouping's provenance pair, so
            // assoc_cluster_id/assoc_cluster_main survive into the PR job's
            // pctree tarball.  C++ default false; key omitted when off.
            [if save_assoc_cluster_id then 'save_assoc_cluster_id']: true,
            // Re-stamp real_cluster_id at save time into ONE globally unique
            // ident epoch (doc 53).  Without it the array mixes the numbering
            // examine_bundles recorded with the numbering enumerate_idents has
            // installed since, so 31% of values name two clusters -- fine for
            // the within-cluster consumers, wrong for anything that joins on it.
            // Group membership is unchanged either way, so the un-merge and TGM
            // are verdict-neutral.  Applies to BOTH the saved pctree and this
            // node's Bee per-blob labels: the re-stamp runs right after the
            // clustering pipeline, before the Bee fills.  Per-step trace_bee
            // layers are mid-pipeline snapshots and keep their own ids.
            // TRISTATE: C++ default TRUE (doc 53, owner decision); null here =
            // inherit, key omitted.  Pass false ONLY to reproduce the two-epoch
            // values for A/B archaeology.  Gated by save_real_cluster_id, so it
            // is a structural no-op wherever that is off.
            [if real_cluster_id_global != null then 'real_cluster_id_global']: real_cluster_id_global,
            // Dump the optical flash / charge-light "op" display into this same
            // mabc-all-apa.zip (reads the merged-root "opflash" PC + per-cluster
            // matched-flash association written by QLMatching). bee_detector
            // above supplies the Bee geom.
            save_opflash: true,
            // Emit one "op" row per flash carrying ALL matched cluster ids
            // (op_cluster_ids array) with summed predicted PE, so a flash matched
            // to several clusters shows them together (MicroBooNE-style).
            bee_flash_per_flash: true,
            // doc pr/94 round 3.  Minimum total predicted light (PE) for a
            // matched cluster to appear in op_cluster_ids.  C++ default 100
            // (the legacy dump_light filter): a genuine match below it is drawn
            // as "matched to no flash", which is what made SBND 18255/73038's
            // 26.5 cm beam-flash-matched activity (3.6 PE predicted of the
            // flash's 602.6 PE) look unmatched while the PR chain reconstructed
            // it.  null => key omitted => byte-identical display; 0 shows every
            // genuine match.
            [if bee_flash_pred_min != null then 'bee_flash_pred_min']: bee_flash_pred_min,
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
                    coords: common_corr_coords(pos_offset_on, use_sce),
                    individual: false,
                },
            ] + (if trace_bee then trace_sets(cm_pipeline, common_corr_coords(pos_offset_on, use_sce)) else []),
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
            // Write a real tensor output when a tensor_outname is set; otherwise
            // keep the historical discard.
            dump_mode: tensor_outname == '',
        },
    }, nin=1, nout=0),
    // This maker produces ONLY the clustering + matching all-APA MABC (optionally
    // terminated by its own TensorFileSink when dump=true).  The follow-up PR
    // tagger pass (clus_pr / the maker's pr() method) and the larwirecell
    // wclsTensorSetLabeler are NOT wired here: the entry configuration
    // (e.g. sbnd/wcls-img-clus-matching-xin.jsonnet) assembles
    //   MABC -> pr() -> wclsTensorSetLabeler -> dump
    // itself, so clus.jsonnet stays clustering+matching only.
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
              // trackfitting_config_file: the SBND TrackFitting parameter JSON.
              // DEFAULT = the canonical in-tree file, resolved through
              // WIRECELL_PATH by TaggerCheckSTM/TaggerCheckNeutrino
              // (Persist::resolve; an absolute path is passed through unchanged).
              // '' selects the C++ preset defaults, which are uBooNE-hard-coded --
              // never right for SBND, which is why this no longer defaults to ''.
              trackfitting_config_file='pgrapher/experiment/sbnd/sbnd_track_fitting.json',
              particle_dataset=null, extra_uses=[],
              // dl_weights: SCN (DL) neutrino-vertex weights, WIRECELL_PATH-resolved.
              // DEFAULT = the uBooNE-trained net, i.e. the DL vertex is ON for SBND
              // (owner adopted 2026-07-30 on nueCC48 evt 18253/1/172230, where the
              // geometric vertex sat at the far end of a proton track and the DL
              // vertex moved it 9.7 cm onto the true interaction point -- docs/pr/4).
              // '' restores the geometric-only vertex (still what every identity
              // gate must pass, CLAUDE.md M4).  The weights remain uBooNE-TRAINED:
              // SBND retraining is docs/pr/2 gap G3 and stays open.
              // REQUIRES libpython preloaded in the job (LD_PRELOAD=libpython3.11.so.1.0);
              // without it the SCN module import fails and the code falls back to the
              // geometric vertex with a single WARN line -- grep the log for
              // "DL vertex failed" (expect none).
              dl_weights='uboone/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth',
              // DL re-rank operating point, threaded for A/B runs (doc pr/79).
              // min_accept: rerank-total acceptance threshold for the DL route
              // (TaggerCheckNeutrino); top_k: DL voxel candidates admitted to
              // the reranker.
              // min_accept 4.0 -> 10.0 adopted 2026-08-15 (owner, doc pr/79):
              // live A/B on the 473-label hand-scan measured +36/473 correct
              // (322 -> 358; 51 fixed / 15 regressed, every regression's
              // baseline DL total in [4,10)).  Pass 4.0 for the pre-flip arm.
              dl_vtx_min_accept_score=10.0,
              dl_vtx_top_k=5,
              // doc pr/105: dl_vtx_rerank=false selects the legacy single-argmax
              // DL branch (top voxel snapped, dl_vtx_cut gate, else the
              // traditional vertex) -- strategy 3.1 of the vertex-selection
              // comparison.  Default true = the pinned production value, so the
              // compiled JSON is byte-identical when unset.
              dl_vtx_rerank=true,
              // beam_window: internal-unit [low, high] on the matched flash time
              // (cluster_t0).  DEFAULT = the SBND BNB gate after the
              // frame_apply_at_caf correction, which is what production passes
              // (run_nusel_evt.sh BEAM_WINDOW="0.2,2.2").  [0,0] disables it and
              // makes beam_window_only inert.
              beam_window=[0.2 * wc.us, 2.2 * wc.us],
              // --- TGM/FC knobs.  Every default below is the SBND PRODUCTION
              // operating point (run_full1k_nusel.sh's flag set: -chord -rescue
              // -rescue-chord -fvz 5 -fvzi 3 -main-pair-real -fvx 2.5 -fvy 3,
              // plus NUCAND/CHORD_MODE defaults), adopted as the module defaults
              // 2026-07-27 (owner; sbnd_xin/docs/64_cfg-sync.md).  They used to
              // sit at the pre-adoption values here and were reached only via
              // runner flags, so cfg/ advertised a configuration nobody ran.
              // Each C++ knob still defaults off/legacy, so other detectors are
              // unaffected; pass the old value here for a pre-adoption A/B:
              //   tgm_neutrino_candidate false, tgm_chord_charge false,
              //   tgm_chord_mode 'chord', tgm_component_extremes false,
              //   tgm_component_rescue false, tgm_rescue_chord false,
              //   tgm_main_pair false, tgm_main_pair_mode 'path',
              //   tgm_fv_zmax_margin 3, tgm_fv_zmax_margin_interior 0,
              //   tgm_fv_x_margin 2, tgm_fv_y_margin 2.5. ---
              tgm_neutrino_candidate=true, tgm_chord_charge=true,
              tgm_chord_mode='path', tgm_component_extremes=true,
              tgm_component_rescue=true, tgm_rescue_chord=true,
              tgm_main_pair=true, tgm_main_pair_mode='real',
              tgm_fv_zmax_margin=5, tgm_fv_zmax_margin_interior=3,
              tgm_fv_x_margin=2.5, tgm_fv_y_margin=3,
              save_stm_fit=false, unmerge_bundle_mode='real',
              // doc pr/34 §10 particle-flow (Bee mc tree) port-fidelity knobs.
              // C++ defaults false; keys omitted when off => byte-identical
              // pre-knob config.  Display-only stage: mc.json is the artifact.
              pf_track_main_cluster_only=false,
              pf_track_bridged_clusters=false,   // doc pr/40 round 9 B2; C++ default false. Key omitted when off => byte-identical.
              pf_shower_vertex_barrier=false,
              pf_shower_parent_precedence=false,
              pf_pi0_node_per_id=false,
              pf_pdg_name_prototype_fallback=false,
              // doc pr/38 Round 4: graph-faithful parentage for barrier-orphaned
              // PF track segments.  C++ default false; key omitted when off =>
              // byte-identical pre-knob config.  Display-only (mc.json).
              pf_orphan_track_parentage=false,
              // doc pr/65 round 3: audit-only orphan net -- log unclaimed
              // segments, fabricate no PF-root node.  C++ default false; key
              // omitted when off => byte-identical.  Display-only (mc.json).
              pf_orphan_audit_only=false,
              // doc pr/84 round 2 (F1/F2): vertex-touching pseudo-parent
              // suppression + remote-gap anchor.  C++ defaults false/3cm/8cm;
              // keys omitted when off/null => byte-identical pre-knob config.
              // Display-only (mc.json).  Scalar params are in cm.
              pf_direct_when_touching=false,
              pf_touch_max=null,
              pf_touch_cross_main=false,
              pf_touch_cross_max=null,
              pf_pseudo_gap_from_main=false,
              // doc pr/84 round 3 (G1): guarantee unique jsTree node ids in
              // mc.json.  C++ default false; key omitted when off =>
              // byte-identical pre-fix config.
              pf_unique_node_ids=false,
              // pf_drop_stray_satellites (doc pr/92): skip showers the kine
              // side's kine_drop_stray_satellites gate dropped from Enu, so
              // the PF tree and Enu describe the same particle set.  Inert
              // unless that kine knob is also on.  C++ default false; key
              // omitted when off => byte-identical pre-pr/92 config.
              pf_drop_stray_satellites=false,
              // pf_orphan_confident_track (doc pr/93 round 4): emit a root PF
              // node for an unclaimed main-cluster segment with a confident
              // non-electron template PID, length > pf_orphan_track_min_cm
              // (null => C++ 50cm), and straight-long (SBND 18255-315167's
              // freed 150.7cm proton).  C++ default false; keys omitted when
              // off => byte-identical.
              pf_orphan_confident_track=false,
              pf_orphan_track_min_cm=null,
              // pf_track_owns_loose_vertex (doc pr/93 round 4): a vertex the
              // track BFS reached via a real segment is not claimable by a
              // root shower whose only tie to it is the loose fill_sets view
              // (SBND 18264-69314's 67 MeV daughter stolen from its muon).
              // C++ default false; key omitted when off => byte-identical.
              pf_track_owns_loose_vertex=false,
              // restore_demoted_mains (doc pr/20 Part I P2; C++ default false,
              // key omitted when null => byte-identical pre-knob config): tag a
              // split-off part that was ITSELF a matched bundle main before the
              // flash-time merge with flag_demoted_main.  Requires the Q/L stage
              // to have run with save_bundle_main_provenance, else the visitor
              // warns and flags nothing.  Only the OUTER (flash) un-merge takes
              // it -- unmerge_assoc undoes the isolated grouping, a different
              // question.
              restore_demoted_mains=null,
              // require_provenance (doc pr/23 sec 4.2; C++ default false, key
              // omitted when null => byte-identical pre-knob config): with
              // restore_demoted_mains on, ABORT on a pctree with no wasmain
              // array (a legacy pre-pr/20 tree) instead of warn-and-skip.
              // Legacy-tree runs opt out with false explicitly.
              require_provenance=null,
              // evaluate_demoted_mains (doc pr/20 Part I P3; C++ default false,
              // key omitted when false => byte-identical pre-knob config): let
              // TaggerCheckTGM / STM / FC evaluate a cluster carrying
              // flag_demoted_main.  Inert unless restore_demoted_mains above put
              // the flag there.
              evaluate_demoted_mains=false,
              // tgm_exempt_demoted_main (doc pr/25, SBND evt 320029; C++ default
              // false, key omitted when false => byte-identical pre-knob
              // config): with tgm_main_pair on, TaggerCheckTGM's
              // main_pair_rejects vetoes EVERY demoted-main pair
              // unconditionally (its own real_cluster_main/path-component
              // provenance is all-zero by construction post-carve) -- this
              // exempts flag_demoted_main clusters from that veto so the
              // usual CASE-A/CASE-B boundary geometry actually runs on them.
              // DESIGNED, NOT the SBND default -- changes cosmic verdicts,
              // owner sign-off pending (escalation rule 1).
              tgm_exempt_demoted_main=false,
              // skip_cosmic_companions / cosmic_companion_min_length (doc pr/20
              // Part I P4; C++ defaults false / 0, keys omitted when off =>
              // byte-identical): act on the verdict P3 produced -- drop a TGM/
              // STM-tagged companion from the neutrino's other_clusters, unless
              // it is shorter than the floor (cm).  Inert without P3.
              skip_cosmic_companions=false, cosmic_companion_min_length=null,
              // nu_fallback_demoted_mains (docs/73 sec 12, round 3): when NO
              // candidate survives the primary loop, consider demoted mains
              // (same gates).  Inert without restore_demoted_mains upstream;
              // pairs with evaluate_demoted_mains (P3) so the candidates carry
              // tagger verdicts.  SBND PRODUCTION ON since 2026-08-17
              // (docs/73 sec 12 owner flip).
              nu_fallback_demoted_mains=true,
              // sp_photon_flag: store the single-photon tagger's verdict in
              // TaggerInfo::photon_flag, as prototype NeutrinoID.cxx:271 does.
              // The port ran singlephoton_tagger() and filled its shw_sp_*
              // features but discarded the verdict (docs/pr/26 sec. 8.2).
              // C++ default false = that gap; key omitted when off =>
              // byte-identical.  Only the uBooNE tagger ntuple's photon_flag
              // branch changes when on -- nothing in the chain reads it.
              sp_photon_flag=false,
              // mip_dqdx: SBND MIP dQ/dx scale in e/cm handed to
              // TaggerCheckSTM AND (since docs pr/7-pr/8) to
              // tagger_check_neutrino as the PR chain's flat-template/cal_4mom
              // amplitude (C++ default 50000 = MicroBooNE).  56000 is the
              // SBND value derived in sbnd_xin/docs/48; pass 50000 to isolate
              // the reference-table change from the MIP-scale change in an A/B.
              mip_dqdx=56000,
              // stm_consistent_fv: TaggerCheckSTM's containment gate uses the
              // same fiducial + margins as tagger_check_tgm / tagger_check_fc.
              // Default TRUE = SBND production; pass false to restore the
              // pre-doc-49 FiducialUtils fallback for an A/B.  Details at the
              // tagger_check_stm call below.
              stm_consistent_fv=true,
              // stm_accept_guards: the doc-63 round-1 STM acceptance guards
              // (charge-desert one-objectness veto, spike-not-ramp nu-vertex
              // veto, eval ratio2 normalization cap).  C++ default false.
              // DEFAULT TRUE = SBND production as of doc 63 (owner
              // 2026-07-26); pass false for a pre-campaign A/B (key then
              // omitted => byte-identical pre-knob config).
              stm_accept_guards=true,
              // stm_proton_muon_guard: the doc-63 round-2 muon-consistency
              // guard on detect_proton (C++ default false; key omitted when
              // off => byte-identical).  DEFAULT TRUE = SBND production as
              // of doc 63 (owner 2026-07-26).
              stm_proton_muon_guard=true,
              // stm_cathode_guard: the doc-63 round-3 cathode-truncation veto
              // (C++ default false; key omitted when off => byte-identical).
              // DEFAULT TRUE = SBND production as of doc 63 (owner
              // 2026-07-26).
              stm_cathode_guard=true,
              // stm_anode_dist_fix: the doc-63 round-4a fix of the inverted
              // face selection in check_stm_conditions' dist_to_anode helper
              // (the anode-clipped-TGM catch fired at the SBND cathode).
              // C++ default false; key omitted when off => byte-identical.
              // DEFAULT TRUE = SBND production as of doc 63 (owner
              // 2026-07-26, after validation).
              stm_anode_dist_fix=true,
              // stm_second_track_guard: the doc-63 round-4b/4c second-track
              // vetoes (long straight leftover past the kink; long MIP-like
              // other-track segment).  C++ default false; key omitted when
              // off => byte-identical.  DEFAULT TRUE = SBND production as of
              // doc 63 (owner 2026-07-26, after validation).
              stm_second_track_guard=true,
              // stm_deficit_guard / stm_vertex_kink_guard: the doc-63
              // round-5 stop-region vetoes (charge-deficient end with no
              // Bragg = truncation; sharp end turn into a hot prong =
              // vertex).  C++ defaults false; keys omitted when off =>
              // byte-identical.  DEFAULT TRUE = SBND production as of doc 63
              // (owner 2026-07-26, after validation).
              stm_deficit_guard=true,
              stm_vertex_kink_guard=true,
              // stm_d66_cuts: the doc-66 sec 12 diffusion-margin cut package
              // (Michel-veto res_length floor 6 -> 6.5 cm, detect_proton
              // track_medium gate 1.0 -> 1.05, block-B ks2 entry 0.05 ->
              // 0.055, C1 peak clause 4.3 -> 4.1).  The 4.0/8.8 diffusion
              // revert moved four owner-adjudicated bundles across the
              // prototype constants by hairline margins; these values were
              // chosen from a full-1000-event margin sweep with zero measured
              // collateral (sbnd_xin/docs/66 sec 12.2).  C++ defaults are the
              // prototype constants; false omits the keys => byte-identical
              // pre-package config.  DEFAULT TRUE = SBND production as of doc
              // 66 (owner 2026-07-27).
              stm_d66_cuts=true,
              stm_michel_res_cm=6.5,
              stm_proton_tm_max=1.05,
              stm_proton_b_ks2_max=0.055,
              stm_proton_c_peak_max=4.1,
              // beam_window_only: run the PR tail (steiner + TGM + STM + FC)
              // ONLY on the beam-coincident bundle -- the main clusters whose
              // matched flash time (cluster_t0) falls in beam_window, plus (for
              // steiner) the companions sharing their matched_flash_gid.
              // DEFAULT TRUE = SBND production as of doc 56: the taggers are
              // validated cosmic/containment verdicts for the beam bundle, and
              // the out-of-time bundles they used to be run on are cosmics by
              // construction (~11 bundles/event, 1 in window).  Pass false to
              // restore the evaluate-every-bundle behavior (C++ knobs default
              // off, keys omitted => compiled config byte-identical to the
              // pre-doc-56 one).  Inert if beam_window is empty.
              beam_window_only=true,
              // nu_skip_cosmic: tagger_check_neutrino skips in-window mains
              // already tagged cosmic upstream (flag_TGM, flag_STM, or Q/L
              // lm_flag>0), so neutrino PR runs only on untagged nu candidates.
              // Per-main: a cosmic-tagged longest bundle does not veto a clean
              // runner-up.  C++ default false (key omitted when off =>
              // byte-identical uBooNE config).  DEFAULT TRUE for SBND (owner
              // 2026-07-30, sbnd_xin/docs/pr/3); zero production impact while
              // tagger_check_neutrino is not in pipeline_names.
              nu_skip_cosmic=true,
              // nu_skip_cosmic_bundle: lift that verdict from the main to the
              // whole flash bundle -- if ANY in-window main sharing a
              // matched_flash_gid is cosmic-tagged, no main of that bundle is
              // eligible, so neutrino PR does not run on the bundle at all.
              // Needed because a bundle can hold a second main: SBND evt
              // 18255/52195 gid 6 holds the TGM-tagged 513 cm cathode-crosser
              // AND a 5-point 1.7 cm shard, and under the per-main rule the
              // shard became the nu candidate and pulled the bundle's untagged
              // 400 cm associated muon through a full PR + track fit.
              // C++ default false (key omitted when off => byte-identical
              // uBooNE config).  DEFAULT TRUE for SBND (owner 2026-08-01,
              // sbnd_xin/docs/pr/3 sec. 8).  NOT bit-identical: it removes PR
              // output on cosmic bundles that previously produced it.
              nu_skip_cosmic_bundle=true,
              // nu_skip_cosmic_bundle_min_length (cm): design-A guard on the
              // bundle veto (docs/pr/16 sec. 7).  An UNTAGGED in-window main at
              // least this long survives the veto: the TGM/STM taggers examined
              // it and declined to tag it, and the cosmic-tagged bundle-mate
              // stays out of the PR ensemble regardless (companions are
              // associated-only).  Unexamined shards (out-of-scope mains carry
              // no verdict; evt 52195's 1.3 cm shard) stay vetoed.  C++ default
              // 0 = veto every bundle-mate (key omitted => byte-identical).
              // DEFAULT 15 for SBND (owner 2026-08-01): restores evt 10550's
              // 18.5 cm candidate, keeps 13 cm and below vetoed.
              nu_skip_cosmic_bundle_min_length=15,
              // dir_weak_use_score: route the PR chain's direction-weakness
              // reads through segment_is_dir_weak() -- the faithful port of the
              // prototype's ONLY public accessor ProtoSegment::is_dir_weak()
              // (score>0.07/0.13 thresholds; sentinel score 100 = weak).  The
              // original port read the raw flag at all ~83 sites (docs/pr/6).
              // C++ default false (key omitted => byte-identical uBooNE
              // config).  DEFAULT TRUE for SBND (owner 2026-07-30); zero
              // production impact while tagger_check_neutrino is not in
              // pipeline_names.
              dir_weak_use_score=true,
              // mip_dqdx_median: the PR chain's median-dQ/dx threshold scale in
              // e/cm (C++ default 43000 = uBooNE; every x1.2/1.3/1.4/1.75...
              // ratio cut is quoted against it -- docs pr/7 sec 5, pr/8).  The
              // 50k-role flat-template amplitude reuses the existing mip_dqdx
              // arg above (56000).  48000 keeps the uBooNE 43/50 internal ratio
              // (owner 2026-07-30, placeholder pending an SBND median-MIP
              // measurement).  null falls back to the uBooNE default.
              mip_dqdx_median=48000,
              // Proton-template direction vote (doc pr/8): when the muon-vs-flat
              // direction gate abstains both ways, a proton-consistent, clearly
              // asymmetric template fit toward a free end declares the direction
              // (guards G1-G4).  C++ default false (key omitted => byte-identical
              // uBooNE config).  DEFAULT TRUE for SBND (owner 2026-07-30);
              // score/asym thresholds stay at the C++ defaults 0.25/1.3 pending
              // the pr/8 sec 6 calibration.
              proton_dir_vote=true,
              // Endpoint-trim retry (doc pr/9 sec 6 F1): on double abstention,
              // retry the dQ/dx PID once with exactly 1 sample excluded at the
              // hypothesized stopping end -- an ill-defined endpoint dilutes or
              // inflates the tip dQ/dx, which is compared against the template's
              // Bragg maximum (evt 172230: tip -37% vetoed a clean proton;
              // trimmed retry passes at score 0.077).  Dynamic: trims 0 samples
              // when the untrimmed comparison already decides, never more than
              // 1.  Runs before the proton_dir_vote fallback.  C++ default
              // false (key omitted => byte-identical uBooNE config).  DEFAULT
              // TRUE for SBND (owner 2026-07-30).
              endpoint_trim_retry=true,
              // fit_vertex short-segment exclusion (doc pr/9 sec 11 F3c), cm.
              // Segments with wcpt-path length below the cut do not count as
              // vertex-fit legs (they stay in the graph/particle flow): a
              // 0.62 cm drift-blur vertex-activity stub dragged the evt 172230
              // main vertex 2.4 cm via the post-stub refit.  >=3 surviving legs
              // => fit on the survivors; <=2 => skip the fit (the two-leg
              // position was already fit).  null = C++ default 0 = legacy
              // include-all (key omitted => byte-identical uBooNE config).
              // DEFAULT 1.0 cm for SBND (owner 2026-07-30, deliberate
              // prototype divergence).
              fit_vertex_min_seg_length=1.0,
              // cathode_x / cathode_kink_xcut (cm): the cathode kink veto,
              // doc pr/20 Part II B0.  segment_search_kink's four accept
              // criteria are each guarded by para_angle > 7.5-15 deg, a guard
              // against isochronous imaging artifacts that is INVERTED at the
              // cathode -- the crossing is drift-x dominated (para_angle
              // 61-78 deg, wide open) while the ~2 cm transverse mismatch
              // across the gap supplies the turn, so one cosmic leaves the PR
              // graph as two segments.  cathode_kink_xcut suppresses kink
              // ACCEPT candidates within that distance of x = cathode_x; the
              // angle arithmetic itself is untouched.
              // null = C++ default 0 = OFF (key omitted => byte-identical).
              cathode_x=null, cathode_kink_xcut=null,
              // cathode_wide_kink_angle (deg) / _skirt / _baseline (cm): the
              // wide-baseline cathode kink ACCEPT, doc pr/47 sec 8 (O1) --
              // the converse of the veto above.  At a cathode-crossing fit
              // index the shipped index-windowed refl_angle statistic is
              // suppressed by the gap/distortion (52085: a genuine 33-38 deg
              // two-prong junction measures only ~23 deg locally), so a fifth
              // accept path fires when the skirt-excluded PCA turn angle
              // across the crossing exceeds the angle cut.  null = C++
              // default 0 = OFF (key omitted => byte-identical).
              cathode_wide_kink_angle=null,
              cathode_wide_kink_skirt=null,
              cathode_wide_kink_baseline=null,
              // shower_topo_demote_len (cm, doc pr/25 sec 3): demote any
              // kShowerTopology segment longer than this to a track.  Owner
              // hand-scan 2026-08-03: 10/10 long shower-topology segments on
              // a selected nu-candidate main cluster are tracks, none showers.
              // null = C++ default 0 = OFF (key omitted => byte-identical).
              shower_topo_demote_len=null,
              // ---- doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs ------------
              // All five default to the pre-pr/30 behaviour, so the compiled JSON is
              // byte-identical until one is set.  See cfg/pgrapher/common/clus.jsonnet
              // for what each does.  fit_exclusion (P1), graph_endpoint_strict (P8) and
              // oov_prototype_parity (F2) turn NEW behaviour on (default off);
              // first_seg_local_pca (P2) and other_seg_relaxed_accept (P4) gate behaviour
              // that is ALREADY production, so null = on = legacy and false restores the
              // prototype's narrower version.
              fit_exclusion=false, graph_endpoint_strict=false, graph_endpoint_tol=null,
              oov_prototype_parity=false, first_seg_local_pca=null, other_seg_relaxed_accept=null,
              // shower_topo_proto_dir (doc pr/31 sec 11, F2 was P2): true skips the
              // stage-3 segment_determine_shower_direction call, leaving the topology
              // shower with the direction segment_is_shower_topology set -- the
              // prototype's state.  C++ default false = today's path = byte-identical.
              shower_topo_proto_dir=false,
              // doc pr/32 sec 11: the four stage-4 vertex-ID port fixes.  C++
              // default false = today's path; keys omitted => byte-identical.
              // The SBND operating point lives in wct-pr-perevt.jsonnet.
              vertex_dir_use_fit_point=false, shower_traj_recheck_parity=false,
              main_vertex_require_descriptor=false, main_vertex_candidate_flag=false,
              // doc pr/31 sec 12: the sec 10.12 topology/PID/direction port
              // fixes (F5 cont-muon 30cm dir3, F6 empty-window abstain, F3
              // shower-topo reset, F1 preserve-4mom, F4 local median, F7
              // vertex-by-index -- F7 deliberately dormant pending pr/30 F4).
              // C++ default false = today's path; keys omitted =>
              // byte-identical.  The SBND operating point lives in
              // wct-pr-perevt.jsonnet.
              cont_muon_dir3_30cm=false, track_comp_empty_abstain=false,
              shower_topo_reset=false, reclass_preserve_4mom=false,
              dir_track_median_local=false, examine_showers_vertex_by_index=false,
              // Steiner TERMINAL filter fidelity (doc pr/29 D1 and D12).  Both
              // OFF here = the historical toolkit behaviour, keys omitted =>
              // byte-identical config.  Turning either on can only ADD Steiner
              // terminals, never remove them, because both restore a way for a
              // terminal to PASS the reference-containment test:
              //   steiner_terminal_wire_tol=1     the prototype's one wire of
              //     slack on both sides of all three planes.  Applies ONLY to
              //     the terminal filter; get_extreme_wcps shares the C++ helper
              //     and correctly keeps the exact test the prototype gives it.
              //   steiner_terminal_adjacent_slice=true   makes the t+-1 slice
              //     fallback actually resolve.  The time_blob_map key is in
              //     TICKS, so the historical step of 1 matches nothing on SBND
              //     (4 ticks per slice) and the whole branch is dead.
              // NOT flipped on by default: that is an unconditional production
              // output change (CLAUDE.md sec.5 rule 1) and is the owner's call.
              steiner_terminal_wire_tol=0,
              steiner_terminal_adjacent_slice=false,
              // Steiner EDGE-WEIGHT charge fidelity (doc pr/29 D2).  OFF here =
              // the historical toolkit behaviour, key omitted => byte-identical.
              //   steiner_edge_charge_forward_dead_mix=true   weights steiner
              //     edges with charges computed under the disable_dead_mix_cell
              //     this chain actually passes (false), as the prototype does;
              //     the toolkit dropped the argument at the call and always
              //     took the true branch of calc_charge_wcp.
              // Unlike the two above this can move a weight in EITHER
              // direction, so it can add or drop tree edges.
              steiner_edge_charge_forward_dead_mix=false,
              // Isochronous first-segment endpoint finding (doc pr/24 round 2,
              // SBND evt 271851): principal-axis endpoints for filled 2-D
              // sheet clusters instead of the wire-footprint boundary metric.
              // false/nulls = C++ defaults (false / 40 cm / 25 cm / 0.35 /
              // 0.02) = OFF (keys omitted => byte-identical).
              iso_endpoint=false,
              iso_endpoint_min_length=null,
              iso_endpoint_max_xext=null,
              iso_endpoint_xext_frac=null,
              iso_endpoint_xext_quantile=null,
              iso_endpoint_tube_radius=null,
              iso_endpoint_min_aspect=null,
              // examine_vertices_3 extension-retraction guard (doc pr/24
              // round 5).  false/null = C++ defaults (false / -1.0 cm) = OFF
              // (keys omitted => byte-identical).
              v3_extension_guard=false,
              v3_extension_min_gain=null,
              // doc pr/67: log-only trajectory-coverage probe + the
              // counterfactual override for find_proto_vertex's hardcoded
              // main-cluster branch-search round budget.  false/null =
              // C++ defaults = OFF (keys omitted => byte-identical).
              traj_cover_probe=false,
              pr_find_other_rounds=null,
              // protect_bundle stage knobs (doc pr/23): the PR-stage
              // overclustering protection (uboone's second graph examination,
              // ProtectOverClustering.cxx).  The stage only acts when
              // 'protect_bundle' is named in pipeline_names -- which the
              // production jobs do by DEFAULT since the doc pr/23 sec 9 flip
              // (owner 2026-08-02).  protect_graph_name: connected_blobs flavor
              // (null => 'relaxed').  protect_cathode_*: the SBND cathode
              // re-join divergence -- INTERNAL units (unlike cathode_kink_xcut
              // one block up, which is cm); nulls => C++ defaults, i.e. the
              // re-join pass disabled (prototype-faithful).
              protect_graph_name=null,
              // C++ default true (doc pr/23 ordering): convicted in-window
              // mains (TGM/STM/lm) do not open their bundle.  null => key
              // omitted => C++ default.
              protect_skip_convicted=null,
              // doc pr/94 round 3: let a convicted main open its bundle so the
              // bundle's unconvicted members are graph-examined.  C++ default
              // false; key omitted when null.
              protect_open_convicted_bundles=null,
              protect_cathode_x=null,
              protect_cathode_rejoin_xcut=null,
              protect_cathode_rejoin_dyz=null,
              protect_cathode_rejoin_dis=null,
              // Direction-agreement fallback for a pair that fails ONLY
              // cathode_rejoin_dyz (doc pr/25, SBND evt 489327): dyz is a
              // frame-aligned transverse-offset bound and rejects an OBLIQUE
              // crosser by construction (the offset across the cathode gap is
              // real track travel, not misalignment).  nulls => C++ default
              // 0 = fallback disabled => byte-identical to the block above.
              // INTERNAL units except protect_cathode_rejoin_angle (degrees).
              protect_cathode_rejoin_perp=null,
              protect_cathode_rejoin_angle=null,
              protect_cathode_rejoin_dir_radius=null,
              protect_cathode_rejoin_dir_npts=null,
              // cosmic_y_* (cm): cosmic_tagger()'s four "reaches the top of the
              // detector" tests, re-anchored from uBooNE's top face (y = +117 cm)
              // to SBND's (y = +200 cm, sbnd_y_top above).  The uBooNE literals
              // 100 / 102 / 80 / 50 cm are 17 / 15 / 37 / 67 cm below its top --
              // an entry-point tolerance for a downward cosmic, which does not
              // scale with detector height -- so SBND keeps the same offsets and
              // gets 183 / 185 / 163 / 133 cm.  At the uBooNE values these cuts
              // are near-meaningless on SBND (100 cm is mid-detector), which is
              // what docs/pr/2 sec 2e(iv) flagged.  DEFAULT ON for SBND (owner
              // 2026-07-30); null on any one of them restores that cut's uBooNE
              // value.  NOT bit-identical: cosmic_tagger verdicts can change.
              cosmic_y_top_main=sbnd_y_top - 17,     // main cluster's own top
              cosmic_y_top_strict=sbnd_y_top - 15,   // event top, 1-cosmic branch
              cosmic_y_top_loose=sbnd_y_top - 37,    // event top, global gate
              cosmic_y_small_piece=sbnd_y_top - 67,  // <3 cm debris, PCA centre
              // vertex_z_prior_scale (cm): denominator of the upstream-z penalty
              // (z - min_z)/scale that ranks main-vertex candidates against the
              // +0.25-per-track bonuses.  uBooNE's 200 cm spans ~5.2 penalty
              // units over its 1037 cm detector; SBND is 501 cm long, so the
              // same 200 cm halves the prior's dynamic range -- most visible in
              // compare_main_vertices_global, which ranks candidates from
              // DIFFERENT clusters of the beam bundle, where separations do run
              // toward the full detector length.  DEFAULT = the length-scaled
              // 200 x 501/1037 = 96.6 -> 100 cm (owner 2026-07-30).  The
              // alternative reading -- that the prior is a per-cm trade-off
              // against track bonuses and transfers unchanged -- keeps 200;
              // pass null for that (docs/pr/2 sec 2e(iv)).
              // NOT bit-identical: vertex ranking can change.
              vertex_z_prior_scale=100.0,
              // ssm_target_dir / ssm_absorber_dir: the SSM tagger's beam-line
              // reference directions [x,y,z] in the detector frame (docs/pr/2
              // sec 2e(i)).  null = the C++ defaults = the prototype's uBooNE
              // BNB-target / NuMI-absorber vectors, so the compiled config is
              // unchanged.  SBND HAS NO VALUE FOR EITHER YET: the BNB-target
              // direction must be re-derived in the SBND frame, and the
              // NuMI-absorber features have no obvious SBND meaning.  Until
              // then the 8 ssm_*_angle_{target,absorber} features carry uBooNE
              // geometry -- they are only reachable now, not fixed.
              ssm_target_dir=null,
              ssm_absorber_dir=null,
              // kine_*: the charge -> kinetic-energy calibration constants of
              // NeutrinoEnergyReco (docs/pr/2 sec 2e(iii)).  null = the C++
              // defaults = the uBooNE-tuned literals they replaced.
              // The three RECOMBINATION survival factors carry SBND values
              // since 2026-07-30 (docs/pr/10 sec 6): the uBooNE empirical
              // 0.7 / 0.5 / 0.35 scaled by the table-integrated survival
              // ratio between the official uBooNE Box at 0.273 kV/cm and the
              // doc-55 free-power SBND fit (C excluded -- degenerate with
              // the fudge factors, which deliberately STAY uBooNE: they
              // absorb gain/lifetime normalization, not field physics).
              // NOT bit-identical: kine-derived BDT features move; the
              // T_kine tree of the pr/3-style tracking-pr.root shows the
              // effect (Enu -12..-14% on nuecc48 172230/235435/444187).
              // null on any one restores that factor's uBooNE value.
              // The plane weights [0.25,0.25,1.0], the 0.04 asymmetry
              // switch and kine_w_value (23.6 eV) still have NO SBND value.
              kine_fudge_factor=null,
              kine_recom_factor=0.87,         // 0.70 x 1.249 (track, docs/pr/10)
              kine_shower_fudge_factor=null,
              kine_shower_recom_factor=0.58,  // 0.50 x 1.169 (shower)
              kine_proton_recom_factor=0.51,  // 0.35 x 1.453 (proton)
              kine_plane_weights=null,
              kine_plane_asym_switch=null,
              kine_w_value=null,
              // doc pr/35 sec 10.2 (F1 = P1+P8): live start-segment PDG at the
              // four fill_kine_tree sites (prototype parity).  C++ default
              // false; key omitted when off => byte-identical pre-knob config.
              kine_shower_pdg_live=false,
              // ---- doc sbnd_xin/docs/pr/36 sec 10 tagger-stage knobs -------
              // F1 (= P1): give the match_isFC recompute the SAME fiducial +
              // margins tagger_check_{stm,tgm,fc} use (sbnd_pr_fv +
              // sbnd_pr_fv_margins), mirroring stm_consistent_fv above.
              // false => keys omitted => the historical FiducialUtils
              // fallback, byte-identical pre-knob config.
              neutrino_consistent_fv=false,
              // sbnd_xin/docs/74 G1/G2: give cosmic_tagger()'s containment
              // tests (its inside_fv lambda + the flag-1 vertex test) the
              // SAME sbnd_pr_fv + sbnd_pr_fv_margins as TGM/STM/FC, instead
              // of the FiducialUtils zero-margin sensitive union (which has
              // no wall inset and excludes the CPA slab |x| < 0.45 cm).
              // C++ default false; key omitted when off => byte-identical.
              cosmic_consistent_fv=false,
              // sbnd_xin/docs/75: give the nue/single-photon taggers'
              // containment tests (angular_cut, shower_to_wall,
              // bad_reconstruction_2/_2_sp) the SAME sbnd_pr_fv +
              // sbnd_pr_fv_margins as cosmic_consistent_fv above -- same
              // zero-margin FiducialUtils inconsistency, same fix pattern.
              // C++ default false; key omitted when off => byte-identical.
              nue_sp_consistent_fv=false,
              // F3 (= P2): single-photon SCE correction gate.  Vacuous on
              // SBND today (clus_geom_helper is ''); kept OFF by owner
              // decision 2026-08-04 so a future SBND SCE helper enables it
              // as its own explicit step.  C++ default false.
              sp_sce_correction=false,
              // F4 (= P3+P5): graph-index-ordered tagger accumulation sets
              // (M4 house-rule determinism fix).  C++ default false.
              tagger_ordered_segment_sets=false,
              // F5 (= P6): prototype wcpt-identity stem-endpoint rule at the
              // 18 seg_endpoint_near sites.  C++ default false.
              stem_endpoint_wcpt_parity=false,
              // F6 (= P8): broken_muon_id counts distinct cluster ids.  C++
              // default false.
              broken_muon_cluster_id_count=false,
              // F7 (= P4): neutrino_type verdict bitmask + its T_tagger
              // branch (threaded to BOTH tagger_check_neutrino and
              // tagger_output).  C++ default false.
              neutrino_type_bitmask=false,
              // doc pr/94 Phase 1: per-bundle identity + per-activity
              // cosmic-flag T_tagger branches (tagger_output only, plumbing
              // only -- nothing populates them yet).  C++ default false.
              nu_per_bundle=false,
              nu_per_bundle_min_length=null,   // doc pr/94 Phase 5b; cm; null => C++ default 0 (no floor)
              // doc pr/94 round 3: the selected candidate gets the main-cluster
              // PR treatment for its own pass even when it is a demoted main.
              // C++ default false; key omitted when off.  NOT gated on
              // nu_per_bundle -- the legacy demoted-main fallback needs it too.
              nu_selected_as_main=false,
              // doc 75: closes the DL-swap flag leak nu_selected_as_main's
              // own guard leaves open (see common/clus.jsonnet comment).
              // C++ default false; key omitted when off.
              nu_selected_as_main_snapshot_all=false,
              // ---- doc sbnd_xin/docs/pr/33 §10 EM-shower-clustering knobs.
              // All C++ default false = keys omitted = byte-identical
              // pre-knob config.
              // F1 (= P1): prototype calculate_num_daughter_tracks callee at
              // the main-vertex proton-skip site / the examine_showers
              // daughter_length site.
              daughter_count_proto_main_vertex=false,
              daughter_count_proto_examine_showers=false,
              // F2 (= P2): read the PDG off the object the prototype reads
              // (4 sites start-segment; 1 inverted site shower-type; 2 sites
              // exact ==13 muon test).  Parity at the :170 site needs
              // from_start_segment AND exact_muon_test together.
              shower_pdg_from_start_segment=false,
              shower_pdg_from_shower_type=false,
              shower_pdg_exact_muon_test=false,
              // F3 (= P3): shared pi0-id allocation stream across the two
              // pi0 finders (prevents pio_id collision in the nue tagger
              // pi0 block and the Bee mc.json grouping).
              pi0_id_shared_allocator=false,
              // F4 (= P6): is_shower gains the prototype's abs(pdg)==11
              // disjunct at the cluster-center-point site.
              shower_flag_pdg_electron=false,
              // F5 (= P12): shower_less same-index tie-break by stable
              // shower id (house-rule determinism fix, prototype n/a).
              shower_less_id_tiebreak=false,
              // doc pr/39: exclude a shower's own start vertex from the
              // end_point farthest-vertex search (prototype map_vtx_segs
              // parity).  Ships OFF pending owner gate review.
              shower_endpoint_exclude_start_vertex=false,
              shower_endpoint_skip_orphan_vtx=false,
              // doc pr/91 round 3: flood-fill frontier = visited, not merely
              // present in the view.  Ships OFF pending owner gate review.
              shower_walk_visited_parity=false,
              // doc sbnd_xin/docs/pr/40 -- track (proton/pion/muon)
              // mis-identified as electron.  F1 (persistence), F2/F3 (dQ/dx
              // guards on wholesale track-to-electron conversion).  All
              // default false = legacy = byte-identical.
              track_pid_persist_dqdx=false,
              shower_reclass_dqdx_guard=false,
              shower_topo_dqdx_guard=false,
              // doc sbnd_xin/docs/pr/40 round 2 -- two follow-on defects from the
              // pr/40 round: F1 zero-KE persistence stub, F2 proton-
              // daughter -> pion, F3 reclass_pinfo negative-KE stub.  All
              // default false = legacy = byte-identical.
              reclass_never_computed_ke_floor=false,
              track_pid_persist_4mom=false,
              shower_proton_daughter_pion=false,
              // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on defects from
              // round 2/3's F5: F7 clears the shower flags a relabelled pion
              // still carried (it was still being wrapped as a Shower); F8
              // relabels a muon segment at a multi-proton (>=2, charge-
              // confirmed) non-neutrino-vertex hadronic vertex to pion.  Both
              // default false = legacy = byte-identical.
              shower_proton_daughter_pion_dissolve=false,
              muon_multi_proton_pion=false,
              // doc sbnd_xin/docs/pr/40 round 5 -- muon mis-identified as
              // electron, three independent mechanisms.  F9 narrows F1 so it
              // no longer rescues an undirected electron guess; F10 excludes
              // a long, straight candidate from the main-vertex EM-shower
              // selection; F11 gives segment_is_shower_trajectory the same
              // straightness veto F3 gave segment_is_shower_topology's
              // dQ/dx.  All three default false = legacy = byte-identical.
              track_pid_persist_dqdx_electron_guard=false,
              shower_connect_main_vertex_straight_guard=false,
              shower_traj_straight_guard=false,
              // doc sbnd_xin/docs/pr/40 round 6 -- boundary-level fixes the
              // round-5 measurement demanded.  F12 keeps the shower flood-
              // fill from absorbing a confident straight non-electron track;
              // F13 keeps connecting_to_main_vertex from force-setting a
              // proton-daughter pion to electron; F14 widens the Michel
              // stopping-muon rescue past its weak-dir degree-2 limits.
              // All three default false = legacy = byte-identical.
              shower_absorb_track_guard=false,
              shower_connect_protected_pion_guard=false,
              michel_stem_muon_rescue=false,
              // doc sbnd_xin/docs/pr/74 round 2 -- P1 cascade guard + P2
              // Michel-terminal check.  C++ defaults false/40cm/1.3/40cm.
              // Keys omitted when off/null => byte-identical pre-pr/74.
              shower_in_cascade_guard=false,
              shower_in_max_len=null,
              shower_in_mip_hi=null,
              // doc sbnd_xin/docs/pr/40 round 9 -- straight-track PID guard
              // family + B2 cross-cluster bridge.  C++ defaults false
              // (scalars 25 deg / 1.8 cm live in C++).  Keys omitted when
              // off/null => byte-identical pre-round-9 config.
              shower_connect_from_vertices_straight_guard=false,
              shower_connect_start_seg_straight_guard=false,
              examine_direction_dirsign_shower_in_guard=false,
              daughter_shower_angle_reclass_straight_guard=false,
              shower_topo_reexam_straight_guard=false,
              sfv_kink_max=null,
              shower_nv_bridge_track=false,
              shower_nv_bridge_max_gap=null,
              // doc pr/97 D1; C++ default false => legacy indeterminate main_pi read.
              shower_nv_main_pi_init=false,
              // doc pr/92 -- stray-satellite drop from kine/PF.  false/null =
              // C++ defaults = OFF (20 MeV / 8 cm / 60 deg / 45 deg / 90 cm /
              // 30 cm / 25 deg); keys omitted => byte-identical pre-pr/92.
              kine_drop_stray_satellites=false,
              kine_sat_min_energy=null,
              kine_sat_prox_max=null,
              kine_sat_angle_bad=null,
              kine_sat_angle_main=null,
              kine_sat_far_dis=null,
              kine_sat_axis_dis_cut=null,
              kine_sat_cont_kink=null,
              kine_sat_track_max_nseg=null,
              kine_sat_em_far_dis=null,
              michel_stem_michel_check=false,
              michel_stem_max_far_len=null,
              shower_stem_backfill=false,
              stem_backfill_max_len=null,
              stem_backfill_mip_lo=null,
              stem_backfill_mip_hi=null,
              stem_backfill_min_shower_len=null,
              shower_conn3_unreachable=false,
              conn3_unreachable_min_len=null,
              // doc pr/84 round 2 (F3): stitch radius in cm; null = C++
              // default 0 = OFF, key omitted => byte-identical.
              conn3_stitch_max=null,
              // doc pr/84 round 3 (S1): collapse showers that share a start
              // segment.  C++ default false = OFF, key omitted when off =>
              // byte-identical pre-fix config.
              shower_dedup_start_seg=false,
              shower_traj_michel_stem=false,
              michel_stem_traj_min_len=null,
              michel_stem_traj_max_len=null,
              michel_stem_traj_mip_lo=null,
              michel_stem_traj_max_far_len=null,
              michel_stem_traj_min_kink_deg=null,
              // doc pr/44: a multi-segment long-muon pseudo-shower keeps its
              // muon start segment (prototype parity; the update_particle_type
              // majority vote at the in_main_cluster seeding site is a
              // toolkit-only addition).  C++ default false; key omitted when
              // off => byte-identical pre-fix config.
              shower_long_muon_keep_type=false,
              // doc pr/40 round 10; false = C++ default = OFF, key omitted =>
              // byte-identical.
              shower_bragg_protect_start_segment=false,
              // doc pr/93 round 3 -- C++ defaults false; key-suppressed when off.
              shower_reclass_case_b_dqdx_guard=false,
              shower_accept_pid_guard=false,
              shower_pid_guard_min_len=null,
              shower_vote_track_pid_counts=false,
              shower_cone_absorb_guard=false,
              // doc pr/93 round 4 -- C++ defaults false / null = C++ defaults
              // (orphan floor 50cm; sccc 5cm/15deg base + 12cm/7.5deg aligned).
              // Keys suppressed when off => byte-identical.
              shower_detach_track_stem=false,
              shower_ghost_member_drop=false,
              shower_ghost_overlap_frac=null,
              shower_ghost_dqdx_ratio=null,
              shower_ghost_min_len=null,
              // doc pr/99 round 3; C++ defaults; keys suppressed when off.
              kine_charge_dedup=false,
              kine_charge_rebuild=false,
              // doc pr/101 Enu accounting round; C++ defaults; keys suppressed when off.
              kine_charge_track_ctx=false,
              kine_mass_rules=false,
              kine_hadronic_dqdx=false,
              kine_long_muon_mode=null,
              kine_long_muon_ratio_lo=null,
              kine_long_muon_ratio_hi=null,
              kine_mainvtx_used_guard=false,
              shower_hadronic_tag=false,
              shower_hadronic_min_len=null,
              shower_hadronic_scan_len=null,
              shower_hadronic_bin=null,
              shower_hadronic_r_cyl=null,
              shower_hadronic_r_core=null,
              shower_hadronic_growth_max=null,
              shower_hadronic_growth_bragg=null,
              shower_hadronic_bragg_ratio=null,
              shower_hadronic_stem_ratio=null,
              kine_count_orphan_tracks=false,
              kine_orphan_track_min=null,
              straight_cont_cross_cluster=false,
              sccc_bridge_body=false,
              sccc_max_gap=null,
              sccc_kink_max=null,
              sccc_gap_aligned=null,
              sccc_kink_tight=null,
              // doc pr/43 round 2 -- C++ defaults false; keys suppressed when off.
              single_muon_proton_chain_veto=false,
              single_muon_long_muon_claim=false,
              pid_flag_reconcile=false,
              // doc pr/45 -- find_other_segments empty-2D-tree sentinel guard
              // (SBND 18255-56463 isochronous tail).  C++ default false; key
              // suppressed when off => byte-identical.
              other_seg_empty_2d_guard=false,
              // doc pr/46 -- long-muon stub bridge in find_cont_muon_segment
              // (18255-55595 broken muon behind a 2.4 cm vertex stub).  C++
              // default false; key suppressed when off => byte-identical.
              long_muon_stub_bridge=false,
              // doc pr/48 -- back-to-back track fixes
              // (18255-51513/56211/57903/59335/57485: nu vertex mid-segment
              // on one unbroken track, dQ/dx rising at BOTH ends, no angular
              // kink).  two_end_break = the two-end residual-range break
              // pass; kink_walk_dqdx_stop / kink_break_protect = the 59335
              // walk-overshoot + EV4-absorption fixes.  C++ defaults false;
              // keys suppressed when off => byte-identical.  teb_* numerics:
              // null = C++ defaults (doc pr/48 sec 9 operating point), inert
              // while two_end_break is off.
              two_end_break=false,
              teb_min_len=null,
              teb_min_arm=null,
              teb_min_arm_pts=null,
              teb_stub_max=null,
              teb_accept_range=null,
              teb_rise_r1=null,
              teb_rise_r2=null,
              teb_abs_end_min=null,
              teb_dip_floor=null,
              teb_score_cap_r1=null,
              teb_score_cap_r2=null,
              teb_turn_angle=null,
              teb_turn_baseline=null,
              teb_turn_skirt=null,
              // doc pr/90 round 2: R2 argmax arm-fill guard + second-prong
              // gate cap.  null = C++ default 0 = legacy, key suppressed.
              teb_turn_min_arm_frac=null,
              teb_second_max=null,
              // doc pr/90 round 4: chain-topology admission (D1), route R3
              // turn/activity (D3), R2 bragg veto (D4).  false/null = C++
              // default = legacy, key suppressed.
              teb_chain_topology=false,
              teb_r3_turn=null,
              teb_r3_hot=null,
              teb_bragg_veto_turn=null,
              kink_walk_dqdx_stop=false,
              kink_break_protect=false,
              kink_dqdx_hot_ratio=null,
              // doc pr/49 -- cross-cluster projection-ghost deweighting in
              // the trajectory fit's 2D charge association (18255-57441
              // V-plane ghost: an unrelated cluster's charge aliases with
              // the fitted track in one view and detours the fit; live cells
              // outside the fitted cluster's own blob coverage that sit
              // inside a 3D-distant foreign cluster's keep their measurement
              // at reduced weight -- cells covered by nobody, or claimed by
              // a genuinely touching neighbor, keep full weight, so
              // no-ghost events are untouched).  C++ default -1 = off; null
              // omits the key => byte-identical.  >= 0 = on, value =
              // wire/slice tolerance in cells (0 = strict; the 57441
              // contamination is ONE cell away, so >= 1 re-admits it).
              fit_blob_coverage=null,
              // doc pr/50 -- suspend the pr/49 deweighting during
              // find_proto_vertex (the partition-forming stage; its recursive
              // kink walk is globally sensitive to fit perturbations --
              // 18255-172230 lost its true-kink main vertex to a 2.7 cm
              // neighbor).  All later fitting stages keep the deweighting.
              // C++ default false = pr/49 behavior; false omits the key =>
              // byte-identical.
              fit_blob_coverage_defer=false,
              // doc pr/50 -- main-vertex kink-consistency snap (172230-class
              // near-vertex robustness; C++ defaults in TaggerCheckNeutrino.h).
              // false/null omit the keys => byte-identical.
              vertex_kink_snap=false,
              vks_radius=null,
              vks_min_dis=null,
              vks_angle=null,
              vks_margin=null,
              vks_collinear=null,
              vks_skirt=null,
              vks_baseline=null,
              vks_min_arm=null,
              vks_fit_miss=null,
              vks_hot_ratio=null,
              // doc pr/85 -- carry the old vertex's arms through the snap
              // residual below this arc (cm).  C++ default 0 = off; null
              // omits the key => byte-identical.
              vks_carry_prong=null,
              // doc pr/104 -- main-vertex junction snap (C++ defaults in
              // TaggerCheckNeutrino.h).  false/null omit the keys => byte-identical.
              vertex_junction_snap=false,
              vjs_radius=null,
              vjs_min_arm=null,
              vjs_min_prongs=null,
              vjs_collinear=null,
              vjs_fit_margin=null,
              vjs_fit_rms=null,
              vjs_override_kink_snap=false,
              vjs_min_move=null,
              // doc pr/51 -- main-vertex graph audit (near-vertex graph-shape
              // repair: dup-corridor merge / charge-less-bridge removal /
              // micro-stub absorb + re-seat / one refit; C++ defaults in
              // TaggerCheckNeutrino.h).  false/null omit the keys =>
              // byte-identical.
              // esva_ignore_empty_2d (docs/73 sec 12, round 3, evt 78242):
              // eliminate_short_vertex_activities case 5 must not read the
              // empty-2D-index sentinel (-1) as "covered" on cathode-crossing
              // clusters.  SBND PRODUCTION ON since 2026-08-17 (docs/73
              // sec 12 owner flip).
              esva_ignore_empty_2d=true,
              main_vertex_graph_audit=false,
              mvga_radius=null,
              mvga_dup_tol=null,
              mvga_dup_frac=null,
              mvga_dup_angle=null,
              mvga_bridge_mip=null,
              mvga_reconnect=null,
              mvga_stub=null,
              mvga_stub_pts=null,
              mvga_reseat_angle=null,
              // doc pr/51 round 3 -- op3 satellite-anchor radius (cm).
              // C++ default 0 (main-vertex-only, round 2); null omits the
              // key => byte-identical.
              mvga_satellite=null,
              // doc pr/85 -- op3 interposed-stub absorb at the main-vertex
              // anchor (angle in deg, C++ default 150).  false/null omit
              // the keys => byte-identical.
              mvga_interposed=false,
              mvga_interposed_angle=null,
              // doc pr/86: interposed-splice candidate ceiling, cm (C++
              // default 0 = use mvga_stub).  null omits the key =>
              // byte-identical.
              mvga_interposed_len=null,
              // doc pr/86 P4: satellite-anchor op3 overlap threshold (C++
              // default 0 = use mvga_dup_frac).  null omits the key =>
              // byte-identical.
              mvga_sat_dup_frac=null,
              // doc pr/86 P1b: interposed splice at degree-1 main anchors
              // (C++ default false).  false omits the key => byte-identical.
              mvga_interposed_deg1=false,
              // doc pr/86 round 2: op3 post-carry straighten reach (cm) and
              // op3.5 junction-collapse radius (cm).  C++ defaults 0 (off).
              // Keys omitted when null => byte-identical.
              mvga_splice_straighten=null,
              mvga_approach_collapse=null,
              mvga_straighten_radius=null,
              // doc pr/83 r3: op1 scope/threshold decouple (cm / fraction;
              // radius -1 = unscoped), post-op3 dup pass, carry cap,
              // abandoned-cluster dup audit.  C++ defaults 0/0/false/0/false
              // (all legacy).  Keys omitted when null/false =>
              // byte-identical.
              mvga_op1_radius=null,
              mvga_op1_dup_frac=null,
              mvga_op1_post=false,
              mvga_carry_max=null,
              swap_orphan_dup_audit=false,
              // doc pr/83 r4 -- projective dup collapse; null omits => byte-identical.
              mvga_proj_dup_frac=null,
              mvga_proj_dqdx_ratio=null,
              mvga_proj_angle=null,
              mvga_ac_veto_radius=null,
              mvga_ac_chord_max=null,
              mvga_ac_no_cascade=false,
              mvga_passthru=null,
              mvga_passthru_tol=null,
              mvga_interposed_fallback=false,
              mvga_interposed_fallback_min_angle=null,
              mvga_dup_starved_asym=null,
              mvga_dup_starved_mip=null,
              mvga_dup_starved_span=null,
              // doc pr/51 (18255-506746) -- DL rerank cross-cluster swap
              // guard.  C++ default false; false omits the key =>
              // byte-identical.
              dl_vtx_swap_guard=false,
              // doc pr/106 sec 10: exclusion-free charge cloud for the DL vertex net. C++ default false.
              dl_vtx_cloud_no_exclusion=false,
              // doc pr/89 Arm C (C2) -- rule-1 topology term in the DL
              // rerank composite (weight, frac center).  C++ defaults 0/0;
              // null omits the keys => byte-identical.
              dl_vtx_topo_weight=null,
              dl_vtx_topo_center=null,
              // doc pr/51 round 3 -- apply the traditional-path swap
              // decision instead of discarding it.  C++ default false;
              // false omits the key => byte-identical.
              main_vertex_swap_apply=false,
              // doc pr/51 round 4 -- diagnostic-only rough-path probe.
              // C++ default false; false omits the key => byte-identical.
              rough_path_probe=false,
              // doc pr/51 round 5 -- steiner gap penalty (H1 short-cut fix):
              // do_rough_path routes on the support-penalized
              // "steiner_graph_gap" flavor when scale > 0.  C++ defaults:
              // scale 0 (off), dead_alpha 0.25, min_edge 0.5 cm,
              // sample_step 0.3 cm, point_radius 0.2 cm.  null omits the
              // keys => byte-identical.
              steiner_gap_penalty=null,
              sgp_dead_alpha=null,
              sgp_min_edge=null,
              sgp_sample_step=null,
              sgp_point_radius=null,
              // doc pr/73: per-edge DEBUG sentinel for the steiner_graph_gap
              // scan (endpoints, midpoint, w, bad, both vertex charges,
              // deficit).  Log-only diagnostic.  C++ default false; false
              // omits the key => byte-identical compiled config.
              sgp_edge_probe=false, vertex_scoreboard=false, dl_vtx_harvest=false,
              // doc pr/51 round 6 -- weak-charge deficit term on the same
              // gap flavor.  C++ defaults: weak_scale 0 (off), weak_qref
              // 2000 (charge units).  null omits the keys => byte-identical.
              sgp_weak_scale=null,
              sgp_weak_qref=null,
              // doc pr/73 round 2 F3a -- do_rough_path route excursion cap,
              // cm.  C++ default -1 = off; null omits the key.
              sgp_max_sep=null,
              // doc pr/83 -- oriented break_segment splits (find_vertices, not
              // boost source/target).  C++ default false; key omitted when
              // off => byte-identical pre-fix config.
              break_seg_orient=false,
              // doc pr/51 round 7 -- robust vertex fit (dynamic per-leg
              // direction windows for MyFCN).  C++ defaults: robust false,
              // main_only true, min_len 10, rin_margin 2, rout_frac 0.5,
              // rout_min 9, rout_max 18, angle 20, min_pts 5, min_aniso 3,
              // prior_range 1 (lengths cm, angle deg).  false/null omit the
              // keys => byte-identical.
              mvfit_robust=false,
              mvfit_main_only=null,
              mvfit_min_len=null,
              mvfit_rin_margin=null,
              mvfit_rout_frac=null,
              mvfit_rout_min=null,
              mvfit_rout_max=null,
              mvfit_angle=null,
              mvfit_min_pts=null,
              mvfit_min_aniso=null,
              mvfit_prior_range=null,
              // doc pr/54 -- keep well-supported isolated residual segments
              // in find_other_segments (18255-142421 separated EM shower with
              // no fitted trajectory).  C++ defaults: keep false, floors
              // 25 points / 3 cm.  false/null omit the keys => byte-identical.
              other_seg_keep_isolated=false,
              other_seg_keep_isolated_min_points=null,
              other_seg_keep_isolated_min_length=null,
              // doc pr/102 P1 -- OR-disjuncts on the pr/54 keep: min_nnf
              // (terminal not-faked floor) and len_admit (cm).  C++ defaults
              // 0 = off.  null omits the keys => byte-identical.
              other_seg_keep_isolated_min_nnf=null,
              other_seg_keep_isolated_len_admit=null,
              // doc pr/102 P2 -- 3-D uncovered-charge radius (cm) for the
              // find_other_segments tagging/nnf seats.  C++ default 0 = off.
              // null omits the key => byte-identical.
              other_seg_uncover_3d=null,
              // doc pr/67 round 3 (S2) -- isochronous-snap size gate, cm.
              // C++ default 10.0 = legacy.  null omits the key.
              iso_snap_min_dir_mag=null,
              // doc pr/65 round 3 -- offer graph-unreachable main-cluster
              // segments (kept-isolated pr/54 residuals) to the shower
              // absorbers (reachability-relaxed guards).  C++ default false;
              // false omits the key => byte-identical.
              shower_absorb_unreachable_main=false,
              // doc pr/59 round 2 -- per-cluster orphaned-associate_points
              // rescue.  C++ default false; false omits the key =>
              // byte-identical.
              assoc_full_recluster=false,
              // doc pr/64 round 7 -- reassign same-cluster association
              // orphans that Stage C of clustering_points_segments would
              // otherwise drop (18259-18625: 12-18 pt blob at PF segment
              // 126042's own fit endpoint, in img charge but absent from
              // shower_track/associate_points).  C++ default false; false
              // omits the key => byte-identical.
              assoc_reassign_orphans=false,
              // doc pr/64 round 8 -- clear a merge survivor's
              // associate_points when examine_structure_final_1/_1p/_3
              // deletes a segment that had non-empty associate_points, so
              // pr/59's reassociate_cluster_orphans any_orphan trigger
              // correctly re-fires.  C++ default false; false omits the key
              // => byte-identical.
              assoc_clear_on_merge=false,
              // doc pr/72 round 2 -- guard examine_structure_3 against
              // merging a genuine near-vertex track stub into an unrelated
              // shower/track trunk (18255-196649).  C++ default false;
              // false/null omits the key => byte-identical.  Numeric
              // defaults null = component keeps its own C++ default
              // (fitted from a 117-event census).
              es3_stub_guard=false,
              es3sg_stub_max=null,
              es3sg_len_ratio=null,
              es3sg_ang3_min=null,
              es3sg_ang_ratio=null,
              es3sg_require_terminal=null,
              // doc pr/45 -- paint muon-typed (+-13) pseudo-showers as track in
              // the Bee shower_track layer + PrDisplayDump (18255-56463: 411 cm
              // muon painted red).  C++ default false; key suppressed when off
              // => byte-identical.
              pseudo_shower_track_paint=false,
              // muon_dqdx_curve [c0, c1, pivot_cm, power]: the muon
              // median-dQ/dx-vs-length envelope used by nine tagger cuts, as a
              // multiple of mip_dqdx_median.  DEFAULT = the SBND fit
              // (2026-07-30, docs/pr/10 sec 4): the SBND stopping-muon table
              // median (0.5 kV/cm, /48000) times the uBooNE empirical/table
              // margin g(L), refit with the prototype's functional form.
              // NOT bit-identical: nine tagger cuts move (looser at short L
              // -- the higher field keeps more of the Bragg rise in dQ/dx).
              // Scales as 1/mip_dqdx_median: re-derive when the 48000
              // placeholder becomes a measurement.  null restores the uBooNE
              // refit [0.8866, 0.9533, 18, 0.4234] (byte-identical pre-knob).
              muon_dqdx_curve=[0.8826, 1.0587, 18, 0.4745],
              // use_power_recomb: hand the taggers (STM + neutrino PR) the
              // free-power Modified-Box recombination fitted to SBND stopping
              // tracks (docs/55 sec 7g canonical, PowerBoxRecombination
              // defaults) instead of the plain sbnd_box_recomb above.
              // DEFAULT ON for SBND (owner 2026-07-30, docs/pr/10 sec 7):
              // every model-driven dQ/dx -> dE/dx conversion now uses the
              // measured SBND recombination.  false restores sbnd_box_recomb
              // (the byte-identical pre-pr/10 config).
              use_power_recomb=true,
              // sp_dedx_use_recomb_model: route the single-photon stem dE/dx
              // through the configured recombination model instead of the
              // inline uBooNE-field (0.273 kV/cm) inverse Box.  DEFAULT ON
              // for SBND (owner 2026-07-30) with sp_mean_dedx_cut=2.23, the
              // physical-scale transfer of the legacy 2.3 (docs/pr/10 sec 5).
              // false/null restore the inline formula and the 2.3 literal.
              sp_dedx_use_recomb_model=true,
              sp_mean_dedx_cut=2.23,
              // dl_vtx_cut: max distance (mm; C++ default 25.0 = 2.5 cm) from
              // the DL SCN prediction to accept a candidate vertex.  Threaded
              // for configurability (docs/pr/2 sec 7.4); null keeps the C++
              // default, which is coupled to the uBooNE-trained net (gap G3).
              dl_vtx_cut=null) = {
    // Only gate when the caller actually supplied a window; beam_window=[0,0]
    // (the arg default, i.e. "no beam window") must not silently drop every
    // cluster's tagger evaluation.
    local beam_gate = beam_window_only && beam_window[1] > beam_window[0],
    local dv = detector_volumes(anodes, '', pos_offset_on),
    local pcts = pctransforms(dv),
    // DetectorVolumes implements IFiducial -- used by MakeFiducialUtils / the
    // taggers' inside_fiducial_volume().  NOT the FV_* metadata below:
    // DetectorVolumes::contained() is contained_by(p).valid(), i.e. the union of
    // the per-face IAnodeFace::sensitive() boxes, which AnodePlane builds as
    // x in [anode_x, cathode_x] and so run to the W plane with no margin
    // (SBND |x| <= 201.45, |y| <= 199.965, z in [0, 501.0]; |x| < 0.45 is a hole).
    // Any tagger given no explicit fiducial therefore tests a volume ~3 cm more
    // permissive at every wall than sbnd_pr_fv + sbnd_pr_fv_margins below --
    // TaggerCheckSTM and TaggerCheckNeutrino's match_isFC still do; see
    // sbnd_xin/docs/49_stm-containment-fv-inconsistency.md.
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
    // Free-power Modified Box fitted to SBND stopping-track dQ/dx vs residual
    // range (docs/55 sec 7g; canonical parameters in stm_ref_dqdx.json:
    // R = ln(A+u)/u, u = k*(dEdx/2.1)^p, times normalization C).  The data
    // block restates the C++ defaults so the operating point is visible here;
    // selected via use_power_recomb (docs/pr/10).
    local sbnd_power_recomb = {
        type: 'PowerBoxRecombination',
        name: 'sbnd_power_recomb',
        data: { A: 0.93, k: 0.282371, p: 1.362179, C: 0.855175, pivot: 2.1, Wi: 23.6e-6, dedx_max: 77.0 },
    },
    // The recombination model the taggers actually receive.  With
    // use_power_recomb=false this compiles to exactly the pre-knob config.
    local sbnd_recomb = if use_power_recomb then sbnd_power_recomb else sbnd_box_recomb,
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
    // Margin vector convention (inside_fv / Clustering_Util fc_check): index 4
    // is applied as contained(z - tv[4]), so with a NEGATIVE entry it insets
    // the DOWNSTREAM (z ~ 500 cm) face; index 5 insets the upstream face.
    // tgm_fv_zmax_margin (cm; default 3 = byte-identical legacy value)
    // parametrizes only the downstream inset -- shared by tagger_check_tgm AND
    // tagger_check_fc below, so "contained" keeps one meaning across both
    // verdicts (docs/27_fc-tgm-consistent-fv.md).
    // tgm_fv_x_margin / tgm_fv_y_margin (cm; defaults 2 / 2.5 = byte-identical
    // legacy values) parametrize the drift-x and vertical-y insets, both faces
    // symmetric.  Shared by tagger_check_tgm AND tagger_check_fc, same as the
    // downstream-z knob, so "contained" keeps one meaning across both verdicts.
    local sbnd_pr_fv_margins = [-tgm_fv_x_margin * wc.cm, -tgm_fv_x_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_zmax_margin * wc.cm, -3 * wc.cm],
    // tgm_fv_zmax_margin_interior (cm; default 0 = OFF, key omitted =>
    // byte-identical): when > 0, check_tgm's CASE-A interior-support tests
    // (chord midpoints + waypoint re-check) use THIS downstream-z inset
    // instead of tgm_fv_zmax_margin, i.e. the doc-32 widening becomes
    // endpoint-only.  Rationale (doc 32 caveat, doc 35): a corner clipper
    // running ALONG the downstream wall inside the widened 3->5 cm band
    // (evt287517 cluster 16, evt289805 cluster 9) keeps its midpoint support
    // at the legacy 3 cm interior.  TGM only -- tagger_check_fc containment
    // and the ENDPOINT exit tests keep sbnd_pr_fv_margins unchanged.
    // x/y track the endpoint vector above: the doc-35 endpoint-only widening
    // applies to the downstream-z inset only.
    local sbnd_pr_fv_margins_interior = [-tgm_fv_x_margin * wc.cm, -tgm_fv_x_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_zmax_margin_interior * wc.cm, -3 * wc.cm],
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
        // Restore the prototype main+associated data product before anything
        // fits or walks a "main cluster".  The Q/L stage's flash-time merge
        // (examine_bundles use_flash_t0) collapses every member of a matched
        // bundle into ONE flag_main_cluster cluster; uBooNE hands the PR chain
        // a single connected main plus a vector of companions
        // (wire-cell-prod-stm.cxx: main_cluster + additional_clusters).  Until
        // this split runs, TaggerCheckSTM fits a bundle of detached cosmics as
        // one track, check_other_clusters() has no companions to count, and
        // TaggerCheckNeutrino's matched_flash_gid companion list is empty.
        // Placed between switch_scope and steiner ON PURPOSE: separate() does
        // not carry node-local PCs, so the split must precede steiner_pc
        // creation (the visitor erases a stale steiner_pc and warns if not).
        // Not in pipeline_names => absent from the compiled config => no
        // behavior change.  Runner flag: -unmerge / -unmerge-comp.
        unmerge_bundle: cm.unmerge_bundle(mode=unmerge_bundle_mode,
                                          restore_demoted_mains=restore_demoted_mains,
                                          require_provenance=require_provenance),
        // Second, INNER un-merge: undo the per-APA isolated GROUPING
        // (clustering_isolated save_assoc_id), which merges a main cluster with
        // the small clusters that are near it but NOT connected to it.  The
        // prototype only groups these (Clustering_isolated returns
        // main -> [(assoc, dis)] and leaves live_clusters untouched); the toolkit
        // physically merges them, so the STM/PR endpoint finder walks into a
        // detached clump across empty space (docs 50, 51).
        //
        // Order matters: this runs AFTER unmerge_bundle so the flash grouping is
        // undone first (outer) and the isolated grouping second (inner) --
        // together they reproduce the prototype's main_cluster +
        // additional_clusters exactly.  Both must precede steiner.
        //
        // A cathode crosser survives both: its two halves are each the MAIN of
        // their own per-APA isolated group (cm.isolated() is per-APA only), and
        // the visitor retains the union of every main-marked member, splitting
        // off only assoc_cluster_main == 0.  Doc 52 4a/4b.
        //
        // Not in pipeline_names => absent from the compiled config.
        unmerge_assoc: cm.unmerge_bundle(name='assoc', mode=unmerge_bundle_mode,
                                         id_aname='assoc_cluster_id',
                                         main_aname='assoc_cluster_main'),
        // PR-stage overclustering protection (doc pr/23): uboone's SECOND
        // graph-examination round (WCPPID::Protect_Over_Clustering ->
        // Examine_graph, run after Q/L matching and before NeutrinoID) --
        // split each beam-bundle cluster at graph-component boundaries so
        // vertex-proximate photon fragments merged by Clustering_neutrino are
        // fit as separate clusters instead of one MST-bridged trajectory
        // (doc pr/22 sec 8, evt 386948).  Position (doc pr/23 ordering
        // decision): AFTER the cosmic taggers (tagger_check_tgm/stm/fc) and
        // BEFORE tagger_check_neutrino, with 'steiner' named a second time
        // right after it so the split clusters' steiner products are rebuilt
        // (the visitor drops the pre-split ones; CreateSteinerGraph
        // replace=false touches nothing else).  This matches uboone -- cosmic
        // verdicts on unsplit clusters (wire-cell-prod-stm.cxx:806-810),
        // Protect_Over_Clustering only in the nue executable
        // (wire-cell-prod-nue.cxx:1322).  The cathode re-join
        // knobs keep it from splitting cathode crossers, an SBND geometry
        // uboone did not have (doc pr/20).  In the production pipeline_names
        // by DEFAULT since doc pr/23 sec 9; dropped from the list => absent
        // from the compiled config => the pre-pr/23 chain.
        protect_bundle: cm.protect_bundle(
            graph_name=if protect_graph_name == null then 'relaxed' else protect_graph_name,
            beam_window_only=beam_gate,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1],
            skip_convicted=protect_skip_convicted,
            open_convicted_bundles=protect_open_convicted_bundles,
            cathode_x=protect_cathode_x,
            cathode_rejoin_xcut=protect_cathode_rejoin_xcut,
            cathode_rejoin_dyz=protect_cathode_rejoin_dyz,
            cathode_rejoin_dis=protect_cathode_rejoin_dis,
            cathode_rejoin_perp=protect_cathode_rejoin_perp,
            cathode_rejoin_angle=protect_cathode_rejoin_angle,
            cathode_rejoin_dir_radius=protect_cathode_rejoin_dir_radius,
            cathode_rejoin_dir_npts=protect_cathode_rejoin_dir_npts),
        // SBND has no beam_flash flag (QLMatching sets main/associated_cluster
        // instead) -- process every scope-passing cluster, narrowed to the
        // beam-coincident bundle when beam_window_only is on (the default; see
        // the clus_pr argument).  Keys omitted when off => byte-identical.
        steiner: cm.steiner(retiler=improve2, perf=true, require_beam_flash=false,
                            beam_window_only=beam_gate,
                            beam_window_low=beam_window[0],
                            beam_window_high=beam_window[1],
                            terminal_wire_tol=steiner_terminal_wire_tol,
                            terminal_adjacent_slice=steiner_terminal_adjacent_slice,
                            edge_charge_forward_dead_mix=steiner_edge_charge_forward_dead_mix),
        // The doc pr/23 second steiner pass, named right after protect_bundle:
        // replace=false rebuilds ONLY the clusters protect_bundle purged
        // (split retained + fragments).  A replace=true second pass would
        // erase+emplace every in-window cluster's steiner_graph and dangle
        // the GraphAlgorithms the STM fit cached in the tagger stage
        // (bad_alloc from garbage num_vertices, SBND evt 54095).  Distinct
        // component name => distinct instance; in the production
        // pipeline_names by DEFAULT since doc pr/23 sec 9, always paired
        // with protect_bundle.
        steiner_refresh: cm.steiner(name='refresh',
                            retiler=improve2, perf=true, require_beam_flash=false,
                            beam_window_only=beam_gate,
                            beam_window_low=beam_window[0],
                            beam_window_high=beam_window[1],
                            // Same terminal-filter settings as the first pass:
                            // the refresh rebuilds only the clusters
                            // protect_bundle purged, and they must be built the
                            // same way as their peers (doc pr/29 D1, D12).
                            terminal_wire_tol=steiner_terminal_wire_tol,
                            terminal_adjacent_slice=steiner_terminal_adjacent_slice,
                            edge_charge_forward_dead_mix=steiner_edge_charge_forward_dead_mix,
                            replace=false),
        fiducialutils: cm.fiducialutils(),
        tagger_check_stm: cm.tagger_check_stm(
            evaluate_demoted_mains=evaluate_demoted_mains,
            trackfitting_config_file=trackfitting_config_file,
            particle_dataset=wc.tn(particle_dataset),
            recombination_model=wc.tn(sbnd_recomb),
            require_in_scope=true,
            // save_stm_fit (C++ default false; key omitted when off =>
            // byte-identical): persist per-pass STM fits as cluster PCs +
            // grouping slot "stm" for the Bee stm_fit layer and the
            // stm_magnify ROOT dump below.  Runner flag: -stm-fit.
            save_stm_fit=save_stm_fit,
            // mip_dqdx: SBND MIP dQ/dx scale in e/cm, replacing the inherited
            // MicroBooNE 50000.  Anchored to the muon reference table the same
            // way 50000 was anchored to MicroBooNE's: the uBooNE table plateau
            // (rr = 59.5 cm) is 48879.4 e/cm and 50000/48879.4 = 1.02293; the
            // SBND table regenerated at 0.5 kV/cm plateaus at 54657.7 e/cm and
            // 56000/54657.7 = 1.02456, so the same 2.3-2.5% headroom is kept
            // and 56000 is as round a number as 50000 was.
            // NOT byte-identical -- see sbnd_xin/docs/48.
            mip_dqdx=mip_dqdx,
            // accept_guards (C++ default false; key omitted when off =>
            // byte-identical): doc-63 round-1 acceptance guards.
            accept_guards=stm_accept_guards,
            // proton_muon_guard (C++ default false; key omitted when off =>
            // byte-identical): doc-63 round-2 muon-consistency guard.
            proton_muon_guard=stm_proton_muon_guard,
            // cathode_guard (C++ default false; key omitted when off =>
            // byte-identical): doc-63 round-3 cathode-truncation veto.
            cathode_guard=stm_cathode_guard,
            // anode_dist_fix (C++ default false; key omitted when off =>
            // byte-identical): doc-63 round-4a dist_to_anode face fix.
            anode_dist_fix=stm_anode_dist_fix,
            // second_track_guard (C++ default false; key omitted when off =>
            // byte-identical): doc-63 round-4b/4c second-track vetoes.
            second_track_guard=stm_second_track_guard,
            // deficit_guard / vertex_kink_guard (C++ defaults false; keys
            // omitted when off => byte-identical): doc-63 round-5 vetoes.
            deficit_guard=stm_deficit_guard,
            vertex_kink_guard=stm_vertex_kink_guard,
            // doc-66 sec 12 cut package (C++ defaults = prototype constants;
            // keys omitted when stm_d66_cuts=false => byte-identical).
            michel_res_length_cut=(if stm_d66_cuts then stm_michel_res_cm * wc.cm else null),
            proton_tm_max=(if stm_d66_cuts then stm_proton_tm_max else null),
            proton_b_ks2_max=(if stm_d66_cuts then stm_proton_b_ks2_max else null),
            proton_c_peak_max=(if stm_d66_cuts then stm_proton_c_peak_max else null),
            // stm_consistent_fv (default TRUE = SBND production): give the
            // containment gate the SAME fiducial + margins tagger_check_tgm and
            // tagger_check_fc get, so "contained" means one thing across all
            // three verdicts.  Without it cluster_fc_check falls back to
            // FiducialUtils -> the union of per-face SENSITIVE volumes, which
            // exceeds this box at every wall even before the margins (SBND 0.40
            // cm x / 0.65 y / 0.85 z) and is holed at the CPA slab: on a
            // 30-event sample 96 of the 147 clusters the STM tagger skipped as
            // "fully contained" were called exiters by tagger_check_fc, so no
            // dQ/dx fit was ever attempted for them.  The prototype has no such
            // split -- check_stm and check_tgm are members of ONE ToyFiducial
            // and share its boundaries, inset once by boundary_dis_cut = 3 cm --
            // so this restores parity rather than inventing a volume.
            // NOT byte-identical.  Set false for the A/B (runner: -no-stm-fv).
            // See sbnd_xin/docs/49_stm-containment-fv-inconsistency.md.
            // Only the ENDPOINT tests exist in cluster_fc_check, so the
            // interior_fv_tolerance vector tagger_check_tgm needs has no
            // counterpart here.
            fiducial=(if stm_consistent_fv then wc.tn(sbnd_pr_fv) else null),
            fv_tolerance=(if stm_consistent_fv then sbnd_pr_fv_margins else []),
            // Beam-window gate on the main-cluster loop; see the clus_pr
            // beam_window_only argument.  Keys omitted when off =>
            // byte-identical.
            beam_window_only=beam_gate,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1]),
        // STM-stage Magnify-tracking ROOT dump (doc sbnd_xin/docs/40): reads
        // the stm_fit/stm_pass cluster PCs and the "stm" TrackFitting slot,
        // writes tracking-stm.root (T_rec_charge/T_proj_data/T_bad_ch/Trun)
        // with the two-TPC concatenated-per-plane channel convention.  Only
        // active when named in pipeline_names (needs -stm-fit; the WireCellRoot
        // plugin must be loaded by the job).
        stm_magnify: {
            type: 'SbndMagnifyTrackingVisitor',
            name: 'pr',
            data: {
                grouping: 'live',
                track_fitting_name: 'stm',
                output_filename: (if output_dir == '' then '' else output_dir + '/') + 'tracking-stm.root',
                runNo: runNo,
                subRunNo: subRunNo,
                eventNo: eventNo,
                anodes: [wc.tn(a) for a in anodes],
                detector_volumes: wc.tn(dv),
                dQdx_scale: 0.1,
                dQdx_offset: -1000.0,
                // Readout length in ticks; only used to clamp T_bad_ch time
                // ranges (a dead region with unbounded x converts to +-1e9
                // ticks and makes the Magnify projection pads unreadable).
                // C++ default 3427 -- same value, stated here for the record.
                nticks: 3427,
            },
        },
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
            evaluate_demoted_mains=evaluate_demoted_mains,
            fiducial=wc.tn(sbnd_pr_fv),
            fv_tolerance=sbnd_pr_fv_margins,
            // Key omitted when the knob is 0 (empty list) => byte-identical.
            interior_fv_tolerance=(if tgm_fv_zmax_margin_interior > 0
                                   then sbnd_pr_fv_margins_interior else []),
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1],
            // Beam-window gate on the main-cluster loop -- the SAME window the
            // in-beam protection above uses; see the clus_pr beam_window_only
            // argument.  Key omitted when off => byte-identical.
            beam_window_only=beam_gate,
            require_in_scope=true,
            check_neutrino_candidate=tgm_neutrino_candidate,
            require_chord_charge=tgm_chord_charge,
            // C++ default "chord". Key omitted then => byte-identical.
            chord_charge_mode=tgm_chord_mode,
            component_extremes=tgm_component_extremes,
            // C++ default false. Key omitted when off => byte-identical.
            component_rescue=tgm_component_rescue,
            // C++ default false. Key omitted when off => byte-identical
            // doc-32 rescue behavior.  Rescued-end pairs must also pass the
            // straight-chord test (evt288727 two-cosmic composite).
            rescue_chord_check=tgm_rescue_chord,
            // C++ default false. Key omitted when off => byte-identical.
            // A pair may tag only when an end lies in the main cluster --
            // the TGM verdict follows the main cluster, not a merged-in
            // through-going fragment (evt289343 cluster 9).
            main_component_pairs=tgm_main_pair,
            // C++ default "path" (largest-path-component proxy). Key omitted
            // then => byte-identical.  "real" = per-blob flash-merge
            // provenance (needs a save_real_cluster_id pctree; falls back to
            // the proxy on old tarballs).
            main_component_mode=tgm_main_pair_mode,
            // C++ default false. Key omitted when off => byte-identical.
            // See the clus_pr arg comment (doc pr/25, SBND evt 320029).
            exempt_demoted_main_pairs=tgm_exempt_demoted_main),
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
            evaluate_demoted_mains=evaluate_demoted_mains,
            fiducial=wc.tn(sbnd_pr_fv),
            fv_tolerance=sbnd_pr_fv_margins,
            require_in_scope=true,
            // Beam-window gate on the main-cluster loop; see the clus_pr
            // beam_window_only argument.  Keys omitted when off =>
            // byte-identical.
            beam_window_only=beam_gate,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1]),
        // Neutrino pattern recognition on the beam-coincident bundle.  The
        // beam_window gate (on cluster_t0 = matched flash time) replaces
        // uBooNE's single-main + beam_flash selection; companions are the
        // associated clusters sharing the main's matched_flash_gid.
        // dl_weights defaults to the uBooNE-trained SCN net = DL vertex ON
        // (docs/pr/4); '' falls back to the geometric vertex.
        tagger_check_neutrino: cm.tagger_check_neutrino(
            trackfitting_config_file=trackfitting_config_file,
            particle_dataset=wc.tn(particle_dataset),
            recombination_model=wc.tn(sbnd_recomb),
            perf=true,
            dl_weights=dl_weights,
            // The DL re-rank sub-knobs, pinned here rather than inherited: they
            // were inert while dl_weights was '' and went LIVE with the doc-pr/4
            // adoption, so SBND records the operating point it was validated at
            // (= the common/clus.jsonnet defaults as of 2026-07-30, hence the
            // compiled JSON is unchanged by this pinning).  min_accept and top_k
            // are threaded from the clus_pr args (defaults identical to the old
            // pinned literals => byte-identical when unset; doc pr/79).
            dl_vtx_rerank=dl_vtx_rerank,
            dl_vtx_top_k=dl_vtx_top_k,
            dl_vtx_min_accept_score=dl_vtx_min_accept_score,
            dl_vtx_score_scale=1000.0,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1],
            nu_skip_cosmic=nu_skip_cosmic,
            nu_skip_cosmic_bundle=nu_skip_cosmic_bundle,
            nu_skip_cosmic_bundle_min_length=nu_skip_cosmic_bundle_min_length,
            skip_cosmic_companions=skip_cosmic_companions,
            cosmic_companion_min_length=cosmic_companion_min_length,
            nu_fallback_demoted_mains=nu_fallback_demoted_mains,
            // doc pr/94 Phase 2: the SAME switch that books the per-bundle
            // T_tagger/T_kine branches on tagger_output below also turns on the
            // per-bundle candidate loop here -- one row per in-beam-window
            // bundle needs both halves or the extra rows have nothing to carry.
            // nu_per_bundle_demoted_acts is wired straight from
            // evaluate_demoted_mains (not from its own TLA) so the
            // act_evaluated column can never drift from the admission gate the
            // cosmic taggers actually used.
            nu_per_bundle=nu_per_bundle,
            nu_per_bundle_demoted_acts=evaluate_demoted_mains,
            nu_per_bundle_min_length=nu_per_bundle_min_length,
            nu_selected_as_main=nu_selected_as_main,
            nu_selected_as_main_snapshot_all=nu_selected_as_main_snapshot_all,
            sp_photon_flag=sp_photon_flag,
            dir_weak_use_score=dir_weak_use_score,
            mip_dqdx=mip_dqdx,
            mip_dqdx_median=mip_dqdx_median,
            proton_dir_vote=proton_dir_vote,
            endpoint_trim_retry=endpoint_trim_retry,
            fit_vertex_min_seg_length=fit_vertex_min_seg_length,
            cathode_x=cathode_x,
            cathode_kink_xcut=cathode_kink_xcut,
            cathode_wide_kink_angle=cathode_wide_kink_angle,
            cathode_wide_kink_skirt=cathode_wide_kink_skirt,
            cathode_wide_kink_baseline=cathode_wide_kink_baseline,
            shower_topo_demote_len=shower_topo_demote_len,
            fit_exclusion=fit_exclusion,
            graph_endpoint_strict=graph_endpoint_strict,
            graph_endpoint_tol=graph_endpoint_tol,
            oov_prototype_parity=oov_prototype_parity,
            first_seg_local_pca=first_seg_local_pca,
            other_seg_relaxed_accept=other_seg_relaxed_accept,
            shower_topo_proto_dir=shower_topo_proto_dir,
            cont_muon_dir3_30cm=cont_muon_dir3_30cm,
            track_comp_empty_abstain=track_comp_empty_abstain,
            shower_topo_reset=shower_topo_reset,
            reclass_preserve_4mom=reclass_preserve_4mom,
            dir_track_median_local=dir_track_median_local,
            examine_showers_vertex_by_index=examine_showers_vertex_by_index,
            vertex_dir_use_fit_point=vertex_dir_use_fit_point,
            shower_traj_recheck_parity=shower_traj_recheck_parity,
            main_vertex_require_descriptor=main_vertex_require_descriptor,
            main_vertex_candidate_flag=main_vertex_candidate_flag,
            iso_endpoint=iso_endpoint,
            iso_endpoint_min_length=iso_endpoint_min_length,
            iso_endpoint_max_xext=iso_endpoint_max_xext,
            iso_endpoint_xext_frac=iso_endpoint_xext_frac,
            iso_endpoint_xext_quantile=iso_endpoint_xext_quantile,
            iso_endpoint_tube_radius=iso_endpoint_tube_radius,
            iso_endpoint_min_aspect=iso_endpoint_min_aspect,
            v3_extension_guard=v3_extension_guard,
            traj_cover_probe=traj_cover_probe,
            pr_find_other_rounds=pr_find_other_rounds,
            v3_extension_min_gain=v3_extension_min_gain,
            cosmic_y_top_main=cosmic_y_top_main,
            cosmic_y_top_strict=cosmic_y_top_strict,
            cosmic_y_top_loose=cosmic_y_top_loose,
            cosmic_y_small_piece=cosmic_y_small_piece,
            vertex_z_prior_scale=vertex_z_prior_scale,
            ssm_target_dir=ssm_target_dir,
            ssm_absorber_dir=ssm_absorber_dir,
            kine_fudge_factor=kine_fudge_factor,
            kine_recom_factor=kine_recom_factor,
            kine_shower_fudge_factor=kine_shower_fudge_factor,
            kine_shower_recom_factor=kine_shower_recom_factor,
            kine_proton_recom_factor=kine_proton_recom_factor,
            kine_plane_weights=kine_plane_weights,
            kine_plane_asym_switch=kine_plane_asym_switch,
            kine_w_value=kine_w_value,
            kine_shower_pdg_live=kine_shower_pdg_live,
            // doc pr/36 sec 10 (F1): same fiducial + margins as
            // tagger_check_{stm,tgm,fc} above -- one containment definition
            // across the stage.  Keys omitted when off.
            fiducial=(if neutrino_consistent_fv || cosmic_consistent_fv || nue_sp_consistent_fv then wc.tn(sbnd_pr_fv) else null),
            fv_tolerance=(if neutrino_consistent_fv || cosmic_consistent_fv || nue_sp_consistent_fv then sbnd_pr_fv_margins else []),
            // sbnd_xin/docs/74 G1/G2: cosmic_tagger() containment on the same
            // fiducial + margins.  Key omitted when off => byte-identical.
            cosmic_consistent_fv=cosmic_consistent_fv,
            // sbnd_xin/docs/75: nue/single-photon tagger containment on the
            // same fiducial (each site's own hardcoded tolerance -- see
            // NeutrinoTaggerNuE.cxx).  Key omitted when off => byte-identical.
            nue_sp_consistent_fv=nue_sp_consistent_fv,
            sp_sce_correction=sp_sce_correction,
            tagger_ordered_segment_sets=tagger_ordered_segment_sets,
            stem_endpoint_wcpt_parity=stem_endpoint_wcpt_parity,
            broken_muon_cluster_id_count=broken_muon_cluster_id_count,
            neutrino_type_bitmask=neutrino_type_bitmask,
            daughter_count_proto_main_vertex=daughter_count_proto_main_vertex,
            daughter_count_proto_examine_showers=daughter_count_proto_examine_showers,
            shower_pdg_from_start_segment=shower_pdg_from_start_segment,
            shower_pdg_from_shower_type=shower_pdg_from_shower_type,
            shower_pdg_exact_muon_test=shower_pdg_exact_muon_test,
            pi0_id_shared_allocator=pi0_id_shared_allocator,
            shower_flag_pdg_electron=shower_flag_pdg_electron,
            shower_less_id_tiebreak=shower_less_id_tiebreak,
            shower_endpoint_exclude_start_vertex=shower_endpoint_exclude_start_vertex,
            shower_endpoint_skip_orphan_vtx=shower_endpoint_skip_orphan_vtx,
            shower_walk_visited_parity=shower_walk_visited_parity,
            track_pid_persist_dqdx=track_pid_persist_dqdx,
            shower_reclass_dqdx_guard=shower_reclass_dqdx_guard,
            shower_topo_dqdx_guard=shower_topo_dqdx_guard,
            reclass_never_computed_ke_floor=reclass_never_computed_ke_floor,
            track_pid_persist_4mom=track_pid_persist_4mom,
            shower_proton_daughter_pion=shower_proton_daughter_pion,
            shower_proton_daughter_pion_dissolve=shower_proton_daughter_pion_dissolve,
            muon_multi_proton_pion=muon_multi_proton_pion,
            track_pid_persist_dqdx_electron_guard=track_pid_persist_dqdx_electron_guard,
            shower_connect_main_vertex_straight_guard=shower_connect_main_vertex_straight_guard,
            shower_traj_straight_guard=shower_traj_straight_guard,
            shower_absorb_track_guard=shower_absorb_track_guard,
            shower_connect_protected_pion_guard=shower_connect_protected_pion_guard,
            michel_stem_muon_rescue=michel_stem_muon_rescue,
            shower_in_cascade_guard=shower_in_cascade_guard,
            shower_in_max_len=shower_in_max_len,
            shower_in_mip_hi=shower_in_mip_hi,
            shower_connect_from_vertices_straight_guard=shower_connect_from_vertices_straight_guard,
            shower_connect_start_seg_straight_guard=shower_connect_start_seg_straight_guard,
            examine_direction_dirsign_shower_in_guard=examine_direction_dirsign_shower_in_guard,
            daughter_shower_angle_reclass_straight_guard=daughter_shower_angle_reclass_straight_guard,
            shower_topo_reexam_straight_guard=shower_topo_reexam_straight_guard,
            sfv_kink_max=sfv_kink_max,
            shower_nv_bridge_track=shower_nv_bridge_track,
            shower_nv_bridge_max_gap=shower_nv_bridge_max_gap,
            shower_nv_main_pi_init=shower_nv_main_pi_init,
            kine_drop_stray_satellites=kine_drop_stray_satellites,
            kine_sat_min_energy=kine_sat_min_energy,
            kine_sat_prox_max=kine_sat_prox_max,
            kine_sat_angle_bad=kine_sat_angle_bad,
            kine_sat_angle_main=kine_sat_angle_main,
            kine_sat_far_dis=kine_sat_far_dis,
            kine_sat_axis_dis_cut=kine_sat_axis_dis_cut,
            kine_sat_cont_kink=kine_sat_cont_kink,
            kine_sat_track_max_nseg=kine_sat_track_max_nseg,
            kine_sat_em_far_dis=kine_sat_em_far_dis,
            michel_stem_michel_check=michel_stem_michel_check,
            michel_stem_max_far_len=michel_stem_max_far_len,
            shower_stem_backfill=shower_stem_backfill,
            stem_backfill_max_len=stem_backfill_max_len,
            stem_backfill_mip_lo=stem_backfill_mip_lo,
            stem_backfill_mip_hi=stem_backfill_mip_hi,
            stem_backfill_min_shower_len=stem_backfill_min_shower_len,
            shower_conn3_unreachable=shower_conn3_unreachable,
            conn3_unreachable_min_len=conn3_unreachable_min_len,
            conn3_stitch_max=conn3_stitch_max,
            shower_dedup_start_seg=shower_dedup_start_seg,
            shower_traj_michel_stem=shower_traj_michel_stem,
            michel_stem_traj_min_len=michel_stem_traj_min_len,
            michel_stem_traj_max_len=michel_stem_traj_max_len,
            michel_stem_traj_mip_lo=michel_stem_traj_mip_lo,
            michel_stem_traj_max_far_len=michel_stem_traj_max_far_len,
            michel_stem_traj_min_kink_deg=michel_stem_traj_min_kink_deg,
            shower_long_muon_keep_type=shower_long_muon_keep_type,
            shower_bragg_protect_start_segment=shower_bragg_protect_start_segment,
            shower_reclass_case_b_dqdx_guard=shower_reclass_case_b_dqdx_guard,
            shower_accept_pid_guard=shower_accept_pid_guard,
            shower_pid_guard_min_len=shower_pid_guard_min_len,
            shower_vote_track_pid_counts=shower_vote_track_pid_counts,
            shower_cone_absorb_guard=shower_cone_absorb_guard,
            shower_detach_track_stem=shower_detach_track_stem,
            shower_ghost_member_drop=shower_ghost_member_drop,
            shower_ghost_overlap_frac=shower_ghost_overlap_frac,
            shower_ghost_dqdx_ratio=shower_ghost_dqdx_ratio,
            shower_ghost_min_len=shower_ghost_min_len,
            kine_charge_dedup=kine_charge_dedup,
            kine_charge_rebuild=kine_charge_rebuild,
            kine_charge_track_ctx=kine_charge_track_ctx,
            kine_mass_rules=kine_mass_rules,
            kine_hadronic_dqdx=kine_hadronic_dqdx,
            kine_long_muon_mode=kine_long_muon_mode,
            kine_long_muon_ratio_lo=kine_long_muon_ratio_lo,
            kine_long_muon_ratio_hi=kine_long_muon_ratio_hi,
            kine_mainvtx_used_guard=kine_mainvtx_used_guard,
            shower_hadronic_tag=shower_hadronic_tag,
            shower_hadronic_min_len=shower_hadronic_min_len,
            shower_hadronic_scan_len=shower_hadronic_scan_len,
            shower_hadronic_bin=shower_hadronic_bin,
            shower_hadronic_r_cyl=shower_hadronic_r_cyl,
            shower_hadronic_r_core=shower_hadronic_r_core,
            shower_hadronic_growth_max=shower_hadronic_growth_max,
            shower_hadronic_growth_bragg=shower_hadronic_growth_bragg,
            shower_hadronic_bragg_ratio=shower_hadronic_bragg_ratio,
            shower_hadronic_stem_ratio=shower_hadronic_stem_ratio,
            kine_count_orphan_tracks=kine_count_orphan_tracks,
            kine_orphan_track_min=kine_orphan_track_min,
            straight_cont_cross_cluster=straight_cont_cross_cluster,
            sccc_bridge_body=sccc_bridge_body,
            sccc_max_gap=sccc_max_gap,
            sccc_kink_max=sccc_kink_max,
            sccc_gap_aligned=sccc_gap_aligned,
            sccc_kink_tight=sccc_kink_tight,
            single_muon_proton_chain_veto=single_muon_proton_chain_veto,
            single_muon_long_muon_claim=single_muon_long_muon_claim,
            pid_flag_reconcile=pid_flag_reconcile,
            other_seg_empty_2d_guard=other_seg_empty_2d_guard,
            long_muon_stub_bridge=long_muon_stub_bridge,
            two_end_break=two_end_break,
            teb_min_len=teb_min_len,
            teb_min_arm=teb_min_arm,
            teb_min_arm_pts=teb_min_arm_pts,
            teb_stub_max=teb_stub_max,
            teb_accept_range=teb_accept_range,
            teb_rise_r1=teb_rise_r1,
            teb_rise_r2=teb_rise_r2,
            teb_abs_end_min=teb_abs_end_min,
            teb_dip_floor=teb_dip_floor,
            teb_score_cap_r1=teb_score_cap_r1,
            teb_score_cap_r2=teb_score_cap_r2,
            teb_turn_angle=teb_turn_angle,
            teb_turn_baseline=teb_turn_baseline,
            teb_turn_skirt=teb_turn_skirt,
            teb_turn_min_arm_frac=teb_turn_min_arm_frac,
            teb_second_max=teb_second_max,
            teb_chain_topology=teb_chain_topology,
            teb_r3_turn=teb_r3_turn,
            teb_r3_hot=teb_r3_hot,
            teb_bragg_veto_turn=teb_bragg_veto_turn,
            kink_walk_dqdx_stop=kink_walk_dqdx_stop,
            kink_break_protect=kink_break_protect,
            kink_dqdx_hot_ratio=kink_dqdx_hot_ratio,
            fit_blob_coverage=fit_blob_coverage,
            fit_blob_coverage_defer=fit_blob_coverage_defer,
            vertex_kink_snap=vertex_kink_snap,
            vks_radius=vks_radius,
            vks_min_dis=vks_min_dis,
            vks_angle=vks_angle,
            vks_margin=vks_margin,
            vks_collinear=vks_collinear,
            vks_skirt=vks_skirt,
            vks_baseline=vks_baseline,
            vks_min_arm=vks_min_arm,
            vks_fit_miss=vks_fit_miss,
            vks_hot_ratio=vks_hot_ratio,
            vks_carry_prong=vks_carry_prong,
            vertex_junction_snap=vertex_junction_snap,
            vjs_radius=vjs_radius,
            vjs_min_arm=vjs_min_arm,
            vjs_min_prongs=vjs_min_prongs,
            vjs_collinear=vjs_collinear,
            vjs_fit_margin=vjs_fit_margin,
            vjs_fit_rms=vjs_fit_rms,
            vjs_override_kink_snap=vjs_override_kink_snap,
            vjs_min_move=vjs_min_move,
            esva_ignore_empty_2d=esva_ignore_empty_2d,
            main_vertex_graph_audit=main_vertex_graph_audit,
            mvga_radius=mvga_radius,
            mvga_dup_tol=mvga_dup_tol,
            mvga_dup_frac=mvga_dup_frac,
            mvga_dup_angle=mvga_dup_angle,
            mvga_bridge_mip=mvga_bridge_mip,
            mvga_reconnect=mvga_reconnect,
            mvga_stub=mvga_stub,
            mvga_stub_pts=mvga_stub_pts,
            mvga_reseat_angle=mvga_reseat_angle,
            mvga_satellite=mvga_satellite,
            mvga_interposed=mvga_interposed,
            mvga_interposed_angle=mvga_interposed_angle,
            mvga_interposed_len=mvga_interposed_len,
            mvga_sat_dup_frac=mvga_sat_dup_frac,
            mvga_interposed_deg1=mvga_interposed_deg1,
            mvga_splice_straighten=mvga_splice_straighten,
            mvga_approach_collapse=mvga_approach_collapse,
            mvga_straighten_radius=mvga_straighten_radius,
            mvga_op1_radius=mvga_op1_radius,
            mvga_op1_dup_frac=mvga_op1_dup_frac,
            mvga_op1_post=mvga_op1_post,
            mvga_carry_max=mvga_carry_max,
            swap_orphan_dup_audit=swap_orphan_dup_audit,
            mvga_proj_dup_frac=mvga_proj_dup_frac,
            mvga_proj_dqdx_ratio=mvga_proj_dqdx_ratio,
            mvga_proj_angle=mvga_proj_angle,
            mvga_ac_veto_radius=mvga_ac_veto_radius,
            mvga_ac_chord_max=mvga_ac_chord_max,
            mvga_ac_no_cascade=mvga_ac_no_cascade,
            mvga_passthru=mvga_passthru,
            mvga_passthru_tol=mvga_passthru_tol,
            mvga_interposed_fallback=mvga_interposed_fallback,
            mvga_interposed_fallback_min_angle=mvga_interposed_fallback_min_angle,
            mvga_dup_starved_asym=mvga_dup_starved_asym,
            mvga_dup_starved_mip=mvga_dup_starved_mip,
            mvga_dup_starved_span=mvga_dup_starved_span,
            dl_vtx_swap_guard=dl_vtx_swap_guard,
            dl_vtx_cloud_no_exclusion=dl_vtx_cloud_no_exclusion,
            dl_vtx_topo_weight=dl_vtx_topo_weight,
            dl_vtx_topo_center=dl_vtx_topo_center,
            main_vertex_swap_apply=main_vertex_swap_apply,
            rough_path_probe=rough_path_probe,
            steiner_gap_penalty=steiner_gap_penalty,
            sgp_dead_alpha=sgp_dead_alpha,
            sgp_min_edge=sgp_min_edge,
            sgp_sample_step=sgp_sample_step,
            sgp_point_radius=sgp_point_radius,
            sgp_edge_probe=sgp_edge_probe, vertex_scoreboard=vertex_scoreboard, dl_vtx_harvest=dl_vtx_harvest,
            sgp_weak_scale=sgp_weak_scale,
            sgp_weak_qref=sgp_weak_qref,
            sgp_max_sep=sgp_max_sep,
            break_seg_orient=break_seg_orient,
            mvfit_robust=mvfit_robust,
            mvfit_main_only=mvfit_main_only,
            mvfit_min_len=mvfit_min_len,
            mvfit_rin_margin=mvfit_rin_margin,
            mvfit_rout_frac=mvfit_rout_frac,
            mvfit_rout_min=mvfit_rout_min,
            mvfit_rout_max=mvfit_rout_max,
            mvfit_angle=mvfit_angle,
            mvfit_min_pts=mvfit_min_pts,
            mvfit_min_aniso=mvfit_min_aniso,
            mvfit_prior_range=mvfit_prior_range,
            other_seg_keep_isolated=other_seg_keep_isolated,
            other_seg_keep_isolated_min_points=other_seg_keep_isolated_min_points,
            other_seg_keep_isolated_min_length=other_seg_keep_isolated_min_length,
            other_seg_keep_isolated_min_nnf=other_seg_keep_isolated_min_nnf,
            other_seg_keep_isolated_len_admit=other_seg_keep_isolated_len_admit,
            other_seg_uncover_3d=other_seg_uncover_3d,
            iso_snap_min_dir_mag=iso_snap_min_dir_mag,
            shower_absorb_unreachable_main=shower_absorb_unreachable_main,
            assoc_full_recluster=assoc_full_recluster,
            assoc_reassign_orphans=assoc_reassign_orphans,
            assoc_clear_on_merge=assoc_clear_on_merge,
            es3_stub_guard=es3_stub_guard,
            es3sg_stub_max=es3sg_stub_max,
            es3sg_len_ratio=es3sg_len_ratio,
            es3sg_ang3_min=es3sg_ang3_min,
            es3sg_ang_ratio=es3sg_ang_ratio,
            es3sg_require_terminal=es3sg_require_terminal,
            muon_dqdx_curve=muon_dqdx_curve,
            sp_dedx_use_recomb_model=sp_dedx_use_recomb_model,
            sp_mean_dedx_cut=sp_mean_dedx_cut,
            dl_vtx_cut=dl_vtx_cut),
        // NuMu / nue BDT scorers (UbooneNumuBDTScorer / UbooneNueBDTScorer,
        // geometry-free TaggerInfo consumers).  The weights are the
        // uBooNE-TRAINED XMLs from wire-cell-data uboone/weights/ -- the same
        // 35 files the qlport gate job books -- so on SBND the scores are
        // UNCALIBRATED (availability + relative ranking only; SBND retraining
        // is docs/pr/2 gap G1).  Must run after tagger_check_neutrino, nue
        // after numu.  Only compiled in when named in pipeline_names.
        local bdt_weights_dir = 'uboone/weights',
        numu_bdt_scorer: cm.numu_bdt_scorer(
            numu1_weights_xml=     bdt_weights_dir + '/numu_tagger1.weights.xml',
            numu2_weights_xml=     bdt_weights_dir + '/numu_tagger2.weights.xml',
            numu3_weights_xml=     bdt_weights_dir + '/numu_tagger3.weights.xml',
            cosmict10_weights_xml= bdt_weights_dir + '/cos_tagger_10.weights.xml',
            numu_xgboost_xml=      bdt_weights_dir + '/numu_scalars_scores_0923.xml'),
        nue_bdt_scorer: cm.nue_bdt_scorer(
            mipid_weights_xml=       bdt_weights_dir + '/mipid_BDT.weights.xml',
            gap_weights_xml=         bdt_weights_dir + '/gap_BDT.weights.xml',
            hol_lol_weights_xml=     bdt_weights_dir + '/hol_lol_BDT.weights.xml',
            cme_anc_weights_xml=     bdt_weights_dir + '/cme_anc_BDT.weights.xml',
            mgo_mgt_weights_xml=     bdt_weights_dir + '/mgo_mgt_BDT.weights.xml',
            br1_weights_xml=         bdt_weights_dir + '/br1_BDT.weights.xml',
            br3_weights_xml=         bdt_weights_dir + '/br3_BDT.weights.xml',
            br3_3_weights_xml=       bdt_weights_dir + '/br3_3_BDT.weights.xml',
            br3_5_weights_xml=       bdt_weights_dir + '/br3_5_BDT.weights.xml',
            br3_6_weights_xml=       bdt_weights_dir + '/br3_6_BDT.weights.xml',
            stemdir_br2_weights_xml= bdt_weights_dir + '/stem_dir_br2_BDT.weights.xml',
            trimuon_weights_xml=     bdt_weights_dir + '/stl_lem_brm_BDT.weights.xml',
            br4_tro_weights_xml=     bdt_weights_dir + '/br4_tro_BDT.weights.xml',
            mipquality_weights_xml=  bdt_weights_dir + '/mipquality_BDT.weights.xml',
            pio_1_weights_xml=       bdt_weights_dir + '/pio_1_BDT.weights.xml',
            pio_2_weights_xml=       bdt_weights_dir + '/pio_2_BDT.weights.xml',
            stw_spt_weights_xml=     bdt_weights_dir + '/stw_spt_BDT.weights.xml',
            vis_1_weights_xml=       bdt_weights_dir + '/vis_1_BDT.weights.xml',
            vis_2_weights_xml=       bdt_weights_dir + '/vis_2_BDT.weights.xml',
            stw_2_weights_xml=       bdt_weights_dir + '/stw_2_BDT.weights.xml',
            stw_3_weights_xml=       bdt_weights_dir + '/stw_3_BDT.weights.xml',
            stw_4_weights_xml=       bdt_weights_dir + '/stw_4_BDT.weights.xml',
            sig_1_weights_xml=       bdt_weights_dir + '/sig_1_BDT.weights.xml',
            sig_2_weights_xml=       bdt_weights_dir + '/sig_2_BDT.weights.xml',
            lol_1_weights_xml=       bdt_weights_dir + '/lol_1_BDT.weights.xml',
            lol_2_weights_xml=       bdt_weights_dir + '/lol_2_BDT.weights.xml',
            tro_1_weights_xml=       bdt_weights_dir + '/tro_1_BDT.weights.xml',
            tro_2_weights_xml=       bdt_weights_dir + '/tro_2_BDT.weights.xml',
            tro_4_weights_xml=       bdt_weights_dir + '/tro_4_BDT.weights.xml',
            tro_5_weights_xml=       bdt_weights_dir + '/tro_5_BDT.weights.xml',
            nue_xgboost_xml=         bdt_weights_dir + '/XGB_nue_seed2_0923.xml'),
        // PR-stage Magnify-tracking ROOT dump (docs/pr/3): fork of the uBooNE
        // writer reading the unnamed TrackFitting slot + PRGraph filled by
        // tagger_check_neutrino, with the two-TPC channel convention and
        // per-point APA from PR::Fit::paf.  Only active when named in
        // pipeline_names (the WireCellRoot plugin must be loaded by the job).
        local tracking_pr_root = (if output_dir == '' then '' else output_dir + '/') + 'tracking-pr.root',
        tracking_visitor: {
            type: 'SbndPrMagnifyTrackingVisitor',
            name: 'pr',
            data: {
                grouping: 'live',
                output_filename: tracking_pr_root,
                runNo: runNo,
                subRunNo: subRunNo,
                eventNo: eventNo,
                anodes: [wc.tn(a) for a in anodes],
                detector_volumes: wc.tn(dv),
                dQdx_scale: 0.1,
                dQdx_offset: -1000.0,
                flag_skip_vertex: false,
                // Readout length in ticks; only clamps T_bad_ch time ranges.
                nticks: 3427,
            },
        },
        // T_tagger/T_kine writer (UbooneTaggerOutputVisitor, reused as-is: it
        // is a pure TaggerInfo/KineInfo dump with no geometry).  Opens the
        // tracking file in UPDATE mode, so it must be named AFTER
        // tracking_visitor (and after both BDT scorers so the scores are set).
        // doc pr/36 sec 10.8 (F7): neutrino_type_bitmask books the
        // neutrino_type/I branch; key omitted when off => schema-identical.
        // doc pr/94 Phase 1: nu_per_bundle books the per-bundle identity +
        // per-activity cosmic-flag branches; key omitted when off =>
        // schema-identical.  Plumbing only -- nothing populates them yet.
        tagger_output: cm.tagger_output(output_filename=tracking_pr_root,
                                        neutrino_type_bitmask=neutrino_type_bitmask,
                                        nu_per_bundle=nu_per_bundle),
        // PR event-display calib dump (docs/pr/26): ONE self-contained JSON per
        // event carrying the PR-graph segments as polylines, the associated
        // track/shower points, the Steiner skeleton with its terminal flag, the
        // fitted 2-D charge per (apa,face,plane) and the dead regions -- i.e.
        // the union of what the Bee PR layers and the Magnify tracking file
        // carry, in one file the Bokeh viewer can read with the stdlib alone.
        // Read-only; mutates nothing.  Only active when named in
        // pipeline_names => compiled config byte-identical otherwise.
        // Must run AFTER tagger_check_neutrino (it reads the TrackFitting slot
        // and the PR graph that stage fills).
        pr_display: {
            type: 'PrDisplayDump',
            name: 'pr',
            data: {
                grouping: 'live',
                output_filename: (if output_dir == '' then '' else output_dir + '/')
                                 + 'calib-pr-evt' + std.toString(eventNo) + '.json',
                runNo: runNo,
                subRunNo: subRunNo,
                eventNo: eventNo,
                anodes: [wc.tn(a) for a in anodes],
                detector_volumes: wc.tn(dv),
                // Same convention as the Bee track_fit layer and the Magnify
                // writer; the dump stores RAW dQ and records these so the
                // viewer can reproduce the Bee colouring if it wants to.
                dQdx_scale: 0.1,
                dQdx_offset: -1000.0,
                nticks: 3427,
                // doc pr/45.  C++ default false; mirrors the Bee shower_track
                // layer knob.  Key omitted when off => byte-identical.
                [if pseudo_shower_track_paint
                 then 'pseudo_shower_track_paint']: true,
            },
        },
    },
    local cm_pipeline = [cm_by_name[n] for n in pipeline_names],
    // The taggers' configs only name the recombination/particle-dataset
    // components; emit them (and the caller's LinterpFunctions etc. via
    // extra_uses) when a tagger is in the pipeline.
    local tagger_uses = (if std.member(pipeline_names, 'tagger_check_stm')
                         || std.member(pipeline_names, 'tagger_check_neutrino')
                         then [sbnd_recomb] + extra_uses else [])
                        + (if std.member(pipeline_names, 'tagger_check_tgm')
                           || std.member(pipeline_names, 'tagger_check_fc')
                           || (stm_consistent_fv
                               && std.member(pipeline_names, 'tagger_check_stm'))
                           // doc pr/36 sec 10 (F1): the neutrino tagger also
                           // names sbnd_pr_fv when its consistent-FV knob is
                           // on.  Redundant in the default pipeline (TGM/FC
                           // already pull it in) => compiled config unchanged
                           // there; load-bearing only for a reduced pipeline.
                           || ((neutrino_consistent_fv || cosmic_consistent_fv || nue_sp_consistent_fv)
                               && std.member(pipeline_names, 'tagger_check_neutrino'))
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
                    // C++ default false.  Append the PR-graph vertex fit points
                    // (real_cluster_id=-1) like the prototype's
                    // fill_skeleton_info_magnify rows (docs/pr/3).  Key present
                    // only when the PR visitor is in the pipeline => default
                    // production compiled config stays byte-identical.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'include_vertex_points']: true,
                    // require_pr_graph: this layer is PR output.  Without it,
                    // an event where TaggerCheckNeutrino selects no candidate
                    // gets the WHOLE clustering dumped here in raw (un-T0-cor)
                    // coordinates -- see docs/pr/3 sec. 9.  C++ default false
                    // (the uBooNE steiner-bound sets rely on that fallback);
                    // key present only when the PR visitor is in the pipeline
                    // => default compiled config byte-identical.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'require_pr_graph']: true,
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
                    // C++ default false.  Prototype per-particle sub-clustering:
                    // non-shower segments get real_cluster_id = cluster*1000+seg
                    // (NeutrinoID::fill_point_info) so Bee can color by particle
                    // (docs/pr/3).  Key present only when the PR visitor is in
                    // the pipeline => default compiled config byte-identical.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'particle_ids']: true,
                    // doc pr/45.  C++ default false.  Muon-typed (+-13)
                    // pseudo-showers paint as track (q=0), matching the PF
                    // tree's "mu-" verdict from the same cached type.  Key
                    // omitted when off => byte-identical.
                    [if pseudo_shower_track_paint
                     then 'pseudo_shower_track_paint']: true,
                    // require_pr_graph: this layer is PR output.  Without it,
                    // an event where TaggerCheckNeutrino selects no candidate
                    // gets the WHOLE clustering dumped here in raw (un-T0-cor)
                    // coordinates -- see docs/pr/3 sec. 9.  C++ default false
                    // (the uBooNE steiner-bound sets rely on that fallback);
                    // key present only when the PR visitor is in the pipeline
                    // => default compiled config byte-identical.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'require_pr_graph']: true,
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
                    // require_pr_graph: this layer is PR output.  Without it,
                    // an event where TaggerCheckNeutrino selects no candidate
                    // gets the WHOLE clustering dumped here in raw (un-T0-cor)
                    // coordinates -- see docs/pr/3 sec. 9.  C++ default false
                    // (the uBooNE steiner-bound sets rely on that fallback);
                    // key present only when the PR visitor is in the pipeline
                    // => default compiled config byte-identical.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'require_pr_graph']: true,
                },
            ]
            // STM fit-trajectory layer, dumped after TaggerCheckSTM runs (doc
            // sbnd_xin/docs/40).  q = dQ*0.1 - 1000, same convention as the
            // PR track_fit layer.  Entry present only when save_stm_fit is on
            // => compiled config byte-identical otherwise.  Name must avoid
            // the substring '-track' (bee3 models.py filters such files).
            + (if save_stm_fit then [{
                name: 'stm_fit',
                visitor: 'TaggerCheckSTM:pr',
                grouping: 'live',
                detector: 'sbnd',
                algorithm: 'stm_fit',
                pcname: 'stm_fit',
                coords: ['x', 'y', 'z'],
                individual: false,
                dQdx_scale: 0.1,
                dQdx_offset: -1000.0,
            }] else []),
            // Particle-flow Bee output ("mc" jsTree JSON), emitted once after
            // TaggerCheckNeutrino runs; inert when the visitor is not in the pipeline.
            bee_pf: [
                {
                    name: 'mc',
                    visitor: 'TaggerCheckNeutrino:pr',
                    grouping: 'live',
                    // C++ defaults false/0 (legacy).  Prototype mc.json parity
                    // (docs/pr/3): TDatabasePDG-style names + integer MeV, and
                    // the WCReader::KeepMC display floors (5 MeV em, 10 MeV
                    // nucleon).  Keys present only when the PR visitor is in
                    // the pipeline => default compiled config byte-identical.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'prototype_names']: true,
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'em_ke_min']: 5 * wc.MeV,
                    // doc pr/38: nucleon floor lowered from the prototype's
                    // 10 MeV (WCReader::KeepMC) to 3 MeV, owner decision
                    // 2026-08-05 -- sub-10-MeV protons attached at the
                    // neutrino vertex (18255-447477 3.2 MeV, 18255-52657
                    // 8 MeV) must show in the particle flow.
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'np_ke_min']: 3 * wc.MeV,
                    // doc pr/34 §10 port-fidelity knobs.  C++ defaults false;
                    // key omitted when off => byte-identical pre-knob config.
                    [if pf_track_main_cluster_only then 'pf_track_main_cluster_only']: true,
                    [if pf_track_bridged_clusters then 'pf_track_bridged_clusters']: true,   // doc pr/40 round 9 B2
                    [if pf_shower_vertex_barrier then 'pf_shower_vertex_barrier']: true,
                    [if pf_shower_parent_precedence then 'pf_shower_parent_precedence']: true,
                    [if pf_pi0_node_per_id then 'pf_pi0_node_per_id']: true,
                    [if pf_pdg_name_prototype_fallback then 'pf_pdg_name_prototype_fallback']: true,
                    // doc pr/38 Round 4.  C++ default false; key omitted when
                    // off => byte-identical pre-knob config.
                    [if pf_orphan_track_parentage then 'pf_orphan_track_parentage']: true,
                    // doc pr/65 round 3.  C++ default false; key omitted when
                    // off => byte-identical pre-knob config.
                    [if pf_orphan_audit_only then 'pf_orphan_audit_only']: true,
                    // doc pr/84 round 2 (F1/F2).  C++ defaults false (scalars
                    // 3 cm / 8 cm); keys omitted when off/null =>
                    // byte-identical pre-knob config.  Params in cm.
                    [if pf_direct_when_touching then 'pf_direct_when_touching']: true,
                    [if pf_touch_max != null then 'pf_touch_max']: pf_touch_max * wc.cm,
                    [if pf_touch_cross_main then 'pf_touch_cross_main']: true,
                    [if pf_touch_cross_max != null then 'pf_touch_cross_max']: pf_touch_cross_max * wc.cm,
                    [if pf_pseudo_gap_from_main then 'pf_pseudo_gap_from_main']: true,
                    [if pf_unique_node_ids then 'pf_unique_node_ids']: true,
                    [if pf_drop_stray_satellites then 'pf_drop_stray_satellites']: true,
                    // doc pr/93 round 4; params in cm.
                    [if pf_orphan_confident_track then 'pf_orphan_confident_track']: true,
                    [if pf_orphan_track_min_cm != null then 'pf_orphan_track_min']: pf_orphan_track_min_cm * wc.cm,
                    [if pf_track_owns_loose_vertex then 'pf_track_owns_loose_vertex']: true,
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
    // Reco-chain reality config -- ONE place grouping every reality-dependent
    // toggle.  use_sce: run the all-APA clustering + Bee in SCE true space
    // (x_sce) vs the T0-corrected reco scope (x_t0cor).  pos_offset_on: per-TPC
    // transverse (y,z) calibration, data-only (see the pos_offset_a0/a1 comment
    // above).
    // NOTE: sim use_sce set to false to match sbnd_xin's MC chain (Xin runs MC
    // in the T0-corrected reco scope, not SCE true space).  Both realities now
    // cluster in x_t0cor; sim differs only by pos_offset_on=false (MC has no
    // per-TPC transverse misalignment).
    local reco = {
        sim:  { use_sce: false, pos_offset_on: false },
        data: { use_sce: false, pos_offset_on: true },
    }[reality],
    local use_sce = reco.use_sce,
    local pos_offset_on = reco.pos_offset_on,
    // bee_sink (default null): when set to a shared IBeeSink node, all Bee
    // output for this node goes into that single shared zip instead of this
    // node's own bee_zip.  Default null -> own zip (production byte-identical).
    per_face(anode, face=0, dump=true, bee_sink=null)::
        clus_per_face(anode, face=face, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on),
    // trace_bee (default false): per-step Bee layers for merge attribution; see
    // trace_sets above.  Diagnostic only, off => byte-identical compiled config.
    per_apa(anode, dump=true, bee_sink=null, trace_bee=false, save_assoc_id=false, sep_vertex_veto=true, nu_iso_band_guard=true, iso_cathode_guard=false, nu_band_veto=true)::
        clus_per_face(anode, face=0, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                      trace_bee=trace_bee, save_assoc_id=save_assoc_id, sep_vertex_veto=sep_vertex_veto,
                      nu_iso_band_guard=nu_iso_band_guard, iso_cathode_guard=iso_cathode_guard,
                      nu_band_veto=nu_band_veto),
    // Production (LArSoft) entry point used by wcls-img-clus.jsonnet.
    per_volume(anode, face=0, dump=true, bee_sink=null)::
        clus_per_face(anode, face=face, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on),
    all_apa(anodes, dump=true, bee_sink=null, premerged=false, tensor_outname='', save_real_cluster_id=false, save_assoc_cluster_id=false,
            trace_bee=false, real_cluster_id_global=null, cathode_rescue_on=true, cathode_rescue_unmatched=true, adopt_nu_fragments=false,
            save_bundle_main_provenance=false,
            rescue_allow_in_beam_far=true, rescue_geom_first=true,
            rescue_pierce_test=true, rescue_pierce_cut=null,
            rescue_dest_beam_for_new=true, rescue_beam_main_only=true,
            bee_flash_pred_min=null)::
        // Clustering + matching ONLY (all-APA MABC).  The follow-up PR tagger
        // pass (pr() below) and the wclsTensorSetLabeler are wired by the entry
        // configuration, not here -- see the note in clus_all_apa.
        clus_all_apa(anodes, dump=dump,
                     output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                     bee_sink=bee_sink, premerged=premerged, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                     tensor_outname=tensor_outname, save_real_cluster_id=save_real_cluster_id,
                     save_assoc_cluster_id=save_assoc_cluster_id,
                     real_cluster_id_global=real_cluster_id_global,
                     trace_bee=trace_bee, cathode_rescue_on=cathode_rescue_on,
                     cathode_rescue_unmatched=cathode_rescue_unmatched,
                     adopt_nu_fragments=adopt_nu_fragments,
                     save_bundle_main_provenance=save_bundle_main_provenance,
                     rescue_allow_in_beam_far=rescue_allow_in_beam_far,
                     rescue_geom_first=rescue_geom_first,
                     rescue_pierce_test=rescue_pierce_test,
                     rescue_pierce_cut=rescue_pierce_cut,
                     rescue_dest_beam_for_new=rescue_dest_beam_for_new,
                     rescue_beam_main_only=rescue_beam_main_only,
                     bee_flash_pred_min=bee_flash_pred_min,
                     use_sce=use_sce, reality=reality),
    // PR job: input is the reloaded post-QL tarball (see clus_pr above).
    // The TGM/FC and beam-window defaults here mirror clus_pr's -- i.e. the SBND
    // production operating point (see the comment block on clus_pr's arg list for
    // the pre-adoption values to pass for an A/B).
    pr(anodes, dump=true, pipeline_names=[], tensor_outname='',
       trackfitting_config_file='pgrapher/experiment/sbnd/sbnd_track_fitting.json',
       particle_dataset=null, extra_uses=[],
       // DL (SCN) vertex ON by default -- see the clus_pr arg comment.
       dl_weights='uboone/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth',
       // DL re-rank operating point -- see the clus_pr arg comment (doc pr/79).
       dl_vtx_min_accept_score=10.0,
       dl_vtx_top_k=5,
       dl_vtx_rerank=true,   // doc pr/105: see the clus_pr arg comment
       beam_window=[0.2 * wc.us, 2.2 * wc.us],
       tgm_neutrino_candidate=true,
       tgm_chord_charge=true, tgm_chord_mode='path',
       tgm_component_extremes=true, tgm_component_rescue=true,
       tgm_rescue_chord=true, tgm_main_pair=true, tgm_main_pair_mode='real',
       tgm_fv_zmax_margin=5, tgm_fv_zmax_margin_interior=3,
       tgm_fv_x_margin=2.5, tgm_fv_y_margin=3,
       save_stm_fit=false, unmerge_bundle_mode='real',
       // doc pr/34 §10 particle-flow port-fidelity knobs; false = C++ default
       // = OFF, key omitted => byte-identical.  See clus_pr.
       pf_track_main_cluster_only=false,
       pf_track_bridged_clusters=false,   // doc pr/40 round 9 B2
       pf_shower_vertex_barrier=false,
       pf_shower_parent_precedence=false,
       pf_pi0_node_per_id=false,
       pf_pdg_name_prototype_fallback=false,
       // doc pr/38 Round 4; false = C++ default = OFF, key omitted =>
       // byte-identical.  See clus_pr.
       pf_orphan_track_parentage=false,
       // doc pr/65 round 3; false = C++ default = OFF.  See clus_pr.
       pf_orphan_audit_only=false,
       // doc pr/84 round 2 (F1/F2); false/null = C++ defaults = OFF, keys
       // omitted => byte-identical.  Params in cm.  See clus_pr.
       pf_direct_when_touching=false,
       pf_touch_max=null,
       pf_touch_cross_main=false,
       pf_touch_cross_max=null,
       pf_pseudo_gap_from_main=false,
       // doc pr/84 round 3 (G1); false = C++ default = OFF.  See clus_pr.
       pf_unique_node_ids=false,
       // doc pr/92; false = C++ default = OFF.  See clus_pr.
       pf_drop_stray_satellites=false,
       // doc pr/93 round 4; C++ defaults; keys suppressed when off.
       pf_orphan_confident_track=false,
       pf_orphan_track_min_cm=null,
       pf_track_owns_loose_vertex=false,
       // doc pr/20 Part I P2; null = C++ default false = OFF.  See clus_pr.
       restore_demoted_mains=null,
       // doc pr/23 sec 4.2; null = C++ default false = warn-and-skip.  See clus_pr.
       require_provenance=null,
       // doc pr/20 Part I P3; false = C++ default = OFF.  See clus_pr.
       evaluate_demoted_mains=false,
       // doc pr/25, SBND evt 320029; false = C++ default = OFF.  See clus_pr.
       tgm_exempt_demoted_main=false,
       // doc pr/20 Part I P4; false / null = C++ defaults = OFF.  See clus_pr.
       skip_cosmic_companions=false, cosmic_companion_min_length=null,
       // docs/73 sec 12 round 3; SBND PRODUCTION ON 2026-08-17.  See clus_pr.
       nu_fallback_demoted_mains=true,
       // sp_photon_flag: store the single-photon tagger's verdict in
       // TaggerInfo::photon_flag, as prototype NeutrinoID.cxx:271 does.
       // The port ran singlephoton_tagger() and filled its shw_sp_*
       // features but discarded the verdict (docs/pr/26 sec. 8.2).
       // C++ default false = that gap; key omitted when off =>
       // byte-identical.  Only the uBooNE tagger ntuple's photon_flag
       // branch changes when on -- nothing in the chain reads it.
       sp_photon_flag=false,
       mip_dqdx=56000, stm_consistent_fv=true, stm_accept_guards=true,
       stm_proton_muon_guard=true, stm_cathode_guard=true,
       stm_anode_dist_fix=true, stm_second_track_guard=true,
       stm_deficit_guard=true, stm_vertex_kink_guard=true,
       stm_d66_cuts=true, stm_michel_res_cm=6.5, stm_proton_tm_max=1.05,
       stm_proton_b_ks2_max=0.055, stm_proton_c_peak_max=4.1,
       beam_window_only=true, nu_skip_cosmic=true,
       // Bundle-level cosmic veto -- see the clus_pr arg comment
       // (docs/pr/3 sec. 8).  NOT bit-identical when on.
       nu_skip_cosmic_bundle=true,
       // Design-A size guard on the veto -- see the clus_pr arg comment
       // (docs/pr/16 sec. 7).  cm; SBND 15 (owner 2026-08-01).
       nu_skip_cosmic_bundle_min_length=15,
       // prototype-faithful is_dir_weak() reads -- see the clus_pr arg comment.
       dir_weak_use_score=true,
       // PR-chain median-dQ/dx scale + proton direction vote -- see the
       // clus_pr arg comments (docs pr/7 sec 5, pr/8).
       mip_dqdx_median=48000, proton_dir_vote=true,
       // Endpoint-trim retry (doc pr/9 sec 6 F1) -- see the clus_pr arg comment.
       endpoint_trim_retry=true,
       // fit_vertex short-segment exclusion (doc pr/9 sec 11 F3c), cm -- see
       // the clus_pr arg comment.  null would give the legacy include-all fit.
       fit_vertex_min_seg_length=1.0,
       // Cathode kink veto (doc pr/20 Part II B0), both cm.  cathode_kink_xcut
       // suppresses segment_search_kink accept candidates within that distance
       // of x = cathode_x, because the kink finder's para_angle guard is
       // inverted at the cathode and breaks crossing cosmics in two.
       // C++ defaults are 0/0 = OFF; null omits both keys.
       // SBND DEFAULT ON (owner 2026-08-02 after the Part VI Bee scan):
       // 5 cm around x = 0.  NOT bit-identical -- 21 of 1000 mcp1k events move,
       // 10 of them relocating the neutrino vertex (doc pr/20 Part VI sec 1).
       // Set both to null for the legacy kink search.
       cathode_x=0, cathode_kink_xcut=5,
       // Wide-baseline cathode kink accept (doc pr/47 sec 8, O1), the
       // converse of the veto: a real kink AT the crossing never reaches the
       // legacy accept thresholds because the gap/distortion suppresses every
       // local window (52085: 33-38 deg junction reads ~23 deg).  Fires when
       // the skirt-excluded PCA turn angle across the crossing is >= this
       // many degrees.  25 sits mid-gap in the measured bimodal census
       // distribution (bulk p90 8.2 deg, tail >= 36.8 deg; doc pr/47 sec 6).
       // SBND DEFAULT ON (owner 2026-08-07).  NOT bit-identical -- footprint
       // measured in doc pr/47 sec 8.  null = legacy (no cathode accept).
       // Skirt/baseline stay null = C++ defaults 3 cm / 15 cm.
       cathode_wide_kink_angle=25,
       cathode_wide_kink_skirt=null,
       cathode_wide_kink_baseline=null,
       // shower_topo_demote_len (cm, doc pr/25 sec 3): demote any
       // kShowerTopology segment longer than this to a track, so it gets real
       // track PID instead of the hard-coded pdg=11/score=100.  Owner
       // hand-scan 2026-08-03: 10/10 long shower-topology segments on a
       // selected nu-candidate main cluster are tracks, none showers.
       // null = C++ default 0 = OFF = byte-identical.  Ships OFF; the SBND
       // operating point lives in wct-pr-perevt.jsonnet (doc 68).
       shower_topo_demote_len=null,
       // ---- doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs ------------
       // All five default to the pre-pr/30 behaviour, so the compiled JSON is
       // byte-identical until one is set.  See cfg/pgrapher/common/clus.jsonnet
       // for what each does.  fit_exclusion (P1), graph_endpoint_strict (P8) and
       // oov_prototype_parity (F2) turn NEW behaviour on (default off);
       // first_seg_local_pca (P2) and other_seg_relaxed_accept (P4) gate behaviour
       // that is ALREADY production, so null = on = legacy and false restores the
       // prototype's narrower version.
       fit_exclusion=false, graph_endpoint_strict=false, graph_endpoint_tol=null,
       oov_prototype_parity=false, first_seg_local_pca=null, other_seg_relaxed_accept=null,
       // shower_topo_proto_dir (doc pr/31 sec 11, F2 was P2): true skips the
       // stage-3 segment_determine_shower_direction call, leaving the topology
       // shower with the direction segment_is_shower_topology set -- the
       // prototype's state.  C++ default false = today's path = byte-identical.
       shower_topo_proto_dir=false,
       // doc pr/32 sec 11: the four stage-4 vertex-ID port fixes.  Ship OFF =
       // byte-identical; the SBND operating point lives in wct-pr-perevt.jsonnet.
       vertex_dir_use_fit_point=false, shower_traj_recheck_parity=false,
       main_vertex_require_descriptor=false, main_vertex_candidate_flag=false,
       // doc pr/31 sec 12: the sec 10.12 topology/PID/direction port fixes.
       // Ship OFF = byte-identical; the SBND operating point lives in
       // wct-pr-perevt.jsonnet.  F7 (examine_showers_vertex_by_index) is
       // deliberately dormant pending pr/30 F4.
       cont_muon_dir3_30cm=false, track_comp_empty_abstain=false,
       shower_topo_reset=false, reclass_preserve_4mom=false,
       dir_track_median_local=false, examine_showers_vertex_by_index=false,
       // Steiner TERMINAL filter fidelity (doc pr/29 D1 and D12).  Both OFF =
       // the historical toolkit behaviour = byte-identical (keys omitted).
       // Ships OFF; the SBND operating point lives in wct-pr-perevt.jsonnet.
       steiner_terminal_wire_tol=0,
       steiner_terminal_adjacent_slice=false,
       // Steiner EDGE-WEIGHT charge fidelity (doc pr/29 D2).  OFF = historical
       // toolkit behaviour = byte-identical (key omitted).  Ships OFF; the SBND
       // operating point lives in wct-pr-perevt.jsonnet.
       steiner_edge_charge_forward_dead_mix=false,
       // Isochronous first-segment endpoint finding (doc pr/24 round 2, SBND
       // evt 271851).  false/nulls = C++ defaults = OFF = byte-identical.
       iso_endpoint=false,
       iso_endpoint_min_length=null,
       iso_endpoint_max_xext=null,
       iso_endpoint_xext_frac=null,
       iso_endpoint_xext_quantile=null,
       iso_endpoint_tube_radius=null,
       iso_endpoint_min_aspect=null,
       // examine_vertices_3 extension-retraction guard (doc pr/24 round 5).
       // false/null = C++ defaults = OFF = byte-identical.
       v3_extension_guard=false,
       v3_extension_min_gain=null,
       // doc pr/67: log-only trajectory-coverage probe + the counterfactual
       // override for find_proto_vertex's hardcoded main-cluster
       // branch-search round budget.  false/null = C++ defaults = OFF =
       // byte-identical.
       traj_cover_probe=false,
       pr_find_other_rounds=null,
       // protect_bundle (doc pr/23): PR-stage overclustering protection.
       // Named in the production pipeline_names by DEFAULT since the sec 9
       // flip (owner 2026-08-02); inert when dropped from the list.  The
       // cathode re-join values are the SBND operating point (INTERNAL units,
       // unlike cathode_kink_xcut above): re-unite graph components whose
       // closest points both sit within 5 cm of the cathode, within 8 cm in
       // 3D and 4 cm transversely, so the stage never re-splits a cathode
       // crosser that cathode_connect/B0 preserved (doc pr/20).  nulls =
       // prototype-faithful (re-join pass disabled).
       protect_graph_name=null,
       protect_skip_convicted=null,
       protect_open_convicted_bundles=null,   // doc pr/94 round 3; null => C++ default false
       protect_cathode_x=0,
       protect_cathode_rejoin_xcut=5 * wc.cm,
       protect_cathode_rejoin_dyz=4 * wc.cm,
       protect_cathode_rejoin_dis=8 * wc.cm,
       // Direction-agreement fallback for a dyz-only failure (doc pr/25,
       // SBND evt 489327): DESIGNED, NOT YET the SBND default (owner has not
       // signed off on turning it on -- escalation rule 1, changes which
       // cathode crossers get re-joined).  nulls => C++ default 0 = fallback
       // disabled => byte-identical to the block above.  Candidate operating
       // point once approved: perp=3cm, angle=20.0, dir_radius=15cm,
       // dir_npts=20 (margins from the pr25_rejoin_census.py 572-event scan:
       // 5 genuine crossers at angle<=7.1 deg / perp<=2.68 cm vs the one
       // same-side junk stub at angle 89.2 deg / perp 6.11 cm, which also
       // fails dir_npts independently).
       protect_cathode_rejoin_perp=null,
       protect_cathode_rejoin_angle=null,
       protect_cathode_rejoin_dir_radius=null,
       protect_cathode_rejoin_dir_npts=null,
       // Detector-extent literals re-anchored to SBND (docs/pr/2 sec 2e(iv)) --
       // see the clus_pr arg comments.  cosmic_y_*: uBooNE's "reaches the top"
       // offsets carried from its y=+117 cm top face to SBND's +200 cm.
       // vertex_z_prior_scale: upstream-z prior scaled by the detector length
       // (200 x 501/1037 -> 100 cm).  null on any of them restores uBooNE.
       cosmic_y_top_main=sbnd_y_top - 17,
       cosmic_y_top_strict=sbnd_y_top - 15,
       cosmic_y_top_loose=sbnd_y_top - 37,
       cosmic_y_small_piece=sbnd_y_top - 67,
       vertex_z_prior_scale=100.0,
       // SSM beam-line references -- null = uBooNE defaults, see the clus_pr
       // arg comment.  No SBND value exists yet (docs/pr/2 sec 2e(i)).
       ssm_target_dir=null, ssm_absorber_dir=null,
       // Charge -> kinetic-energy calibration (docs/pr/2 sec 2e(iii)).  The
       // three recombination survival factors carry the docs/pr/10 SBND
       // values (table-integrated ratio transfer); fudge factors, plane
       // weights, asymmetry switch and W-value remain uBooNE (null) -- see
       // the clus_pr arg comment.
       kine_fudge_factor=null, kine_recom_factor=0.87,
       kine_shower_fudge_factor=null, kine_shower_recom_factor=0.58,
       kine_proton_recom_factor=0.51, kine_plane_weights=null,
       kine_plane_asym_switch=null, kine_w_value=null,
       // doc pr/35 sec 10.2 (F1): live shower PDG at fill_kine_tree.
       // C++ default false; key omitted when off => byte-identical.
       kine_shower_pdg_live=false,
       // doc pr/36 sec 10 tagger-stage knobs -- see the clus_pr arg
       // comments.  All false = keys omitted = byte-identical pre-knob
       // config (and, for neutrino_type_bitmask, an identical T_tagger
       // schema).
       neutrino_consistent_fv=false,
       cosmic_consistent_fv=false,
       nue_sp_consistent_fv=false,
       sp_sce_correction=false,
       tagger_ordered_segment_sets=false,
       stem_endpoint_wcpt_parity=false,
       broken_muon_cluster_id_count=false,
       neutrino_type_bitmask=false,
       // doc pr/94 Phase 1: per-bundle identity + per-activity cosmic-flag
       // T_tagger branches (tagger_output only, plumbing only).  C++
       // default false; key omitted when off => byte-identical.
       nu_per_bundle=false,
       nu_per_bundle_min_length=null,   // doc pr/94 Phase 5b; cm; null => C++ default 0 (no floor)
       nu_selected_as_main=false,       // doc pr/94 round 3; C++ default false; key omitted when off
       nu_selected_as_main_snapshot_all=false,  // doc 75; C++ default false; key omitted when off
       // doc pr/33 sec 10 EM-shower-clustering knobs -- see the clus_pr arg
       // comments.  All false = keys omitted = byte-identical pre-knob
       // config.
       daughter_count_proto_main_vertex=false,
       daughter_count_proto_examine_showers=false,
       shower_pdg_from_start_segment=false,
       shower_pdg_from_shower_type=false,
       shower_pdg_exact_muon_test=false,
       pi0_id_shared_allocator=false,
       shower_flag_pdg_electron=false,
       shower_less_id_tiebreak=false,
       // doc pr/39: exclude a shower's own start vertex from the end_point
       // farthest-vertex search.  Ships OFF pending owner gate review.
       shower_endpoint_exclude_start_vertex=false,
       shower_endpoint_skip_orphan_vtx=false,
       // doc pr/91 round 3: flood-fill frontier = visited, not merely
       // present in the view.  Ships OFF pending owner gate review.
       shower_walk_visited_parity=false,
       // doc sbnd_xin/docs/pr/40 -- track (proton/pion/muon) mis-identified
       // as electron.  All default false = legacy = byte-identical.
       track_pid_persist_dqdx=false,
       shower_reclass_dqdx_guard=false,
       shower_topo_dqdx_guard=false,
       // doc sbnd_xin/docs/pr/40 round 2 -- two follow-on defects from the pr/40
       // round.  All default false = legacy = byte-identical.
       reclass_never_computed_ke_floor=false,
       track_pid_persist_4mom=false,
       shower_proton_daughter_pion=false,
       // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on defects from round
       // 2/3's F5.  All default false = legacy = byte-identical.  See the
       // clus_pr arg comment.
       shower_proton_daughter_pion_dissolve=false,
       muon_multi_proton_pion=false,
       // doc sbnd_xin/docs/pr/40 round 5 -- muon mis-identified as electron,
       // three independent mechanisms.  All default false = legacy =
       // byte-identical.  See the clus_pr arg comment.
       track_pid_persist_dqdx_electron_guard=false,
       shower_connect_main_vertex_straight_guard=false,
       shower_traj_straight_guard=false,
       // doc sbnd_xin/docs/pr/40 round 6 -- boundary-level fixes.  All
       // default false = legacy = byte-identical.  See the clus_pr arg
       // comment.
       shower_absorb_track_guard=false,
       shower_connect_protected_pion_guard=false,
       michel_stem_muon_rescue=false,
       // doc sbnd_xin/docs/pr/74 round 2 -- P1 + P2 (see clus_pr above).
       shower_in_cascade_guard=false,
       shower_in_max_len=null,
       shower_in_mip_hi=null,
       // doc sbnd_xin/docs/pr/40 round 9 -- see clus_pr above.
       shower_connect_from_vertices_straight_guard=false,
       shower_connect_start_seg_straight_guard=false,
       examine_direction_dirsign_shower_in_guard=false,
       daughter_shower_angle_reclass_straight_guard=false,
       shower_topo_reexam_straight_guard=false,
       sfv_kink_max=null,
       shower_nv_bridge_track=false,
       shower_nv_bridge_max_gap=null,
       shower_nv_main_pi_init=false,   // doc pr/97 D1; false => legacy indeterminate main_pi read
       // doc pr/92; false/null = C++ defaults = OFF.  See clus_pr.
       kine_drop_stray_satellites=false,
       kine_sat_min_energy=null,
       kine_sat_prox_max=null,
       kine_sat_angle_bad=null,
       kine_sat_angle_main=null,
       kine_sat_far_dis=null,
       kine_sat_axis_dis_cut=null,
       kine_sat_cont_kink=null,
       kine_sat_track_max_nseg=null,
       kine_sat_em_far_dis=null,
       michel_stem_michel_check=false,
       michel_stem_max_far_len=null,
       shower_stem_backfill=false,
       stem_backfill_max_len=null,
       stem_backfill_mip_lo=null,
       stem_backfill_mip_hi=null,
       stem_backfill_min_shower_len=null,
       shower_conn3_unreachable=false,
       conn3_unreachable_min_len=null,
       // doc pr/84 round 2 (F3); null = C++ default 0 = OFF.  See clus_pr.
       conn3_stitch_max=null,
       // doc pr/84 round 3 (S1); false = C++ default = OFF.  See clus_pr.
       shower_dedup_start_seg=false,
       shower_traj_michel_stem=false,
       michel_stem_traj_min_len=null,
       michel_stem_traj_max_len=null,
       michel_stem_traj_mip_lo=null,
       michel_stem_traj_max_far_len=null,
       michel_stem_traj_min_kink_deg=null,
       // doc pr/44; false = C++ default = OFF, key omitted => byte-identical.
       shower_long_muon_keep_type=false,
       // doc pr/40 round 10; false = C++ default = OFF, key omitted => byte-identical.
       shower_bragg_protect_start_segment=false,
       // doc pr/93 round 3 -- C++ defaults false; key-suppressed when off.
       shower_reclass_case_b_dqdx_guard=false,
       shower_accept_pid_guard=false,
       shower_pid_guard_min_len=null,
       shower_vote_track_pid_counts=false,
       shower_cone_absorb_guard=false,
       // doc pr/93 round 4; C++ defaults; keys suppressed when off.
       shower_detach_track_stem=false,
       shower_ghost_member_drop=false,
       shower_ghost_overlap_frac=null,
       shower_ghost_dqdx_ratio=null,
       shower_ghost_min_len=null,
       // doc pr/99 round 3; C++ defaults; keys suppressed when off.
       kine_charge_dedup=false,
       kine_charge_rebuild=false,
       // doc pr/101 Enu accounting round; C++ defaults; keys suppressed when off.
       kine_charge_track_ctx=false,
       kine_mass_rules=false,
       kine_hadronic_dqdx=false,
       kine_long_muon_mode=null,
       kine_long_muon_ratio_lo=null,
       kine_long_muon_ratio_hi=null,
       kine_mainvtx_used_guard=false,
       shower_hadronic_tag=false,
       shower_hadronic_min_len=null,
       shower_hadronic_scan_len=null,
       shower_hadronic_bin=null,
       shower_hadronic_r_cyl=null,
       shower_hadronic_r_core=null,
       shower_hadronic_growth_max=null,
       shower_hadronic_growth_bragg=null,
       shower_hadronic_bragg_ratio=null,
       shower_hadronic_stem_ratio=null,
       kine_count_orphan_tracks=false,
       kine_orphan_track_min=null,
       straight_cont_cross_cluster=false,
       sccc_bridge_body=false,
       sccc_max_gap=null,
       sccc_kink_max=null,
       sccc_gap_aligned=null,
       sccc_kink_tight=null,
       // doc pr/43 round 2 -- C++ defaults false.
       single_muon_proton_chain_veto=false,
       single_muon_long_muon_claim=false,
       pid_flag_reconcile=false,
       // doc pr/45 -- C++ defaults false.
       other_seg_empty_2d_guard=false,
       // doc pr/46 -- C++ default false.
       long_muon_stub_bridge=false,
       // doc pr/48 -- back-to-back track fixes.  C++ defaults false / null =
       // C++ operating point; OFF here -- the SBND operating point (ALL THREE
       // ON, owner round doc pr/48 sec 9) lives in wct-pr-perevt.jsonnet.
       two_end_break=false,
       teb_min_len=null,
       teb_min_arm=null,
       teb_min_arm_pts=null,
       teb_stub_max=null,
       teb_accept_range=null,
       teb_rise_r1=null,
       teb_rise_r2=null,
       teb_abs_end_min=null,
       teb_dip_floor=null,
       teb_score_cap_r1=null,
       teb_score_cap_r2=null,
       teb_turn_angle=null,
       teb_turn_baseline=null,
       teb_turn_skirt=null,
       // doc pr/90 round 2: R2 argmax arm-fill guard + second-prong gate
       // cap.  null = C++ default 0 = legacy, key suppressed.
       teb_turn_min_arm_frac=null,
       teb_second_max=null,
       // doc pr/90 round 4: chain-topology admission (D1), route R3
       // turn/activity (D3), R2 bragg veto (D4).  false/null = C++ default
       // = legacy, key suppressed.
       teb_chain_topology=false,
       teb_r3_turn=null,
       teb_r3_hot=null,
       teb_bragg_veto_turn=null,
       kink_walk_dqdx_stop=false,
       kink_break_protect=false,
       kink_dqdx_hot_ratio=null,
       // doc pr/49 cross-cluster projection-ghost fit filter: null = C++
       // default -1 (off, key suppressed => byte-identical); >= 0 = on with
       // tolerance cells.
       fit_blob_coverage=null,
              // doc pr/50 -- suspend the pr/49 deweighting during
              // find_proto_vertex (the partition-forming stage; its recursive
              // kink walk is globally sensitive to fit perturbations --
              // 18255-172230 lost its true-kink main vertex to a 2.7 cm
              // neighbor).  All later fitting stages keep the deweighting.
              // C++ default false = pr/49 behavior; false omits the key =>
              // byte-identical.
              fit_blob_coverage_defer=false,
              // doc pr/50 -- main-vertex kink-consistency snap (172230-class
              // near-vertex robustness; C++ defaults in TaggerCheckNeutrino.h).
              // false/null omit the keys => byte-identical.
              vertex_kink_snap=false,
              vks_radius=null,
              vks_min_dis=null,
              vks_angle=null,
              vks_margin=null,
              vks_collinear=null,
              vks_skirt=null,
              vks_baseline=null,
              vks_min_arm=null,
              vks_fit_miss=null,
              vks_hot_ratio=null,
              // doc pr/85 -- carry the old vertex's arms through the snap
              // residual below this arc (cm).  C++ default 0 = off; null
              // omits the key => byte-identical.
              vks_carry_prong=null,
              // doc pr/104 -- main-vertex junction snap (C++ defaults in
              // TaggerCheckNeutrino.h).  false/null omit the keys => byte-identical.
              vertex_junction_snap=false,
              vjs_radius=null,
              vjs_min_arm=null,
              vjs_min_prongs=null,
              vjs_collinear=null,
              vjs_fit_margin=null,
              vjs_fit_rms=null,
              vjs_override_kink_snap=false,
              vjs_min_move=null,
              // doc pr/51 -- main-vertex graph audit + DL rerank
              // cross-cluster swap guard (506746).  false/null omit the
              // keys => byte-identical (see the clus_pr arg comments).
              // docs/73 sec 12 round 3; SBND PRODUCTION ON 2026-08-17.
              // See clus_pr.
              esva_ignore_empty_2d=true,
              main_vertex_graph_audit=false,
              mvga_radius=null,
              mvga_dup_tol=null,
              mvga_dup_frac=null,
              mvga_dup_angle=null,
              mvga_bridge_mip=null,
              mvga_reconnect=null,
              mvga_stub=null,
              mvga_stub_pts=null,
              mvga_reseat_angle=null,
              mvga_satellite=null,
              // doc pr/85 -- op3 interposed-stub absorb at the main-vertex
              // anchor (angle in deg, C++ default 150).  false/null omit
              // the keys => byte-identical.
              mvga_interposed=false,
              mvga_interposed_angle=null,
              // doc pr/86: interposed-splice candidate ceiling, cm (C++
              // default 0 = use mvga_stub).  null omits the key =>
              // byte-identical.
              mvga_interposed_len=null,
              // doc pr/86 P4: satellite-anchor op3 overlap threshold (C++
              // default 0 = use mvga_dup_frac).  null omits the key =>
              // byte-identical.
              mvga_sat_dup_frac=null,
              // doc pr/86 P1b: interposed splice at degree-1 main anchors
              // (C++ default false).  false omits the key => byte-identical.
              mvga_interposed_deg1=false,
              // doc pr/86 round 2: op3 post-carry straighten reach (cm) and
              // op3.5 junction-collapse radius (cm).  C++ defaults 0 (off).
              // Keys omitted when null => byte-identical.
              mvga_splice_straighten=null,
              mvga_approach_collapse=null,
              mvga_straighten_radius=null,
              // doc pr/83 r3: op1 scope/threshold decouple (cm / fraction;
              // radius -1 = unscoped), post-op3 dup pass, carry cap,
              // abandoned-cluster dup audit.  C++ defaults 0/0/false/0/false
              // (all legacy).  Keys omitted when null/false =>
              // byte-identical.
              mvga_op1_radius=null,
              mvga_op1_dup_frac=null,
              mvga_op1_post=false,
              mvga_carry_max=null,
              swap_orphan_dup_audit=false,
              // doc pr/83 r4 -- projective dup collapse; null omits => byte-identical.
              mvga_proj_dup_frac=null,
              mvga_proj_dqdx_ratio=null,
              mvga_proj_angle=null,
              mvga_ac_veto_radius=null,
              mvga_ac_chord_max=null,
              mvga_ac_no_cascade=false,
              mvga_passthru=null,
              mvga_passthru_tol=null,
              mvga_interposed_fallback=false,
              mvga_interposed_fallback_min_angle=null,
              mvga_dup_starved_asym=null,
              mvga_dup_starved_mip=null,
              mvga_dup_starved_span=null,
              dl_vtx_swap_guard=false,
              // doc pr/106 sec 10: exclusion-free charge cloud for the DL vertex net. C++ default false.
              dl_vtx_cloud_no_exclusion=false,
              // doc pr/89 Arm C (C2) -- rule-1 topology term in the DL
              // rerank composite (weight, frac center).  C++ defaults 0/0;
              // null omits the keys => byte-identical.
              dl_vtx_topo_weight=null,
              dl_vtx_topo_center=null,
              main_vertex_swap_apply=false,
              // doc pr/51 round 4 -- diagnostic-only rough-path probe.
              // C++ default false; false omits the key => byte-identical.
              rough_path_probe=false,
              // doc pr/51 round 5 -- steiner gap penalty (H1 short-cut fix).
              // C++ defaults: scale 0 (off), dead_alpha 0.25, min_edge
              // 0.5 cm, sample_step 0.3 cm, point_radius 0.2 cm.  null omits
              // the keys => byte-identical.
              steiner_gap_penalty=null,
              sgp_dead_alpha=null,
              sgp_min_edge=null,
              sgp_sample_step=null,
              sgp_point_radius=null,
              // doc pr/73: per-edge DEBUG sentinel for the steiner_graph_gap
              // scan (endpoints, midpoint, w, bad, both vertex charges,
              // deficit).  Log-only diagnostic.  C++ default false; false
              // omits the key => byte-identical compiled config.
              sgp_edge_probe=false, vertex_scoreboard=false, dl_vtx_harvest=false,
              // doc pr/51 round 6 -- weak-charge deficit term on the same
              // gap flavor.  C++ defaults: weak_scale 0 (off), weak_qref
              // 2000 (charge units).  null omits the keys => byte-identical.
              sgp_weak_scale=null,
              sgp_weak_qref=null,
              // doc pr/73 round 2 F3a -- do_rough_path route excursion cap,
              // cm.  C++ default -1 = off; null omits the key.
              sgp_max_sep=null,
              // doc pr/83 -- oriented break_segment splits (find_vertices, not
              // boost source/target).  C++ default false; key omitted when
              // off => byte-identical pre-fix config.
              break_seg_orient=false,
              // doc pr/51 round 7 -- robust vertex fit.  C++ defaults as in
              // clus_pr above.  false/null omit the keys => byte-identical.
              mvfit_robust=false,
              mvfit_main_only=null,
              mvfit_min_len=null,
              mvfit_rin_margin=null,
              mvfit_rout_frac=null,
              mvfit_rout_min=null,
              mvfit_rout_max=null,
              mvfit_angle=null,
              mvfit_min_pts=null,
              mvfit_min_aniso=null,
              mvfit_prior_range=null,
              // doc pr/54 -- keep well-supported isolated residual segments
              // in find_other_segments.  C++ defaults: keep false, floors
              // 25 points / 3 cm.  false/null omit the keys => byte-identical.
              other_seg_keep_isolated=false,
              other_seg_keep_isolated_min_points=null,
              other_seg_keep_isolated_min_length=null,
              // doc pr/102 P1 -- OR-disjuncts on the pr/54 keep: min_nnf
              // (terminal not-faked floor) and len_admit (cm).  C++ defaults
              // 0 = off.  null omits the keys => byte-identical.
              other_seg_keep_isolated_min_nnf=null,
              other_seg_keep_isolated_len_admit=null,
              // doc pr/102 P2 -- 3-D uncovered-charge radius (cm) for the
              // find_other_segments tagging/nnf seats.  C++ default 0 = off.
              // null omits the key => byte-identical.
              other_seg_uncover_3d=null,
              // doc pr/67 round 3 (S2) -- isochronous-snap size gate, cm.
              // C++ default 10.0 = legacy.  null omits the key.
              iso_snap_min_dir_mag=null,
              // doc pr/65 round 3 -- offer graph-unreachable main-cluster
              // segments (kept-isolated pr/54 residuals) to the shower
              // absorbers (reachability-relaxed guards).  C++ default false;
              // false omits the key => byte-identical.
              shower_absorb_unreachable_main=false,
              // doc pr/59 round 2 -- per-cluster orphaned-associate_points
              // rescue.  C++ default false; false omits the key =>
              // byte-identical.
              assoc_full_recluster=false,
              // doc pr/64 round 7 -- reassign same-cluster association
              // orphans that Stage C of clustering_points_segments would
              // otherwise drop.  C++ default false; false omits the key =>
              // byte-identical.  See the clus_pr arg comment.
              assoc_reassign_orphans=false,
              // doc pr/64 round 8 -- clear a merge survivor's
              // associate_points on structural-merge deletion.  C++ default
              // false; false omits the key => byte-identical.  See the
              // clus_pr arg comment.
              assoc_clear_on_merge=false,
              // doc pr/72 round 2 -- examine_structure_3 stub guard.  C++
              // default false/null; false/null omits the key =>
              // byte-identical.  See the clus_pr arg comment.
              es3_stub_guard=false,
              es3sg_stub_max=null,
              es3sg_len_ratio=null,
              es3sg_ang3_min=null,
              es3sg_ang_ratio=null,
              es3sg_require_terminal=null,
       pseudo_shower_track_paint=false,
       // Muon dQ/dx-vs-length envelope: DEFAULT = the docs/pr/10 SBND fit
       // (see the clus_pr arg comment; null restores the uBooNE refit).
       // Recombination-model selection + single-photon dE/dx routing:
       // DEFAULT ON (owner 2026-07-30, docs/pr/10 sec 7) -- see the clus_pr
       // arg comments; false/null restore the pre-pr/10 legacy behavior.
       muon_dqdx_curve=[0.8826, 1.0587, 18, 0.4745],
       use_power_recomb=true,
       sp_dedx_use_recomb_model=true, sp_mean_dedx_cut=2.23,
       // dl_vtx_cut (mm) is threaded for configurability only (docs/pr/2
       // sec 7.4); null keeps the C++ 25.0 (= 2.5 cm) default.
       dl_vtx_cut=null)::
        clus_pr(anodes, dump=dump,
                output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                pipeline_names=pipeline_names, tensor_outname=tensor_outname,
                trackfitting_config_file=trackfitting_config_file,
                particle_dataset=particle_dataset, extra_uses=extra_uses,
                dl_weights=dl_weights,
                dl_vtx_min_accept_score=dl_vtx_min_accept_score,
                dl_vtx_top_k=dl_vtx_top_k,
                dl_vtx_rerank=dl_vtx_rerank,
                beam_window=beam_window,
                tgm_neutrino_candidate=tgm_neutrino_candidate,
                tgm_chord_charge=tgm_chord_charge,
                tgm_chord_mode=tgm_chord_mode,
                tgm_component_extremes=tgm_component_extremes,
                tgm_component_rescue=tgm_component_rescue,
                tgm_rescue_chord=tgm_rescue_chord,
                tgm_main_pair=tgm_main_pair,
                tgm_main_pair_mode=tgm_main_pair_mode,
                tgm_fv_zmax_margin=tgm_fv_zmax_margin,
                tgm_fv_zmax_margin_interior=tgm_fv_zmax_margin_interior,
                tgm_fv_x_margin=tgm_fv_x_margin,
                tgm_fv_y_margin=tgm_fv_y_margin,
                save_stm_fit=save_stm_fit,
                pf_track_main_cluster_only=pf_track_main_cluster_only,
                pf_track_bridged_clusters=pf_track_bridged_clusters,
                pf_shower_vertex_barrier=pf_shower_vertex_barrier,
                pf_shower_parent_precedence=pf_shower_parent_precedence,
                pf_pi0_node_per_id=pf_pi0_node_per_id,
                pf_pdg_name_prototype_fallback=pf_pdg_name_prototype_fallback,
                pf_orphan_track_parentage=pf_orphan_track_parentage,
                pf_orphan_audit_only=pf_orphan_audit_only,
                pf_direct_when_touching=pf_direct_when_touching,
                pf_touch_max=pf_touch_max,
                pf_touch_cross_main=pf_touch_cross_main,
                pf_touch_cross_max=pf_touch_cross_max,
                pf_pseudo_gap_from_main=pf_pseudo_gap_from_main,
                pf_unique_node_ids=pf_unique_node_ids,
                pf_drop_stray_satellites=pf_drop_stray_satellites,
                pf_orphan_confident_track=pf_orphan_confident_track,
                pf_orphan_track_min_cm=pf_orphan_track_min_cm,
                pf_track_owns_loose_vertex=pf_track_owns_loose_vertex,
                unmerge_bundle_mode=unmerge_bundle_mode,
                restore_demoted_mains=restore_demoted_mains,
                require_provenance=require_provenance,
                evaluate_demoted_mains=evaluate_demoted_mains,
                tgm_exempt_demoted_main=tgm_exempt_demoted_main,
                skip_cosmic_companions=skip_cosmic_companions,
                cosmic_companion_min_length=cosmic_companion_min_length,
                nu_fallback_demoted_mains=nu_fallback_demoted_mains,
                sp_photon_flag=sp_photon_flag,
                mip_dqdx=mip_dqdx,
                stm_consistent_fv=stm_consistent_fv,
                stm_accept_guards=stm_accept_guards,
                stm_proton_muon_guard=stm_proton_muon_guard,
                stm_cathode_guard=stm_cathode_guard,
                stm_anode_dist_fix=stm_anode_dist_fix,
                stm_second_track_guard=stm_second_track_guard,
                stm_deficit_guard=stm_deficit_guard,
                stm_vertex_kink_guard=stm_vertex_kink_guard,
                stm_d66_cuts=stm_d66_cuts,
                stm_michel_res_cm=stm_michel_res_cm,
                stm_proton_tm_max=stm_proton_tm_max,
                stm_proton_b_ks2_max=stm_proton_b_ks2_max,
                stm_proton_c_peak_max=stm_proton_c_peak_max,
                beam_window_only=beam_window_only,
                nu_skip_cosmic=nu_skip_cosmic,
                nu_skip_cosmic_bundle=nu_skip_cosmic_bundle,
                nu_skip_cosmic_bundle_min_length=nu_skip_cosmic_bundle_min_length,
                dir_weak_use_score=dir_weak_use_score,
                mip_dqdx_median=mip_dqdx_median,
                proton_dir_vote=proton_dir_vote,
                endpoint_trim_retry=endpoint_trim_retry,
                fit_vertex_min_seg_length=fit_vertex_min_seg_length,
                cathode_x=cathode_x,
                cathode_kink_xcut=cathode_kink_xcut,
                cathode_wide_kink_angle=cathode_wide_kink_angle,
                cathode_wide_kink_skirt=cathode_wide_kink_skirt,
                cathode_wide_kink_baseline=cathode_wide_kink_baseline,
                shower_topo_demote_len=shower_topo_demote_len,
                fit_exclusion=fit_exclusion,
                graph_endpoint_strict=graph_endpoint_strict,
                graph_endpoint_tol=graph_endpoint_tol,
                oov_prototype_parity=oov_prototype_parity,
                first_seg_local_pca=first_seg_local_pca,
                other_seg_relaxed_accept=other_seg_relaxed_accept,
                shower_topo_proto_dir=shower_topo_proto_dir,
            cont_muon_dir3_30cm=cont_muon_dir3_30cm,
            track_comp_empty_abstain=track_comp_empty_abstain,
            shower_topo_reset=shower_topo_reset,
            reclass_preserve_4mom=reclass_preserve_4mom,
            dir_track_median_local=dir_track_median_local,
            examine_showers_vertex_by_index=examine_showers_vertex_by_index,
                vertex_dir_use_fit_point=vertex_dir_use_fit_point,
                shower_traj_recheck_parity=shower_traj_recheck_parity,
                main_vertex_require_descriptor=main_vertex_require_descriptor,
                main_vertex_candidate_flag=main_vertex_candidate_flag,
                steiner_terminal_wire_tol=steiner_terminal_wire_tol,
                steiner_terminal_adjacent_slice=steiner_terminal_adjacent_slice,
                steiner_edge_charge_forward_dead_mix=steiner_edge_charge_forward_dead_mix,
                iso_endpoint=iso_endpoint,
                iso_endpoint_min_length=iso_endpoint_min_length,
                iso_endpoint_max_xext=iso_endpoint_max_xext,
                iso_endpoint_xext_frac=iso_endpoint_xext_frac,
                iso_endpoint_xext_quantile=iso_endpoint_xext_quantile,
                iso_endpoint_tube_radius=iso_endpoint_tube_radius,
                iso_endpoint_min_aspect=iso_endpoint_min_aspect,
                v3_extension_guard=v3_extension_guard,
                traj_cover_probe=traj_cover_probe,
                pr_find_other_rounds=pr_find_other_rounds,
                v3_extension_min_gain=v3_extension_min_gain,
                protect_graph_name=protect_graph_name,
                protect_skip_convicted=protect_skip_convicted,
                protect_open_convicted_bundles=protect_open_convicted_bundles,
                protect_cathode_x=protect_cathode_x,
                protect_cathode_rejoin_xcut=protect_cathode_rejoin_xcut,
                protect_cathode_rejoin_dyz=protect_cathode_rejoin_dyz,
                protect_cathode_rejoin_dis=protect_cathode_rejoin_dis,
                protect_cathode_rejoin_perp=protect_cathode_rejoin_perp,
                protect_cathode_rejoin_angle=protect_cathode_rejoin_angle,
                protect_cathode_rejoin_dir_radius=protect_cathode_rejoin_dir_radius,
                protect_cathode_rejoin_dir_npts=protect_cathode_rejoin_dir_npts,
                cosmic_y_top_main=cosmic_y_top_main,
                cosmic_y_top_strict=cosmic_y_top_strict,
                cosmic_y_top_loose=cosmic_y_top_loose,
                cosmic_y_small_piece=cosmic_y_small_piece,
                vertex_z_prior_scale=vertex_z_prior_scale,
                ssm_target_dir=ssm_target_dir,
                ssm_absorber_dir=ssm_absorber_dir,
                kine_fudge_factor=kine_fudge_factor,
                kine_recom_factor=kine_recom_factor,
                kine_shower_fudge_factor=kine_shower_fudge_factor,
                kine_shower_recom_factor=kine_shower_recom_factor,
                kine_proton_recom_factor=kine_proton_recom_factor,
                kine_plane_weights=kine_plane_weights,
                kine_plane_asym_switch=kine_plane_asym_switch,
                kine_w_value=kine_w_value,
                kine_shower_pdg_live=kine_shower_pdg_live,
                neutrino_consistent_fv=neutrino_consistent_fv,
                cosmic_consistent_fv=cosmic_consistent_fv,
                nue_sp_consistent_fv=nue_sp_consistent_fv,
                sp_sce_correction=sp_sce_correction,
                tagger_ordered_segment_sets=tagger_ordered_segment_sets,
                stem_endpoint_wcpt_parity=stem_endpoint_wcpt_parity,
                broken_muon_cluster_id_count=broken_muon_cluster_id_count,
                neutrino_type_bitmask=neutrino_type_bitmask,
                nu_per_bundle=nu_per_bundle,
                nu_per_bundle_min_length=nu_per_bundle_min_length,
                nu_selected_as_main=nu_selected_as_main,
                nu_selected_as_main_snapshot_all=nu_selected_as_main_snapshot_all,
                daughter_count_proto_main_vertex=daughter_count_proto_main_vertex,
                daughter_count_proto_examine_showers=daughter_count_proto_examine_showers,
                shower_pdg_from_start_segment=shower_pdg_from_start_segment,
                shower_pdg_from_shower_type=shower_pdg_from_shower_type,
                shower_pdg_exact_muon_test=shower_pdg_exact_muon_test,
                pi0_id_shared_allocator=pi0_id_shared_allocator,
                shower_flag_pdg_electron=shower_flag_pdg_electron,
                shower_less_id_tiebreak=shower_less_id_tiebreak,
                shower_endpoint_exclude_start_vertex=shower_endpoint_exclude_start_vertex,
                shower_endpoint_skip_orphan_vtx=shower_endpoint_skip_orphan_vtx,
                shower_walk_visited_parity=shower_walk_visited_parity,
                track_pid_persist_dqdx=track_pid_persist_dqdx,
                shower_reclass_dqdx_guard=shower_reclass_dqdx_guard,
                shower_topo_dqdx_guard=shower_topo_dqdx_guard,
                reclass_never_computed_ke_floor=reclass_never_computed_ke_floor,
                track_pid_persist_4mom=track_pid_persist_4mom,
                shower_proton_daughter_pion=shower_proton_daughter_pion,
                shower_proton_daughter_pion_dissolve=shower_proton_daughter_pion_dissolve,
                muon_multi_proton_pion=muon_multi_proton_pion,
                track_pid_persist_dqdx_electron_guard=track_pid_persist_dqdx_electron_guard,
                shower_connect_main_vertex_straight_guard=shower_connect_main_vertex_straight_guard,
                shower_traj_straight_guard=shower_traj_straight_guard,
                shower_absorb_track_guard=shower_absorb_track_guard,
                shower_connect_protected_pion_guard=shower_connect_protected_pion_guard,
                michel_stem_muon_rescue=michel_stem_muon_rescue,
                shower_in_cascade_guard=shower_in_cascade_guard,
                shower_in_max_len=shower_in_max_len,
                shower_in_mip_hi=shower_in_mip_hi,
                shower_connect_from_vertices_straight_guard=shower_connect_from_vertices_straight_guard,
                shower_connect_start_seg_straight_guard=shower_connect_start_seg_straight_guard,
                examine_direction_dirsign_shower_in_guard=examine_direction_dirsign_shower_in_guard,
                daughter_shower_angle_reclass_straight_guard=daughter_shower_angle_reclass_straight_guard,
                shower_topo_reexam_straight_guard=shower_topo_reexam_straight_guard,
                sfv_kink_max=sfv_kink_max,
                shower_nv_bridge_track=shower_nv_bridge_track,
                shower_nv_bridge_max_gap=shower_nv_bridge_max_gap,
                shower_nv_main_pi_init=shower_nv_main_pi_init,
                kine_drop_stray_satellites=kine_drop_stray_satellites,
                kine_sat_min_energy=kine_sat_min_energy,
                kine_sat_prox_max=kine_sat_prox_max,
                kine_sat_angle_bad=kine_sat_angle_bad,
                kine_sat_angle_main=kine_sat_angle_main,
                kine_sat_far_dis=kine_sat_far_dis,
                kine_sat_axis_dis_cut=kine_sat_axis_dis_cut,
                kine_sat_cont_kink=kine_sat_cont_kink,
                kine_sat_track_max_nseg=kine_sat_track_max_nseg,
                kine_sat_em_far_dis=kine_sat_em_far_dis,
                michel_stem_michel_check=michel_stem_michel_check,
                michel_stem_max_far_len=michel_stem_max_far_len,
                shower_stem_backfill=shower_stem_backfill,
                stem_backfill_max_len=stem_backfill_max_len,
                stem_backfill_mip_lo=stem_backfill_mip_lo,
                stem_backfill_mip_hi=stem_backfill_mip_hi,
                stem_backfill_min_shower_len=stem_backfill_min_shower_len,
                shower_conn3_unreachable=shower_conn3_unreachable,
                conn3_unreachable_min_len=conn3_unreachable_min_len,
                conn3_stitch_max=conn3_stitch_max,
                shower_dedup_start_seg=shower_dedup_start_seg,
                shower_traj_michel_stem=shower_traj_michel_stem,
                michel_stem_traj_min_len=michel_stem_traj_min_len,
                michel_stem_traj_max_len=michel_stem_traj_max_len,
                michel_stem_traj_mip_lo=michel_stem_traj_mip_lo,
                michel_stem_traj_max_far_len=michel_stem_traj_max_far_len,
                michel_stem_traj_min_kink_deg=michel_stem_traj_min_kink_deg,
                shower_long_muon_keep_type=shower_long_muon_keep_type,
                shower_bragg_protect_start_segment=shower_bragg_protect_start_segment,
                shower_reclass_case_b_dqdx_guard=shower_reclass_case_b_dqdx_guard,
                shower_accept_pid_guard=shower_accept_pid_guard,
                shower_pid_guard_min_len=shower_pid_guard_min_len,
                shower_vote_track_pid_counts=shower_vote_track_pid_counts,
                shower_cone_absorb_guard=shower_cone_absorb_guard,
                shower_detach_track_stem=shower_detach_track_stem,
                shower_ghost_member_drop=shower_ghost_member_drop,
                shower_ghost_overlap_frac=shower_ghost_overlap_frac,
                shower_ghost_dqdx_ratio=shower_ghost_dqdx_ratio,
                shower_ghost_min_len=shower_ghost_min_len,
                kine_charge_dedup=kine_charge_dedup,
                kine_charge_rebuild=kine_charge_rebuild,
                kine_charge_track_ctx=kine_charge_track_ctx,
                kine_mass_rules=kine_mass_rules,
                kine_hadronic_dqdx=kine_hadronic_dqdx,
                kine_long_muon_mode=kine_long_muon_mode,
                kine_long_muon_ratio_lo=kine_long_muon_ratio_lo,
                kine_long_muon_ratio_hi=kine_long_muon_ratio_hi,
                kine_mainvtx_used_guard=kine_mainvtx_used_guard,
                shower_hadronic_tag=shower_hadronic_tag,
                shower_hadronic_min_len=shower_hadronic_min_len,
                shower_hadronic_scan_len=shower_hadronic_scan_len,
                shower_hadronic_bin=shower_hadronic_bin,
                shower_hadronic_r_cyl=shower_hadronic_r_cyl,
                shower_hadronic_r_core=shower_hadronic_r_core,
                shower_hadronic_growth_max=shower_hadronic_growth_max,
                shower_hadronic_growth_bragg=shower_hadronic_growth_bragg,
                shower_hadronic_bragg_ratio=shower_hadronic_bragg_ratio,
                shower_hadronic_stem_ratio=shower_hadronic_stem_ratio,
                kine_count_orphan_tracks=kine_count_orphan_tracks,
                kine_orphan_track_min=kine_orphan_track_min,
                straight_cont_cross_cluster=straight_cont_cross_cluster,
                sccc_bridge_body=sccc_bridge_body,
                sccc_max_gap=sccc_max_gap,
                sccc_kink_max=sccc_kink_max,
                sccc_gap_aligned=sccc_gap_aligned,
                sccc_kink_tight=sccc_kink_tight,
                single_muon_proton_chain_veto=single_muon_proton_chain_veto,
                single_muon_long_muon_claim=single_muon_long_muon_claim,
                pid_flag_reconcile=pid_flag_reconcile,
                other_seg_empty_2d_guard=other_seg_empty_2d_guard,
                long_muon_stub_bridge=long_muon_stub_bridge,
                two_end_break=two_end_break,
                teb_min_len=teb_min_len,
                teb_min_arm=teb_min_arm,
                teb_min_arm_pts=teb_min_arm_pts,
                teb_stub_max=teb_stub_max,
                teb_accept_range=teb_accept_range,
                teb_rise_r1=teb_rise_r1,
                teb_rise_r2=teb_rise_r2,
                teb_abs_end_min=teb_abs_end_min,
                teb_dip_floor=teb_dip_floor,
                teb_score_cap_r1=teb_score_cap_r1,
                teb_score_cap_r2=teb_score_cap_r2,
                teb_turn_angle=teb_turn_angle,
                teb_turn_baseline=teb_turn_baseline,
                teb_turn_skirt=teb_turn_skirt,
                teb_turn_min_arm_frac=teb_turn_min_arm_frac,
                teb_second_max=teb_second_max,
                teb_chain_topology=teb_chain_topology,
                teb_r3_turn=teb_r3_turn,
                teb_r3_hot=teb_r3_hot,
                teb_bragg_veto_turn=teb_bragg_veto_turn,
                kink_walk_dqdx_stop=kink_walk_dqdx_stop,
                kink_break_protect=kink_break_protect,
                kink_dqdx_hot_ratio=kink_dqdx_hot_ratio,
                fit_blob_coverage=fit_blob_coverage,
                fit_blob_coverage_defer=fit_blob_coverage_defer,
                vertex_kink_snap=vertex_kink_snap,
                vks_radius=vks_radius,
                vks_min_dis=vks_min_dis,
                vks_angle=vks_angle,
                vks_margin=vks_margin,
                vks_collinear=vks_collinear,
                vks_skirt=vks_skirt,
                vks_baseline=vks_baseline,
                vks_min_arm=vks_min_arm,
                vks_fit_miss=vks_fit_miss,
                vks_hot_ratio=vks_hot_ratio,
                vks_carry_prong=vks_carry_prong,
                vertex_junction_snap=vertex_junction_snap,
                vjs_radius=vjs_radius,
                vjs_min_arm=vjs_min_arm,
                vjs_min_prongs=vjs_min_prongs,
                vjs_collinear=vjs_collinear,
                vjs_fit_margin=vjs_fit_margin,
                vjs_fit_rms=vjs_fit_rms,
                vjs_override_kink_snap=vjs_override_kink_snap,
                vjs_min_move=vjs_min_move,
                esva_ignore_empty_2d=esva_ignore_empty_2d,
                main_vertex_graph_audit=main_vertex_graph_audit,
                mvga_radius=mvga_radius,
                mvga_dup_tol=mvga_dup_tol,
                mvga_dup_frac=mvga_dup_frac,
                mvga_dup_angle=mvga_dup_angle,
                mvga_bridge_mip=mvga_bridge_mip,
                mvga_reconnect=mvga_reconnect,
                mvga_stub=mvga_stub,
                mvga_stub_pts=mvga_stub_pts,
                mvga_reseat_angle=mvga_reseat_angle,
                mvga_satellite=mvga_satellite,
                mvga_interposed=mvga_interposed,
                mvga_interposed_angle=mvga_interposed_angle,
                mvga_interposed_len=mvga_interposed_len,
            mvga_sat_dup_frac=mvga_sat_dup_frac,
            mvga_interposed_deg1=mvga_interposed_deg1,
            mvga_splice_straighten=mvga_splice_straighten,
            mvga_approach_collapse=mvga_approach_collapse,
            mvga_straighten_radius=mvga_straighten_radius,
            mvga_op1_radius=mvga_op1_radius,
            mvga_op1_dup_frac=mvga_op1_dup_frac,
            mvga_op1_post=mvga_op1_post,
            mvga_carry_max=mvga_carry_max,
            swap_orphan_dup_audit=swap_orphan_dup_audit,
            mvga_proj_dup_frac=mvga_proj_dup_frac,
            mvga_proj_dqdx_ratio=mvga_proj_dqdx_ratio,
            mvga_proj_angle=mvga_proj_angle,
            mvga_ac_veto_radius=mvga_ac_veto_radius,
            mvga_ac_chord_max=mvga_ac_chord_max,
            mvga_ac_no_cascade=mvga_ac_no_cascade,
            mvga_passthru=mvga_passthru,
            mvga_passthru_tol=mvga_passthru_tol,
            mvga_interposed_fallback=mvga_interposed_fallback,
            mvga_interposed_fallback_min_angle=mvga_interposed_fallback_min_angle,
            mvga_dup_starved_asym=mvga_dup_starved_asym,
            mvga_dup_starved_mip=mvga_dup_starved_mip,
            mvga_dup_starved_span=mvga_dup_starved_span,
                dl_vtx_swap_guard=dl_vtx_swap_guard,
                dl_vtx_cloud_no_exclusion=dl_vtx_cloud_no_exclusion,
                dl_vtx_topo_weight=dl_vtx_topo_weight,
                dl_vtx_topo_center=dl_vtx_topo_center,
                main_vertex_swap_apply=main_vertex_swap_apply,
                rough_path_probe=rough_path_probe,
                steiner_gap_penalty=steiner_gap_penalty,
                sgp_dead_alpha=sgp_dead_alpha,
                sgp_min_edge=sgp_min_edge,
                sgp_sample_step=sgp_sample_step,
                sgp_point_radius=sgp_point_radius,
                sgp_edge_probe=sgp_edge_probe, vertex_scoreboard=vertex_scoreboard, dl_vtx_harvest=dl_vtx_harvest,
                sgp_weak_scale=sgp_weak_scale,
                sgp_weak_qref=sgp_weak_qref,
                sgp_max_sep=sgp_max_sep,
                break_seg_orient=break_seg_orient,
                mvfit_robust=mvfit_robust,
                mvfit_main_only=mvfit_main_only,
                mvfit_min_len=mvfit_min_len,
                mvfit_rin_margin=mvfit_rin_margin,
                mvfit_rout_frac=mvfit_rout_frac,
                mvfit_rout_min=mvfit_rout_min,
                mvfit_rout_max=mvfit_rout_max,
                mvfit_angle=mvfit_angle,
                mvfit_min_pts=mvfit_min_pts,
                mvfit_min_aniso=mvfit_min_aniso,
                mvfit_prior_range=mvfit_prior_range,
                other_seg_keep_isolated=other_seg_keep_isolated,
                other_seg_keep_isolated_min_points=other_seg_keep_isolated_min_points,
                other_seg_keep_isolated_min_length=other_seg_keep_isolated_min_length,
                other_seg_keep_isolated_min_nnf=other_seg_keep_isolated_min_nnf,
                other_seg_keep_isolated_len_admit=other_seg_keep_isolated_len_admit,
                other_seg_uncover_3d=other_seg_uncover_3d,
                iso_snap_min_dir_mag=iso_snap_min_dir_mag,
                shower_absorb_unreachable_main=shower_absorb_unreachable_main,
                assoc_full_recluster=assoc_full_recluster,
                assoc_reassign_orphans=assoc_reassign_orphans,
                assoc_clear_on_merge=assoc_clear_on_merge,
                es3_stub_guard=es3_stub_guard,
                es3sg_stub_max=es3sg_stub_max,
                es3sg_len_ratio=es3sg_len_ratio,
                es3sg_ang3_min=es3sg_ang3_min,
                es3sg_ang_ratio=es3sg_ang_ratio,
                es3sg_require_terminal=es3sg_require_terminal,
                pseudo_shower_track_paint=pseudo_shower_track_paint,
                muon_dqdx_curve=muon_dqdx_curve,
                use_power_recomb=use_power_recomb,
                sp_dedx_use_recomb_model=sp_dedx_use_recomb_model,
                sp_mean_dedx_cut=sp_mean_dedx_cut,
                dl_vtx_cut=dl_vtx_cut),
    detector_volumes(anodes, face=0):: detector_volumes(anodes=anodes, face=face, pos_offset_on=pos_offset_on),
    // Primitives the entry configuration needs to build the wclsTensorSetLabeler
    // node itself (it is no longer wired inside clus_all_apa).  All are the exact
    // objects the labeler used to receive, with the reality-correct pos_offset_on.
    pc_transforms(dv):: pctransforms(dv),
    // The corrected coordinate-array names the all-APA MABC uses for its Bee
    // (clustering_global): data ['x_t0cor','y_cor'/'y','z_cor'/'z'], sim
    // ['x_sce','y_sce','z_sce'].  The entry hands these to the labeler so the
    // tagger_stm/tgm/fc sets overlay clustering_global.
    bee_coords:: common_corr_coords(pos_offset_on, use_sce),
    sce_field_fwd:: sce_field_fwd,
    drift_speed:: drift_speed,
    time_offset:: time_offset,
    // FV box (spans both TPCs) for the labeler's particle-flow fiducial cut.
    fiducial_box()::
        local fvm = dvm(pos_offset_on).overall;
        {
            type: 'BoxFiducial',
            name: 'all-overall-fv',
            data: { bounds: {
                tail: { x: fvm.FV_xmin + fvm.FV_xmin_margin,
                        y: fvm.FV_ymin + fvm.FV_ymin_margin,
                        z: fvm.FV_zmin + fvm.FV_zmin_margin },
                head: { x: fvm.FV_xmax - fvm.FV_xmax_margin,
                        y: fvm.FV_ymax - fvm.FV_ymax_margin,
                        z: fvm.FV_zmax - fvm.FV_zmax_margin },
            } },
        },
}
