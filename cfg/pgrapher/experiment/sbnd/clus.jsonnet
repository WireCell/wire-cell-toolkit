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
local clus_per_face(anode, face, dump, output_dir, runNo, subRunNo, eventNo, bee_sink=null, rse_from_ident=false, pos_offset_on=true, trace_bee=false, save_assoc_id=false, sep_vertex_veto=true, nu_iso_band_guard=true, iso_cathode_guard=false) = {
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
        cm.neutrino(protect_iso_band=nu_iso_band_guard,
                    protect_iso_band_xext=(if nu_iso_band_guard then 20 * wc.cm else null)),
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
local clus_all_apa(anodes, dump, output_dir, runNo, subRunNo, eventNo, bee_sink=null, premerged=false, rse_from_ident=false, pos_offset_on=true, tensor_outname='', save_real_cluster_id=false, trace_bee=false, save_assoc_cluster_id=false, real_cluster_id_global=null, cathode_rescue_on=true, cathode_rescue_unmatched=true, adopt_nu_fragments=false) = {
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
        adopt_dis=(if adopt_nu_fragments then 13*wc.cm else null))] else [])
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
            ] + (if trace_bee then trace_sets(cm_pipeline, common_corr_coords(pos_offset_on)) else []),
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
        unmerge_bundle: cm.unmerge_bundle(mode=unmerge_bundle_mode),
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
        // SBND has no beam_flash flag (QLMatching sets main/associated_cluster
        // instead) -- process every scope-passing cluster, narrowed to the
        // beam-coincident bundle when beam_window_only is on (the default; see
        // the clus_pr argument).  Keys omitted when off => byte-identical.
        steiner: cm.steiner(retiler=improve2, perf=true, require_beam_flash=false,
                            beam_window_only=beam_gate,
                            beam_window_low=beam_window[0],
                            beam_window_high=beam_window[1]),
        fiducialutils: cm.fiducialutils(),
        tagger_check_stm: cm.tagger_check_stm(
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
            main_component_mode=tgm_main_pair_mode),
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
            // compiled JSON is unchanged by this pinning).
            dl_vtx_rerank=true,
            dl_vtx_top_k=5,
            dl_vtx_min_accept_score=4.0,
            dl_vtx_score_scale=1000.0,
            beam_window_low=beam_window[0],
            beam_window_high=beam_window[1],
            nu_skip_cosmic=nu_skip_cosmic,
            nu_skip_cosmic_bundle=nu_skip_cosmic_bundle,
            nu_skip_cosmic_bundle_min_length=nu_skip_cosmic_bundle_min_length,
            dir_weak_use_score=dir_weak_use_score,
            mip_dqdx=mip_dqdx,
            mip_dqdx_median=mip_dqdx_median,
            proton_dir_vote=proton_dir_vote,
            endpoint_trim_retry=endpoint_trim_retry,
            fit_vertex_min_seg_length=fit_vertex_min_seg_length,
            cathode_x=cathode_x,
            cathode_kink_xcut=cathode_kink_xcut,
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
        tagger_output: cm.tagger_output(output_filename=tracking_pr_root),
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
                    [if std.member(pipeline_names, 'tagger_check_neutrino')
                     then 'np_ke_min']: 10 * wc.MeV,
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
    // trace_bee (default false): per-step Bee layers for merge attribution; see
    // trace_sets above.  Diagnostic only, off => byte-identical compiled config.
    per_apa(anode, dump=true, bee_sink=null, trace_bee=false, save_assoc_id=false, sep_vertex_veto=true, nu_iso_band_guard=true, iso_cathode_guard=false)::
        clus_per_face(anode, face=0, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                      trace_bee=trace_bee, save_assoc_id=save_assoc_id, sep_vertex_veto=sep_vertex_veto,
                      nu_iso_band_guard=nu_iso_band_guard, iso_cathode_guard=iso_cathode_guard),
    // Production (LArSoft) entry point used by wcls-img-clus.jsonnet.
    per_volume(anode, face=0, dump=true, bee_sink=null)::
        clus_per_face(anode, face=face, dump=dump,
                      output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                      bee_sink=bee_sink, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on),
    all_apa(anodes, dump=true, bee_sink=null, premerged=false, tensor_outname='', save_real_cluster_id=false, save_assoc_cluster_id=false,
            trace_bee=false, real_cluster_id_global=null, cathode_rescue_on=true, cathode_rescue_unmatched=true, adopt_nu_fragments=false)::
        clus_all_apa(anodes, dump=dump,
                     output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                     bee_sink=bee_sink, premerged=premerged, rse_from_ident=rse_from_ident, pos_offset_on=pos_offset_on,
                     tensor_outname=tensor_outname, save_real_cluster_id=save_real_cluster_id,
                     save_assoc_cluster_id=save_assoc_cluster_id,
                     real_cluster_id_global=real_cluster_id_global,
                     trace_bee=trace_bee, cathode_rescue_on=cathode_rescue_on,
                     cathode_rescue_unmatched=cathode_rescue_unmatched,
                     adopt_nu_fragments=adopt_nu_fragments),
    // PR job: input is the reloaded post-QL tarball (see clus_pr above).
    // The TGM/FC and beam-window defaults here mirror clus_pr's -- i.e. the SBND
    // production operating point (see the comment block on clus_pr's arg list for
    // the pre-adoption values to pass for an A/B).
    pr(anodes, dump=true, pipeline_names=[], tensor_outname='',
       trackfitting_config_file='pgrapher/experiment/sbnd/sbnd_track_fitting.json',
       particle_dataset=null, extra_uses=[],
       // DL (SCN) vertex ON by default -- see the clus_pr arg comment.
       dl_weights='uboone/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth',
       beam_window=[0.2 * wc.us, 2.2 * wc.us],
       tgm_neutrino_candidate=true,
       tgm_chord_charge=true, tgm_chord_mode='path',
       tgm_component_extremes=true, tgm_component_rescue=true,
       tgm_rescue_chord=true, tgm_main_pair=true, tgm_main_pair_mode='real',
       tgm_fv_zmax_margin=5, tgm_fv_zmax_margin_interior=3,
       tgm_fv_x_margin=2.5, tgm_fv_y_margin=3,
       save_stm_fit=false, unmerge_bundle_mode='real',
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
       // C++ defaults are 0/0 = OFF; null here omits both keys, so the
       // compiled JSON is byte-identical to the pre-B0 config.
       cathode_x=null, cathode_kink_xcut=null,
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
                dl_weights=dl_weights, beam_window=beam_window,
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
                unmerge_bundle_mode=unmerge_bundle_mode,
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
                muon_dqdx_curve=muon_dqdx_curve,
                use_power_recomb=use_power_recomb,
                sp_dedx_use_recomb_model=sp_dedx_use_recomb_model,
                sp_mean_dedx_cut=sp_mean_dedx_cut,
                dl_vtx_cut=dl_vtx_cut),
    detector_volumes(anodes, face=0):: detector_volumes(anodes=anodes, face=face, pos_offset_on=pos_offset_on),
}
