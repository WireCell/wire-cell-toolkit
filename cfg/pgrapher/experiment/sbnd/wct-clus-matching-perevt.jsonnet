// CANONICAL in-tree SBND clustering + Q/L matching job — single source of truth.
// Promoted from sbnd_xin 2026-07-27 (sbnd_xin/docs/64_cfg-sync.md);
// sbnd_xin/wct-clus-matching-perevt.jsonnet is now a one-line re-export of this
// file.  Runnable directly:
//   wire-cell -c pgrapher/experiment/sbnd/wct-clus-matching-perevt.jsonnet --tla-...
// (WIRECELL_PATH must contain toolkit/cfg + wire-cell-data).
//
// Per-event SBND charge-light (Q/L) matching — standalone, self-contained.
//
// Reads the toolkit's OWN per-event imaging output (run_img_evt.sh →
// work/evt<ID>/icluster-apa<N>-{active,masked}.npz) plus that event's opflash,
// runs per-APA clustering + charge-light matching, and writes one
// work/evt<ID>/mabc-all-apa.zip (img + clustering + 2-view dead + op/Q-L layers).
//
// Unlike wct-clus-matching-standalone.jsonnet (all-10, yuhw's larsoft active
// clusters, masked imaged in-graph), this reads active AND masked from the
// per-event npz — no in-graph imaging — so it is self-contained per event and
// parallelizable.  Charge clusters are the toolkit's own (run_img_evt.sh), not
// yuhw's larsoft dumps.  Driven by run_ql_evt.sh.
//
// Usage (called from run_ql_evt.sh):
//   wire-cell \
//     --tla-str  input=work/evt2 \
//     --tla-code anode_indices='[0,1]' \
//     --tla-str  output_dir=work/evt2 \
//     --tla-code run=0 --tla-code subrun=0 --tla-code event=2 \
//     --tla-str  reality=sim \
//     --tla-str  semimodel_file=semi-analytical-sbnd.json \
//     --tla-code DL=4.0 --tla-code DT=8.8 \
//     --tla-code lifetime=35 --tla-code driftSpeed=1.563 \
//     -c wct-clus-matching-perevt.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

function(
    input          = '.',
    anode_indices  = [0, 1],
    output_dir     = '.',
    run            = 0,
    subrun         = 0,
    event          = 0,
    reality        = 'sim',
    semimodel_file = 'semi-analytical-sbnd.json',
    // DL/DT: sbndcode's production diffusion (wcsimsp_sbnd.fcl), restored
    // 2026-07-27, reverting the 6.5781/13.1349 retune of 2026-07-25 (docs/66).
    // Inert here in the same sense as lifetime below -- compiling this job at
    // either pair gives byte-identical JSON, zero 'DL' keys.
    DL             = 4.0,
    DT             = 8.8,
    // lifetime: electron lifetime in ms.  35 = the SBND simparams.jsonnet value,
    // so this chain and the simulation now state the same number (owner,
    // 2026-07-27).  It was 6.0, an inherited placeholder from the first standalone
    // Q/L test (wcp-porting-img 655bd6a: "DL=6.2 DT=9.8 lifetime=6
    // driftSpeed=1.565") whose DL/DT half was corrected to the SBND physical
    // values in 9f498089 while lifetime was not.
    // The change is provably inert: this knob only feeds the sim Drifter's
    // attenuation and NO reco component reads it -- compiling this job with
    // lifetime=6 and lifetime=35 gives byte-identical JSON (zero 'lifetime' keys).
    // NO electron-lifetime / charge-attenuation correction is applied anywhere in
    // imaging, clustering, Q/L matching or PR; adding one would be a real change
    // needing a measured data lifetime and a dQ/dx revisit, not this knob.
    // See sbnd_xin/docs/64_cfg-sync.md sec 4.
    lifetime       = 35.0,
    driftSpeed     = 1.563,
    // Joint multi-APA matching toggle. false (default) = historical per-APA path
    // (one QLMatching per APA -> PointTreeMerging -> all-APA MABC). true = one
    // joint QLMatching node matches both APAs and merges, feeding MABC directly.
    joint          = false,
    // Hand-scan calibration dump path. '' (default) = off, production-identical.
    // When set, QLMatching writes one per-event JSON (both TPCs) for the Q/L
    // hand-scan viewer (sbnd_xin/ql_scan). run_ql_evt.sh -calib points it at
    // work/evt<ID>/calib-evt<ID>.json.
    calib_dump     = '',
    // Cathode-crossing TPC0/TPC1 offset diagnostic. '' (default) = off,
    // production-identical. When set (any non-empty string), QLMatching logs the
    // three-vector (dir0, dir1, conn) decomposition per cross-TPC cathode-crossing
    // pair (grep "QLCATHODE" in the run log). run_ql_evt.sh -cathode-diag enables it.
    cathode_diag   = '',
    // Per-PMT predicted-light non-linearity correction. true (default, SBND going
    // forward) maps each PMT's accumulated predicted PE into the saturated (observed)
    // space (qlmatching.jsonnet -> pmt_nonlinearity_params.jsonnet). false = identity.
    // run_ql_evt.sh threads it via PMT_NL / --tla-code pmt_nl.
    pmt_nl         = true,
    // Per-event dynamic dead-PMT auto-mask. false (default) = off, production-identical.
    // true masks, per event, a PMT that never fires while its live neighbours do (a
    // run-dead channel absent from the static ch_mask). run_ql_evt.sh -auto-mask enables it.
    auto_mask      = false,
    // Beam-window flash preference OVERRIDE overlay. Since the validation-round-2
    // adoption the preference is ON in the production config
    // (cfg/.../sbnd/qlmatching.jsonnet: weight 0.5, rescue 0.2, gate ks 0.3 /
    // pred 2%; docs/22_ql-beam-flash-preference.md), so false (default) here
    // simply inherits production. true re-asserts the same keys via the shim
    // overlay -- useful only with beam_pref_weight/beam_pref_rescue below to
    // scan a different operating point (run_ql_evt.sh -beam-pref +
    // BEAMPREF_WEIGHT/BEAMPREF_RESCUE).
    beam_pref      = false,
    // Beam-preference operating point (inert while beam_pref=false; keys are
    // suppressed then, so the compiled config stays byte-identical). weight =
    // LASSO L1 multiplier for beam-window bundles (validated 0.5; 0.2
    // over-collects), rescue = empty-flash rescue steal guard (a non-beam flash
    // must beat a beam-window match by 1/rescue to re-steal). run_ql_evt.sh
    // threads BEAMPREF_WEIGHT / BEAMPREF_RESCUE env into these for scans.
    beam_pref_weight = 0.5,
    beam_pref_rescue = 0.2,
    // Stamp flag_main_cluster on every matched bundle main (QLMatching
    // flag_matched_mains).  DEFAULT TRUE = SBND production (run_ql_evt.sh's
    // MAINFLAG is 1 unless -no-main-flag).  false = legacy: only the mains that
    // decompose_cluster_groups SPLIT carry the flag, so a compact single-component
    // match is skipped by TaggerCheckTGM/STM/FC and shows up as "no-bundle" in the
    // nusel table (SBND evt286021, 1.158 us beam flash).  C++ default false; the
    // key is omitted when off => compiled config byte-identical to pre-fix.
    main_flag      = true,
    // LM (light-mismatch) tagger (QLMatching lm_tagger).  DEFAULT TRUE = SBND
    // production (run_full1k_nusel.sh passes -lm).  Every FINAL matched bundle is
    // judged by per-drift-side KS shape + pred/meas normalization; the verdict
    // (0 pass / 1 low-energy / 2 light mismatch) is stamped as cluster scalar
    // "lm_flag" (read by nusel_extract.py's lm column) and dumped per bundle
    // into the calib JSON (lm, lm_ks, lm_pred, lm_meas, lm_length_cm).  Pass
    // false for the pre-LM baseline (C++ default false; key omitted when off =>
    // compiled config byte-identical to pre-LM).
    lm             = true,
    // Persistent post-QL intermediate output. '' (default) = off, production-identical
    // (the terminal TensorFileSink stays a dump_mode no-op). When set to a path, the
    // all-APA MABC output point-cloud tree (live+dead, cluster_t0/flash annotations,
    // opflash PC) is written there as a TensorDM tar.gz -- the input of the
    // downstream pattern-recognition job (sbnd/docs/sbnd-pattern-recognition.md).
    // run_ql_evt.sh -save-pctree points it at work/ql_evt<ID>/pctree-evt<ID>.tar.gz.
    save_tensors   = '',
    // Persist the flash-merge per-blob provenance (real_cluster_id /
    // real_cluster_main "perblob" arrays) through the save_tensors tarball:
    // the tensor serializer drops heterogeneous PC keys, so without this the
    // arrays exist only in-memory and the PR job cannot tell which points
    // were the bundle's main cluster (doc 38).  Only meaningful WITH
    // save_tensors.  false (default) = byte-identical legacy tarball (C++
    // default false; key omitted when off).  Runner flag: -save-rcid.
    save_rcid      = false,
    // DIAGNOSTIC (default false): dump one Bee "clustering" layer per clustering
    // step (per-APA AND all-APA), named tr<NN>_<Type>, so a merge can be
    // attributed to the pass that made it instead of guessed.  See
    // cfg/.../sbnd/clus.jsonnet trace_sets() and docs/51.  Off => the
    // bee_points_sets lists are unchanged => compiled config byte-identical.
    // run_ql_evt.sh -trace-bee / SBND_TRACE_BEE=1.
    trace_bee      = false,
    // save_assoc (default false): turn on doc 52's isolated-grouping provenance.
    // clustering_isolated writes the per-blob pair assoc_cluster_id /
    // assoc_cluster_main recording which pre-merge cluster each blob came from
    // and which member was the MAIN; merge_clusters carries that pair across
    // every later merge (including cathode_connect and the flash merge) and
    // MABC homogenizes it into the save_tensors tarball.  The PR job then runs
    // a second ClusteringUnmergeBundle (pipeline visitor 'unmerge_assoc') to
    // undo the grouping, so TaggerCheckSTM / the PR chain fit the main alone --
    // the prototype's main_cluster + additional_clusters layout.
    // Only meaningful WITH save_tensors.  Off => both keys omitted => compiled
    // config byte-identical.  Runner flag: -save-assoc / SBND_SAVE_ASSOC=1.
    save_assoc     = false,
    // rcid_global (default null = inherit the C++ default, which is TRUE since
    // doc 53): re-stamp real_cluster_id at save time into ONE globally unique
    // ident epoch.  Without it the array mixes the numbering examine_bundles
    // recorded with the numbering enumerate_idents has installed since -- 31% of
    // values name two clusters on the d52ron 30-event set -- which is harmless
    // for the within-cluster consumers but wrong for anything that joins on the
    // value (including the Bee per-blob label).  The re-stamp runs right after
    // the clustering pipeline, before the Bee fills, so the saved pctree and the
    // mabc Bee zips carry the SAME ids; per-step trace_bee layers are
    // mid-pipeline snapshots and keep their own.
    // Group membership is identical either way, so
    // the un-merge and TGM are verdict-neutral (measured: 179 clusters, same
    // partition; nusel verdict tables row-for-row identical).
    // Only meaningful WITH save_tensors + save_rcid.  Pass false ONLY to
    // reproduce the two-epoch values for A/B archaeology.
    // Runner flag: -no-rcid-global / SBND_QL_RCID_GLOBAL=0.
    rcid_global    = null,
    // realign (default null = inherit the C++ default, which is TRUE since
    // doc 52 §12.8): QLMatching realign_perblob.  Pass false ONLY to
    // reproduce the pre-fix misaligned behavior for A/B archaeology (the
    // recomposed "isolated" rows then stay stale, changing the all-APA
    // examine_bundles main-overlap vote inputs back to the old ones).
    // Runner flag: -no-realign / SBND_REALIGN=0.
    realign        = null,
    // cathode_rescue (SBND default TRUE since the doc pr/14 §7.4 validation,
    // owner decision 2026-08-01): the cathode BUNDLE rescue between
    // cathode_connect and examine_bundles -- the isolated patch for the
    // flash-reco absorbing-window defect that leaves the two halves of a
    // cathode crosser in different flash bundles (sbnd_xin/docs/pr/14).
    // false restores the pre-pr/14 pipeline (compiled config byte-identical
    // to before the knob existed; runner: SBND_CATHODE_RESCUE=0).  Retire
    // when the light reconstruction is fixed upstream.
    cathode_rescue = true,
    // cathode_rescue_unmatched (SBND default TRUE since doc pr/17, validated
    // 2026-08-01): second rescue pass adopting a NON-MATCHED cluster into the
    // beam bundle when it geometrically continues a beam-window cluster
    // across the cathode (run 18255 evt 56463 with the vertex veto ON: the
    // rejoined nu is flashless and invisible downstream).  Needs
    // cathode_rescue.  Fires 1/1000 mcp1k, 0/48 nueCC48 (doc pr/17 sec 7).
    // false omits the key => compiled config byte-identical to pre-pr/17
    // (runner: SBND_RESCUE_UNMATCHED=0 to escape).
    cathode_rescue_unmatched = true,
    // sep_vertex_veto (SBND default TRUE since doc pr/15, owner decision
    // 2026-08-01): per-APA separate() un-splits a neutrino-vertex "V" whose two
    // dominant pieces both END at their mutual closest approach (run 18255 evt
    // 56463: the nu was cut in two at its vertex by the top-cosmic angle
    // ladder).  false omits the key => compiled config byte-identical to
    // before the knob existed (runner: SBND_SEP_VVETO=0).
    sep_vertex_veto = true,
    // nu_iso_band_guard (SBND default TRUE since doc pr/18, owner decision
    // 2026-08-01): the per-APA neutrino stage may not merge an isochronous
    // band with a non-band cluster spanning > 20 cm of drift, even on touch
    // (run 18255 evt 10550: separate correctly split the nu candidate off the
    // cosmic band, neutrino re-merged them at 0.31 cm).  false omits the keys
    // => byte-identical pre-knob config (runner: SBND_NU_ISO_GUARD=0).
    nu_iso_band_guard = true,
)
    // Build params inside the function so all physics values are TLAs.  These
    // are the documented Q/L drift/diffusion values (matching run_clust_QL_evt.sh),
    // not wct-clustering's 4.0/8.8/35.
    local base = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
    local params = base {
        lar: super.lar {
            DL:          DL * wc.cm2 / wc.s,
            DT:          DT * wc.cm2 / wc.s,
            lifetime:    lifetime * wc.ms,
            drift_speed: driftSpeed * wc.mm / wc.us,
        },
    };

    local tools_all = tools_maker(params);
    local anodes  = [tools_all.anodes[i] for i in anode_indices];
    local nanodes = std.length(anodes);

    // --- Charge sources: the toolkit's own per-event imaging output ---
    local cluster_source(fname) = g.pnode({
        type: 'ClusterFileSource',
        name: fname,
        data: { inname: fname, anodes: [wc.tn(a) for a in anodes] },
    }, nin=0, nout=1, uses=anodes);

    local active_clusters = [cluster_source('%s/icluster-apa%d-active.npz' % [input, a.data.ident]) for a in anodes];
    local masked_clusters = [cluster_source('%s/icluster-apa%d-masked.npz' % [input, a.data.ident]) for a in anodes];

    local clus_mod = import 'clus.jsonnet';
    local clus_maker = clus_mod(
        output_dir=output_dir,
        runNo=run,
        subRunNo=subrun,
        eventNo=event,
        reality=reality);
    local clus_pipes = [clus_maker.per_apa(anodes[n], dump=false, trace_bee=trace_bee, save_assoc_id=save_assoc, sep_vertex_veto=sep_vertex_veto, nu_iso_band_guard=nu_iso_band_guard)
                        for n in std.range(0, nanodes - 1)];

    // --- Q/L matching nodes ---
    // opflash TensorFileSource is built inline so it reads from the per-event
    // `input` dir (qlm.opflash_source hardcodes a bare cwd-relative filename).
    local qlm = (import 'qlmatching.jsonnet')(params);
    // SBND CPA structure-exclusion fiducial (cushion is the live tuning knob).
    local cathode_fv = (import 'cathode_fiducial.jsonnet')(cx=0.5*wc.cm, cy=0.5*wc.cm, cz=0.5*wc.cm);
    local opflash_sources = [g.pnode({
        type: 'TensorFileSource',
        name: 'opflash_src_apa%d' % anodes[n].data.ident,
        data: {
            inname: '%s/opflash_apa%d.tar.gz' % [input, anodes[n].data.ident],
            prefix: 'opflash_',
        },
    }, nin=0, nout=1) for n in std.range(0, nanodes - 1)];
    local flash_attach   = [qlm.flash_attach(n) for n in std.range(0, nanodes - 1)];
    local matching_pipes = [qlm.matching(anodes[n], clus_maker.detector_volumes([anodes[n]]),
                                         n, reality, semimodel_file,
                                         cathode_fiducial=cathode_fv.tn,
                                         calib_dump=calib_dump, cathode_diag=cathode_diag,
                                         pmt_nl=pmt_nl,
                                         // auto_mask / beam_pref are already ON at the SBND
                                         // production operating point inside qlmatching.jsonnet,
                                         // so these TLAs only RE-ASSERT them (that was also the
                                         // old sbnd_xin wrapper's behavior): false => null =>
                                         // key omitted => inherit production, NOT "off".  To
                                         // genuinely disable either one, pass auto_mask=false /
                                         // beam_pref=false to qlmatching.jsonnet directly.
                                         auto_mask=(if auto_mask then true else null),
                                         beam_pref=(if beam_pref then true else null),
                                         beam_pref_weight=(if beam_pref then beam_pref_weight else null),
                                         beam_pref_rescue=(if beam_pref then beam_pref_rescue else null),
                                         main_flag=main_flag, lm=lm, realign_perblob=realign)
                            for n in std.range(0, nanodes - 1)];

    // --- Graph: per-APA matching (default) or joint multi-APA matching ---
    local graph =
        if joint then
            // active ─┐
            // masked ─┼─► clus.per_apa ─► FlashTensorToOpticalPCs ─┐
            // opflash ┘                                            ├─► QLMatching(joint) ─► MABC
            //   (other APA's flash_attach) ───────────────────────┘     (merges per-APA trees)
            local jointql = qlm.matching_joint(anodes, clus_maker.detector_volumes(anodes),
                                               reality, semimodel_file,
                                               cathode_fiducial=cathode_fv.tn,
                                               calib_dump=calib_dump, cathode_diag=cathode_diag,
                                               pmt_nl=pmt_nl,
                                               // See the per-APA call above: false => inherit
                                               // production, not "off".
                                               auto_mask=(if auto_mask then true else null),
                                               beam_pref=(if beam_pref then true else null),
                                               beam_pref_weight=(if beam_pref then beam_pref_weight else null),
                                               beam_pref_rescue=(if beam_pref then beam_pref_rescue else null),
                                               main_flag=main_flag, lm=lm, realign_perblob=realign);
            // MABC takes the single pre-merged tree directly (no PointTreeMerging).
            local clus_all = clus_maker.all_apa(anodes, dump=true, premerged=true, tensor_outname=save_tensors, save_real_cluster_id=save_rcid, save_assoc_cluster_id=save_assoc, trace_bee=trace_bee, real_cluster_id_global=rcid_global, cathode_rescue_on=cathode_rescue, cathode_rescue_unmatched=cathode_rescue_unmatched);
            local per_apa_pre = [g.intern(
                innodes=[active_clusters[n], masked_clusters[n], opflash_sources[n]],
                centernodes=[clus_pipes[n]],
                outnodes=[flash_attach[n]],
                edges=[
                    g.edge(active_clusters[n], clus_pipes[n], 0, 0),   // port 0 = /live
                    g.edge(masked_clusters[n], clus_pipes[n], 0, 1),   // port 1 = /dead (2-view)
                    g.edge(clus_pipes[n], flash_attach[n], 0, 0),      // port 0 = pctree
                    g.edge(opflash_sources[n], flash_attach[n], 0, 1), // port 1 = opflash
                ]
            ) for n in std.range(0, nanodes - 1)];
            g.intern(
                innodes=per_apa_pre,
                centernodes=[jointql],
                outnodes=[clus_all],
                edges=[g.edge(per_apa_pre[i], jointql, 0, i) for i in std.range(0, nanodes - 1)]
                      + [g.edge(jointql, clus_all, 0, 0)]
            )
        else
            // active ─┐
            // masked ─┼─► clus.per_apa ─► FlashTensorToOpticalPCs ─► QLMatching ─► (PointTreeMerging -> MABC)
            // opflash ┘
            local per_apa = [g.intern(
                innodes=[active_clusters[n], masked_clusters[n], opflash_sources[n]],
                centernodes=[clus_pipes[n], flash_attach[n]],
                outnodes=[matching_pipes[n]],
                edges=[
                    g.edge(active_clusters[n], clus_pipes[n], 0, 0),   // port 0 = /live
                    g.edge(masked_clusters[n], clus_pipes[n], 0, 1),   // port 1 = /dead (2-view)
                    g.edge(clus_pipes[n], flash_attach[n], 0, 0),      // port 0 = pctree
                    g.edge(opflash_sources[n], flash_attach[n], 0, 1), // port 1 = opflash
                    g.edge(flash_attach[n], matching_pipes[n], 0, 0),
                ]
            ) for n in std.range(0, nanodes - 1)];
            local clus_all = clus_maker.all_apa(anodes, dump=true, tensor_outname=save_tensors, save_real_cluster_id=save_rcid, save_assoc_cluster_id=save_assoc, trace_bee=trace_bee, real_cluster_id_global=rcid_global, cathode_rescue_on=cathode_rescue, cathode_rescue_unmatched=cathode_rescue_unmatched);
            g.intern(
                innodes=per_apa,
                outnodes=[clus_all],
                edges=[g.edge(per_apa[i], clus_all, 0, i)
                       for i in std.range(0, nanodes - 1)]
            );

    local app = {
        type: 'Pgrapher',
        data: { edges: g.edges(graph) },
    };

    local cmdline = {
        type: 'wire-cell',
        data: {
            plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellAux', 'WireCellSio',
                      'WireCellSigProc', 'WireCellImg', 'WireCellRoot', 'WireCellTbb',
                      'WireCellClus', 'WireCellMatch'],
            apps: ['Pgrapher'],
        },
    };

    [cmdline] + g.uses(graph) + cathode_fv.configs + [app]
