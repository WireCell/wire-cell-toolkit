// Canonical SBND charge-light (Q/L) matching helper.
//
// Single source of truth for the SBND matching graph nodes and the matching-only
// detector constants.  Mirrors clus.jsonnet / img.jsonnet: a factory function whose
// methods build typed pnodes.  Re-exported by the standalone dev chain via a thin
// sbnd_xin/qlmatching.jsonnet shim (like sbnd_xin/clus.jsonnet).
//
// Pipeline per APA:  TensorFileSource -> FlashTensorToOpticalPCs -> QLMatching
//   - TensorFileSource   reads the opflash archive into an opflash matrix tensor set.
//   - FlashTensorToOpticalPCs   (Aux fan-in) expands that matrix into the canonical
//                         "flash"/"light"/"flashlight" point clouds on the live root.
//   - QLMatching          reads the canonical flash PCs and writes back a per-cluster
//                         matched-flash scalar (so Cluster::get_flash() works).
//
// The matching constants (nchan, ch_mask) live HERE, not in params.jsonnet, on
// purpose: params.jsonnet is imported by many sim/production configs, and these are
// matching-only knobs.  drift_speed is taken from the caller's params so the common
// SBND value (and any DL/DT/lifetime overrides) flow through unchanged.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params) {
    // --- SBND matching constants (matching-only) ---
    local nchan = 312,
    // Dead PMT channels excluded from matching.  67/92/170/218/248 added after the
    // sbnd_xin saturation-PE study (sbnd_xin/docs/saturation-pe.md): they read 0 PE in
    // 100% of flashes in BOTH data and MC but were previously unmasked, so they only got
    // dropped in bright MC flashes by the (now-disabled) mc_saturation_pe gate.  Masking
    // them here handles them consistently in both modes and at all flash brightnesses.
    local ch_mask = [39, 64, 66, 67, 71, 85, 86, 87, 92, 115, 138, 141, 170, 197, 217,
                     218, 221, 222, 223, 226, 245, 248, 249, 302],

    // Per-PMT predicted-PE non-linearity overlay (maps each PMT's accumulated
    // predicted PE into the saturated/observed space; params fitted by
    // sbnd_xin/pmt_nonlinearity_curve.py --emit-qlmatching). Applied by default for
    // SBND (pmt_nl=true on matching()/matching_joint()); pass pmt_nl=false to disable.
    local nlp = import 'pgrapher/experiment/sbnd/pmt_nonlinearity_params.jsonnet',
    local nl_on = {
        pmt_nonlinearity: true,
        pmt_nl_knee: nlp.pmt_nl_knee,
        pmt_nl_beta: nlp.pmt_nl_beta,
        pmt_nl_gamma: nlp.pmt_nl_gamma,
    },

    // Opflash archive reader for APA n.
    opflash_source(n):: g.pnode({
        type: 'TensorFileSource',
        name: 'opflash_src_apa%d' % n,
        data: {
            inname: 'opflash_apa%d.tar.gz' % n,
            prefix: 'opflash_',
        },
    }, nin=0, nout=1),

    // Opflash matrix -> canonical flash/light/flashlight PCs (2->1 fan-in).
    // port 0 = cluster pctree, port 1 = opflash matrix.
    flash_attach(n):: g.pnode({
        type: 'FlashTensorToOpticalPCs',
        name: 'flash_attach_apa%d' % n,
        data: {
            nchan: nchan,
            // Add the per-frame "frame_apply_at_caf" (ns) offset from the opflash
            // tensorset metadata to every flash time (no-op if the key is absent).
            // Set false to read the raw, uncorrected flash time.
            correct_flash_time: true,
        },
    }, nin=2, nout=1),

    // DATA predicted-light scale (sim stays 1.0); see QtoL note below. 1.0 disables.
    local data_qtol = 0.86,

    // Common QLMatching `data` block (everything except the anode binding).
    // Single source of truth so the per-APA matching() and the joint
    // matching_joint() below cannot drift apart when these knobs are retuned.
    local match_data(dv, reality, semimodel_file, cathode_fiducial) = {
            detector_volumes: wc.tn(dv),
            // CPA structure-exclusion fiducial tn (cfg/.../cathode_fiducial.jsonnet);
            // '' => disabled, cathode-end flag_at_x_boundary uses the flat-cathode
            // 1-D test. SBND passes the cpa-exclusion CompositeFiducial tn here and
            // injects its component configs at the top level of the job.
            cathode_fiducial: cathode_fiducial,
            beamonly: false,
            data: if reality == 'data' then true else false,
            // DATA-only predicted-light scale (Q/L PE-error study, match/docs): the
            // data prediction over-predicts ~16%, so scale its light yield down (sim is
            // left at 1.0 -- it under-predicts). Set data_qtol=1.0 to disable. QtoL feeds
            // pred_flash = q*QtoL*vis*eff, so this is exactly a per-event prediction scale.
            QtoL: if reality == 'data' then data_qtol else 1.0,
            drift_speed: params.lar.drift_speed,
            nchan: nchan,  // must match FlashTensorToOpticalPCs.nchan (writer/reader coupling)
            ch_mask: ch_mask,
            flash_minPE: 50,
            semimodel_file: semimodel_file,

            // --- Matching tuning constants (see match/docs/improve_progress.md).
            // These were inline literals in QLMatching; surfaced here as the single
            // source of truth.  Values equal the C++ defaults, so behavior is
            // unchanged unless deliberately retuned. ---
            // §A active-volume cushions.  The raw active-volume bounds now come
            // from the DetectorVolumes service (m_dv->inner_bounds) rather than
            // hard-coded SBND literals; these signed cushions adjust the
            // effective PE-inclusion / boundary-flag windows in the per-TPC
            // anode->cathode drift coordinate u (outward-positive).  Defaults
            // follow the MicroBooNE prototype convention and shift the inclusion
            // window slightly vs the old ±2000/0-5000 literals (intended);
            // override them to recover the old bounds bit-identically.
            // (Validation: anode_ext1=1.45, cathode_ext1=0.45, y_cushion=-0.03,
            // z_cushion=0 reproduce the old [-200,0]/[0,200], |y|<=200, z[0,500]
            // gate bit-for-bit — used to confirm this refactor is physics-neutral.)
            anode_ext1: -2.0*wc.cm,
            anode_ext2: 4.0*wc.cm,
            cathode_ext1: 1.2*wc.cm,
            cathode_ext2: -2.0*wc.cm,
            y_cushion: 0.0*wc.cm,
            z_cushion: 0.0*wc.cm,
            // §D pre-selection / bad-match gates.
            // mc_saturation_pe DISABLED (set above any real flash PE).  The sbnd_xin
            // saturation-PE study (sbnd_xin/docs/saturation-pe.md) found no genuine PMT
            // saturation in MC: the zero-PE PMTs in bright flashes are dead channels
            // (now in ch_mask above), identical in MC and data, not simulated saturation.
            // The gate is MC-only, so it just introduced an MC/data asymmetry; retired.
            // (C++ default stays 5000; set 1e12 here to keep the gate inert.)
            mc_saturation_pe: 1e12,
            // §E out-of-beam QA cuts.
            outbeam_ks_max: 0.2,
            outbeam_chi2ndf_max: 20,
            outbeam_pe_frac: 0.5,
            // §C LASSO weights.
            lasso_lambda: 0.1,
            delta_charge: 0.01,
            delta_light: 0.025,
            delta_shape: 0.01,
            bkg_weight: 0.5,
            pe_mismatch_knee: 0.3,
            pe_mismatch_floor: 0.3,
            // §G per-PMT light-error model (Q/L PE-error study, match/docs): a constant
            // floor (PE) below the knee, fractional above -- larger fractional error at low
            // PE, ->frac at high PE. Applied to BOTH data and sim (conservative for sim, kept
            // consistent so cuts derived on data carry over). Old values: 0.3 / 0.3 / 1.0.
            pe_err_floor: 5.0,
            pe_err_frac: 0.25,
            pe_err_knee: 20.0,
            // Compute the bundle chi2's per-PMT error from the PREDICTED PE (not measured).
            // SBND-on; default false keeps the measured-based chi2 (other detectors, and
            // the LASSO solve, are unchanged -- the LASSO weight stays per-flash measured).
            pe_err_on_pred: true,
            flash_pe_threshold: 0.0,
            // §F bundle-quality thresholds.
            bundle_ks_merge_max: 0.2,
            bundle_chi2ndf_merge_max: 20,
            bundle_addmerge_exponent: 0.8,
            highconsist_ks_max: 0.06,
            highconsist_min_ndf: 3,
            bundle_pe_ndf_knee: 1.0,
            // Apply ch_mask / per-TPC / saturation masking to the KS shape metric
            // too (the chi2 and LASSO paths already mask). Without this, masked
            // channels keep their raw measured PE in the KS while the prediction is
            // 0 there, distorting the shape comparison (worst for the per-TPC mask,
            // whose opposite-TPC PMTs carry real PE). C++ default OFF; SBND-on.
            bundle_mask_ks: true,
            // §H raw readout-window truncation flag is always computed by
            // QLMatching (T0-independent, APA-agnostic) and is currently inert
            // (no consumer). edge threshold = 24 ticks (6 live slices, rebin 4).
            // readout_window_ticks is the EXCLUSIVE window end used for the
            // trailing-edge test: SBND daq.nticks is 3427, but with rebin 4 the
            // final 4-tick slice's slice_index_max reaches 3428, so 3428 is the
            // correct reference. A cluster is flagged truncated when its leading
            // slice is in [0,24] or its trailing slice is in [3404,3427]
            // (= 3428-24 .. window end).
            window_edge_ticks: 24,
            readout_window_ticks: 3428,
            // Discard any (flash, cluster) bundle whose cluster is not contained
            // in the TPC drift box once the flash T0 x-offset is applied — the
            // prototype flag_good_bundle gate (ToyMatching.cxx 272-275). See the
            // 4-part in-window guard in compute_endpoint_flags (match/docs
            // qlmatching-code.md §4.1a). Default OFF in C++; enabled here for SBND.
            require_containment: true,
    },

    // Charge-light matching for APA n.  `dv` is the DetectorVolumes node for this
    // anode (clus_maker.detector_volumes([anode])); it is emitted by the clustering
    // graph, here we only reference it by type:name.
    // `pmt_nl` (default true) bakes the per-PMT predicted-PE non-linearity overlay
    // (nl_on) into the node; pass pmt_nl=false to disable it. `extra` is an optional
    // data overlay merged last (default {} => no-op) for other per-call tweaks.
    matching(anode, dv, n, reality, semimodel_file, cathode_fiducial='', calib_dump='', pmt_nl=true, extra={}):: g.pnode({
        type: 'QLMatching',
        name: 'matching%d' % n,
        data: { anode: wc.tn(anode), calib_dump: calib_dump }
              + match_data(dv, reality, semimodel_file, cathode_fiducial)
              + (if pmt_nl then nl_on else {})
              + extra,
    }, nin=1, nout=1),

    // Joint multi-APA charge-light matching: ONE QLMatching node with one input
    // port per anode (multiplicity = #anodes).  It matches each APA independently
    // in its own isolated run (bit-identical to the per-APA matching() result) and
    // then merges the per-APA cluster trees into one output, reproducing the
    // standalone clus_all_apa PointTreeMerging it replaces.  `dv` is the all-anode
    // DetectorVolumes (clus_maker.detector_volumes(anodes)).  Same tuning as
    // matching(); adds the anodes list and the opflash root-PC concatenation.
    matching_joint(anodes, dv, reality, semimodel_file, cathode_fiducial='', calib_dump='', pmt_nl=true, extra={}):: g.pnode({
        type: 'QLMatching',
        name: 'matching_joint',
        data: {
            anodes: [wc.tn(a) for a in anodes],
            // Concatenate the per-APA optical-flash display PC into the merged
            // grouping (mirrors clus_all_apa PointTreeMerging.root_pcs_to_merge).
            root_pcs_to_merge: ['opflash'],
            // Hand-scan calibration dump path ('' => off, production-identical).
            // The joint node sees BOTH APAs in one operator() call, so a single
            // dump file holds both TPCs (sbnd_xin/ql_scan).
            calib_dump: calib_dump,
        } + match_data(dv, reality, semimodel_file, cathode_fiducial)
          + (if pmt_nl then nl_on else {})  // PMT non-linearity ON by default for SBND (pmt_nl=false disables)
          + extra,  // optional overlay (default {} => no-op) for other per-call tweaks
        // The all-anode DetectorVolumes is referenced only here (the per-APA path's
        // clustering pulls in the per-APA DVs; this all-anode one would otherwise be
        // dangling), so declare it as a dependency to get it into the job config.
    }, nin=std.length(anodes), nout=1, uses=[dv]),
}
