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
    local ch_mask = [39, 64, 66, 71, 85, 86, 87, 115, 138, 141, 197, 217, 221,
                     222, 223, 226, 245, 249, 302],

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
        },
    }, nin=2, nout=1),

    // Charge-light matching for APA n.  `dv` is the DetectorVolumes node for this
    // anode (clus_maker.detector_volumes([anode])); it is emitted by the clustering
    // graph, here we only reference it by type:name.
    matching(anode, dv, n, reality, semimodel_file):: g.pnode({
        type: 'QLMatching',
        name: 'matching%d' % n,
        data: {
            anode: wc.tn(anode),
            detector_volumes: wc.tn(dv),
            beamonly: false,
            data: if reality == 'data' then true else false,
            QtoL: 1.0,
            drift_speed: params.lar.drift_speed,
            nchan: nchan,  // must match FlashTensorToOpticalPCs.nchan (writer/reader coupling)
            ch_mask: ch_mask,
            flash_minPE: 50,
            semimodel_file: semimodel_file,

            // --- Matching tuning constants (see match/docs/improve_progress.md).
            // These were inline literals in QLMatching; surfaced here as the single
            // source of truth.  Values equal the C++ defaults, so behavior is
            // unchanged unless deliberately retuned. ---
            // §A active-volume / drift bounds (cathode seam stays the origin x=0).
            x_bound: 2000*wc.mm,
            y_bound: 2000*wc.mm,
            z_min: 0*wc.mm,
            z_max: 5000*wc.mm,
            pmt_dist: 1950*wc.mm,
            // §D pre-selection / bad-match gates.
            mc_saturation_pe: 5000,
            drift_out_frac: 0.25,
            min_pred_pe: 10,
            preselect_chi2ndf_max: 1e4,
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
            // §G flash PE-error model.
            pe_err_floor: 0.3,
            pe_err_frac: 0.3,
            pe_err_knee: 1.0,
            flash_pe_threshold: 0.0,
            // §F bundle-quality thresholds.
            bundle_ks_merge_max: 0.2,
            bundle_chi2ndf_merge_max: 20,
            bundle_addmerge_exponent: 0.8,
            highconsist_ks_max: 0.06,
            highconsist_min_ndf: 3,
            bundle_pe_ndf_knee: 1.0,
        },
    }, nin=1, nout=1),
}
