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
            mc_saturation_pe: 5000,
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
    }, nin=1, nout=1),
}
