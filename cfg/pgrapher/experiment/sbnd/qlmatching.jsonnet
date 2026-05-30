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
        },
    }, nin=1, nout=1),
}
