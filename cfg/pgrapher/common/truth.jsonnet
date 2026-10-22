// Node builders of the WireCellTruth plugin: native MC-truth plumbing
// (wcp-porting-img pr_ml/docs/19).  Every builder takes an `extra` object
// merged into the component's data so that knobs are injected where the
// inode literal is built (M6).  Nothing here is wired into a production job
// by default; the detector jobs add these nodes behind default-OFF
// arguments with the key-suppression idiom.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";

{
    plugin: "WireCellTruth",

    /// BlobDepoFill with per-track bookkeeping, same interface (ICluster,
    /// IDepoSet) -> ICluster.  Drop-in for the tru0/tru catchers; bit-identical
    /// output with the smearing knobs off.  speed is required (as BlobDepoFill).
    blob_depo_label(name, speed, time_offset=0.0, extra={}) :: g.pnode({
        type: "BlobDepoLabel",
        name: name,
        data: {
            speed: speed,
            time_offset: time_offset,
        } + extra,
    }, nin=2, nout=1),

    /// The per-(blob, track) charge table as tensors: (ICluster, IDepoSet) -> ITensorSet.
    /// datapath pattern "%d" = cluster ident; with_cluster adds the filled cluster.
    blob_depo_label_tensor(name, speed, time_offset=0.0, datapath="clusters/%d", with_cluster=true, extra={}) :: g.pnode({
        type: "BlobDepoLabelTensor",
        name: name,
        data: {
            speed: speed,
            time_offset: time_offset,
            datapath: datapath,
            with_cluster: with_cluster,
        } + extra,
    }, nin=2, nout=1),

    /// The truth_nu / truth_pf tensors of truth files written by
    /// truth/python/dump_truth_table.py: a source node.  inname "%d" = ident.
    truth_table_source(name, inname, idents=[0], prefix="truth/%d", extra={}) :: g.pnode({
        type: "TruthTableSource",
        name: name,
        data: {
            inname: inname,
            idents: idents,
            prefix: prefix,
        } + extra,
    }, nin=0, nout=1),

    /// The per-point "truth" point cloud on a live point tree:
    /// (ITensorSet from PointTreeBuilding, ICluster, IDepoSet) -> ITensorSet.
    /// anode is the AnodePlane object, face the face index of the tree,
    /// detector_volumes the DetectorVolumes object (both as PointTreeBuilding);
    /// speed/time_offset are the BlobDepoFill association knobs (the truth
    /// tier's values); truth_table "" = particle table from the depos.
    point_labeler(name, anode, face, detector_volumes, speed, time_offset=0.0,
                  inpath="pointtrees/%d/live", truth_table="", extra={}) :: g.pnode({
        type: "TruthPointLabeler",
        name: name,
        data: {
            anode: wc.tn(anode),
            face: face,
            detector_volumes: wc.tn(detector_volumes),
            speed: speed,
            time_offset: time_offset,
            inpath: inpath,
            [if truth_table != "" then "truth_table"]: truth_table,
        } + extra,
    }, nin=3, nout=1, uses=[anode, detector_volumes]),
}
