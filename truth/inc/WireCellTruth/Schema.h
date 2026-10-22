/** Schema of the native MC-truth products.

    One place for the names, column orders, sentinel values and class
    codes shared by every component of WireCellTruth and by the
    consumers (clus visitors, python readers).  Changing anything here
    bumps schema_version and is recorded in labels_meta.json.

    Design: wcp-porting-img pr_ml/docs/19 (sec 3.4-3.6).
 */

#ifndef WIRECELLTRUTH_SCHEMA
#define WIRECELLTRUTH_SCHEMA

#include <string>
#include <vector>

namespace WireCell::Truth {

    /// Bumped on any change of a name, order or sentinel below.
    constexpr int schema_version = 2;
    std::string schema_version_string();

    /// Sentinels of the per-point / per-blob track key "tid".
    constexpr int tid_ghost = -1;   ///< no true charge in the blob
    constexpr int tid_shared = -2;  ///< largest fraction below 0.5

    /// Semantic class of a particle (SPINE convention; shape
    /// precedence for shared points is track > michel > shower >
    /// delta > le).
    enum class Cls : int { track = 0, shower = 1, michel = 2, delta = 3, le = 4, unknown = -1 };
    std::string cls_name(Cls c);

    /// Column names, in order, of the 2-D double tensors.

    /// truth_nu: one row per generator interaction (SBND
    /// wclsTruthInformationAttacher columns, in order).
    const std::vector<std::string>& truth_nu_columns();

    /// truth_pf: one row per G4 particle kept.  The SBND columns come
    /// first and in order; the native columns follow.  Readers index
    /// by md["columns"].
    const std::vector<std::string>& truth_pf_columns();

    /// blobtruth: sparse (blob, track) charge table.
    const std::vector<std::string>& blobtruth_columns();

    /// blobtruth_sum: one row per blob.
    const std::vector<std::string>& blobtruth_sum_columns();

    /// Keys of the per-blob-node local point cloud named "truth"
    /// (same size and order as the "3d" cloud).  The key set is
    /// fixed so that every node is key-homogeneous.
    const std::vector<std::string>& truth_pc_keys();

    /// Name of that point cloud.
    constexpr const char* truth_pc_name = "truth";

    /// Datatype strings of the aux tensors.
    constexpr const char* truth_nu_datatype = "truth_nu";
    constexpr const char* truth_pf_datatype = "truth_pf";
    constexpr const char* blobtruth_datatype = "blobtruth";

    /// Members of the depo file written by DepoFileSink when
    /// truth_sidecar is on: depo_truth_<ident>.npy columns.
    const std::vector<std::string>& depo_truth_columns();

    /// Index of a column name in a column list, or -1.
    int column_index(const std::vector<std::string>& cols, const std::string& name);

    /// Geant4 creation-process names, in code order.  The "process"
    /// column of truth_pf is the index into this list; an unlisted
    /// process gets process_unknown (0).  python/truth_table.py
    /// carries the same list (checked by the dumper against the
    /// "process_names" member it writes).
    const std::vector<std::string>& process_names();
    constexpr int process_unknown = 0;
    int process_code(const std::string& name);

    /// Columns of the sed_ids table of the truth file (one row per
    /// SimEnergyDeposit, in vector order): the depo id in Stage A
    /// index mode is the row number.
    const std::vector<std::string>& sed_ids_columns();

    /// Members of the truth file (truth_<ident>.npz) written by
    /// python/dump_truth_table.py and read by TruthTable.
    constexpr const char* truth_file_nu = "truth_nu";
    constexpr const char* truth_file_pf = "truth_pf";
    constexpr const char* truth_file_sed_ids = "sed_ids";
    constexpr const char* truth_file_sed_trackq = "sed_trackq";
}

#endif
