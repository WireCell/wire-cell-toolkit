#include "WireCellTruth/Schema.h"

using namespace WireCell;

std::string Truth::schema_version_string()
{
    return "wct-truth-schema-" + std::to_string(schema_version);
}

std::string Truth::cls_name(Cls c)
{
    switch (c) {
    case Cls::track: return "track";
    case Cls::shower: return "shower";
    case Cls::michel: return "michel";
    case Cls::delta: return "delta";
    case Cls::le: return "le";
    default: return "unknown";
    }
}

const std::vector<std::string>& Truth::truth_nu_columns()
{
    static const std::vector<std::string> cols = {
        "nu_idx", "pdg", "ccnc", "mode", "int_type", "flavor",
        "E", "vtx_x", "vtx_y", "vtx_z", "t", "edep"};
    return cols;
}

const std::vector<std::string>& Truth::truth_pf_columns()
{
    static const std::vector<std::string> cols = {
        // SBND attacher columns, in order
        "nu_row", "trackid", "parent_trackid", "mother_trackid", "pdg", "process",
        "E", "KE",
        "start_x", "start_y", "start_z", "start_t",
        "end_x", "end_y", "end_z", "end_t",
        "start_px", "start_py", "start_pz",
        // native additions
        "orig_trackid", "gen_label", "pf_node", "cls", "merged_into", "nvis_q"};
    return cols;
}

const std::vector<std::string>& Truth::blobtruth_columns()
{
    static const std::vector<std::string> cols = {
        "blob_ident", "face", "slice_ident", "trackid", "q_bk", "e_bk"};
    return cols;
}

const std::vector<std::string>& Truth::blobtruth_sum_columns()
{
    static const std::vector<std::string> cols = {
        "blob_ident", "q_true", "q_top", "tid_top", "ntracks"};
    return cols;
}

const std::vector<std::string>& Truth::truth_pc_keys()
{
    static const std::vector<std::string> keys = {
        "tid", "orig_tid", "pid_frac", "true_q_frac", "instance_valid",
        "pdg", "cls", "iid", "pf_node",
        "d_traj", "dir_x", "dir_y", "dir_z", "s_t",
        "v_heat_vtx", "v_heat_end", "v_heat_kink", "v_heat_conv",
        "end", "dedx_true"};
    return keys;
}

const std::vector<std::string>& Truth::depo_truth_columns()
{
    static const std::vector<std::string> cols = {"energy", "orig_id", "reserved"};
    return cols;
}

int Truth::column_index(const std::vector<std::string>& cols, const std::string& name)
{
    for (size_t i = 0; i < cols.size(); ++i) {
        if (cols[i] == name) return (int)i;
    }
    return -1;
}

const std::vector<std::string>& Truth::process_names()
{
    // Index = code.  Append only; never reorder (schema_version).
    static const std::vector<std::string> names = {
        "unknown", "primary",
        "eIoni", "muIoni", "hIoni", "ionIoni",
        "eBrem", "muBrems", "hBrems",
        "annihil", "conv", "compt", "phot", "Rayl",
        "Decay", "muMinusCaptureAtRest", "muMinusAtomicCapture",
        "hadElastic", "CoulombScat", "msc",
        "neutronInelastic", "protonInelastic", "pi+Inelastic", "pi-Inelastic",
        "kaon+Inelastic", "kaon-Inelastic", "kaon0LInelastic", "kaon0SInelastic",
        "lambdaInelastic", "sigma+Inelastic", "sigma-Inelastic",
        "anti_protonInelastic", "anti_neutronInelastic",
        "dInelastic", "tInelastic", "He3Inelastic", "alphaInelastic", "ionInelastic",
        "nCapture", "hFritiofCaptureAtRest", "hBertiniCaptureAtRest",
        "muPairProd", "hPairProd", "photonNuclear", "electronNuclear", "muonNuclear",
        "RadioactiveDecay", "RadioactiveDecayBase", "Radioactivation",
        "Transportation", "StepLimiter"};
    return names;
}

int Truth::process_code(const std::string& name)
{
    const int ind = column_index(process_names(), name);
    return ind < 0 ? process_unknown : ind;
}

const std::vector<std::string>& Truth::sed_ids_columns()
{
    static const std::vector<std::string> cols = {"trackid", "orig_trackid"};
    return cols;
}
