// Comparison app: prototype T_tagger/T_kine vs toolkit T_tagger/T_kine.
// Opens both ROOT files, walks all scalar and vector branches, groups diffs
// by originating tagger function, and writes per-category histograms + a
// T_summary tree to report.root.  Exit code 0 means all categories agree.
//
// Usage:
//   wire-cell-uboone-tagger-compare -p<proto.root> -t<toolkit.root>
//       [-o<report.root>] [-v]
//
// Matching is positional (entry 0 vs 0, etc.).

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"

// ============================================================
//  Category table
// ============================================================
struct Category {
    std::string name;
    std::vector<std::string> prefixes;
};

static std::vector<Category> make_categories()
{
    return {
        // nu_x/y/z are exact names (no further suffix), so prefix match works safely.
        {"nu_vertex",  {"nu_x", "nu_y", "nu_z"}},
        // cosmic tagger flags and sub-taggers (cosmict_2..10 also start with "cosmict_")
        {"cosmic",     {"cosmic_", "cosmict_"}},
        {"gap",        {"gap_"}},
        // mip_quality_ must come before mip_ since mip_quality_ ⊂ mip_
        {"mip",        {"mip_quality_", "mip_"}},
        {"ssm",        {"ssm_"}},
        {"shw_sp",     {"shw_sp_"}},
        {"stem",       {"stem_"}},
        // pi0 and related: pio_, sig_, mgo_, mgt_, stw_, spt_ all live in the pi0/shower-pi0 family
        {"pio_family", {"pio_", "sig_", "mgo_", "mgt_", "stw_", "spt_"}},
        {"lem",        {"lem_"}},
        // bad-reconstruction taggers share a br_ parent
        {"br",         {"brm_", "cme_", "anc_", "br_", "br1_", "br2_", "br3_", "br4_"}},
        {"tro",        {"tro_"}},
        {"hol_lol",    {"hol_", "lol_"}},
        {"vis",        {"vis_"}},
        {"numu",       {"numu_"}},
        {"nue",        {"nue_"}},
        {"match",      {"match_"}},
        // kine_ covers all T_kine branches (comparison loop uses them separately)
        {"kine",       {"kine_"}},
    };
}

static std::string classify_branch(const std::string& bname,
                                   const std::vector<Category>& cats)
{
    for (const auto& cat : cats) {
        for (const auto& pfx : cat.prefixes) {
            if (bname.size() >= pfx.size() &&
                bname.substr(0, pfx.size()) == pfx) {
                return cat.name;
            }
        }
    }
    return "_uncategorized";
}

// ============================================================
//  Branch-level type detection
// ============================================================
enum BrType { SCALAR_F, SCALAR_I, VEC_F, VEC_I, SKIP };

static BrType detect_type(TBranch* br)
{
    const char* cls = br->GetClassName();
    if (cls && std::strlen(cls) > 0) {
        std::string cn(cls);
        if (cn.find("vector<float>") != std::string::npos)  return VEC_F;
        if (cn.find("vector<int>")   != std::string::npos)  return VEC_I;
        return SKIP;  // unsupported object type
    }
    TObjArray* leaves = br->GetListOfLeaves();
    if (!leaves || leaves->GetEntries() == 0) return SKIP;
    TLeaf* leaf = (TLeaf*)leaves->At(0);
    std::string tn(leaf->GetTypeName());
    if (tn == "Float_t" || tn == "float") return SCALAR_F;
    if (tn == "Int_t"   || tn == "int")   return SCALAR_I;
    return SKIP;
}

// ============================================================
//  Per-branch comparison state
// ============================================================
struct BranchStat {
    std::string name;
    std::string category;
    BrType      type = SKIP;
    bool        in_kine = false;  // true for T_kine branches

    // scalar buffers (addresses remain stable because BranchStat lives in a list)
    float proto_f = 0, tool_f = 0;
    int   proto_i = 0, tool_i = 0;
    // object buffers (ROOT sets these pointers on first GetEntry)
    std::vector<float>* proto_vf = nullptr;
    std::vector<float>* tool_vf  = nullptr;
    std::vector<int>*   proto_vi = nullptr;
    std::vector<int>*   tool_vi  = nullptr;

    TBranch* proto_br = nullptr;
    TBranch* tool_br  = nullptr;

    // cumulative stats
    long   n_compared  = 0;
    long   n_diff      = 0;
    long   n_sentinel  = 0;
    double max_abs_diff = 0;
    double sum_abs_diff = 0;

    // histogram of normalized diff (toolkit − proto) / scale
    TH1F* h_norm_diff = nullptr;
};

static void setup_branch_stat(BranchStat& bs, TDirectory* dir)
{
    BrType t = detect_type(bs.proto_br);
    // must agree with toolkit branch type
    if (t != detect_type(bs.tool_br)) {
        bs.type = SKIP;
        return;
    }
    bs.type = t;

    switch (t) {
    case SCALAR_F:
        bs.proto_br->SetAddress(&bs.proto_f);
        bs.tool_br->SetAddress(&bs.tool_f);
        break;
    case SCALAR_I:
        bs.proto_br->SetAddress(&bs.proto_i);
        bs.tool_br->SetAddress(&bs.tool_i);
        break;
    case VEC_F:
        bs.proto_br->SetAddress(&bs.proto_vf);
        bs.tool_br->SetAddress(&bs.tool_vf);
        break;
    case VEC_I:
        bs.proto_br->SetAddress(&bs.proto_vi);
        bs.tool_br->SetAddress(&bs.tool_vi);
        break;
    default:
        return;
    }

    if (dir) {
        dir->cd();
        std::string hname = "h_" + bs.name + "_ndiff";
        bs.h_norm_diff = new TH1F(
            hname.c_str(),
            (bs.name + ";(toolkit-proto)/scale;entries").c_str(),
            100, -2.0, 2.0);
    }
}

// ============================================================
//  Compare one event for one branch
// ============================================================
static void compare_one(BranchStat& bs, long proto_entry, long tool_entry)
{
    if (!bs.proto_br || !bs.tool_br || bs.type == SKIP) return;
    bs.proto_br->GetEntry(proto_entry);
    bs.tool_br->GetEntry(tool_entry);

    constexpr float  SENTINEL_F = -999.0f;
    constexpr double EPS        = 1e-6;

    switch (bs.type) {
    case SCALAR_F: {
        bs.n_compared++;
        if (bs.proto_f == SENTINEL_F || bs.tool_f == SENTINEL_F) {
            bs.n_sentinel++;
            return;
        }
        double diff  = bs.tool_f - bs.proto_f;
        double scale = std::max({std::abs((double)bs.proto_f),
                                 std::abs((double)bs.tool_f), EPS});
        double nd = diff / scale;
        if (std::abs(nd) > 1e-3) {
            bs.n_diff++;
            bs.max_abs_diff  = std::max(bs.max_abs_diff, std::abs(diff));
            bs.sum_abs_diff += std::abs(diff);
        }
        if (bs.h_norm_diff) bs.h_norm_diff->Fill(nd);
        break;
    }
    case SCALAR_I: {
        bs.n_compared++;
        if (bs.proto_i != bs.tool_i) {
            bs.n_diff++;
            bs.max_abs_diff  = std::max(bs.max_abs_diff,
                                        std::abs((double)(bs.tool_i - bs.proto_i)));
            bs.sum_abs_diff += std::abs((double)(bs.tool_i - bs.proto_i));
        }
        break;
    }
    case VEC_F: {
        if (!bs.proto_vf || !bs.tool_vf) return;
        bs.n_compared++;
        if (bs.proto_vf->size() != bs.tool_vf->size()) {
            bs.n_diff++;
            return;
        }
        bool event_differs = false;
        for (size_t k = 0; k < bs.proto_vf->size(); ++k) {
            float pv = (*bs.proto_vf)[k];
            float tv = (*bs.tool_vf)[k];
            if (pv == SENTINEL_F || tv == SENTINEL_F) { bs.n_sentinel++; continue; }
            double diff  = tv - pv;
            double scale = std::max({std::abs((double)pv), std::abs((double)tv), EPS});
            double nd = diff / scale;
            if (std::abs(nd) > 1e-3) {
                event_differs = true;
                bs.max_abs_diff  = std::max(bs.max_abs_diff, std::abs(diff));
                bs.sum_abs_diff += std::abs(diff);
                if (bs.h_norm_diff) bs.h_norm_diff->Fill(nd);
            }
        }
        if (event_differs) bs.n_diff++;
        break;
    }
    case VEC_I: {
        if (!bs.proto_vi || !bs.tool_vi) return;
        bs.n_compared++;
        if (bs.proto_vi->size() != bs.tool_vi->size()) {
            bs.n_diff++;
            return;
        }
        for (size_t k = 0; k < bs.proto_vi->size(); ++k) {
            if ((*bs.proto_vi)[k] != (*bs.tool_vi)[k]) {
                bs.n_diff++;
                bs.max_abs_diff = std::max(bs.max_abs_diff,
                    std::abs((double)((*bs.tool_vi)[k] - (*bs.proto_vi)[k])));
                bs.sum_abs_diff += std::abs((double)((*bs.tool_vi)[k] - (*bs.proto_vi)[k]));
                break;
            }
        }
        break;
    }
    default: break;
    }
}

// ============================================================
//  Build branch stats for one tree pair
// ============================================================
static void build_stats(TTree* ptree, TTree* ttree, bool is_kine,
                        const std::vector<Category>& cats,
                        std::map<std::string, TDirectory*>& cat_dirs,
                        TFile* out_file,
                        std::list<BranchStat>& stat_list)
{
    TObjArray* proto_branches = ptree->GetListOfBranches();
    for (int bi = 0; bi < proto_branches->GetEntries(); ++bi) {
        TBranch* pb = (TBranch*)proto_branches->At(bi);
        std::string bname(pb->GetName());
        TBranch* tb = ttree->GetBranch(bname.c_str());
        if (!tb) {
            std::cout << "[MISSING in toolkit] " << bname << "\n";
            continue;
        }

        std::string cat = classify_branch(bname, cats);
        if (cat_dirs.find(cat) == cat_dirs.end()) {
            TDirectory* d = out_file->mkdir(cat.c_str());
            cat_dirs[cat] = d;
        }

        stat_list.push_back(BranchStat{});   // stable address in list
        BranchStat& bs  = stat_list.back();
        bs.name         = bname;
        bs.category     = cat;
        bs.in_kine      = is_kine;
        bs.proto_br     = pb;
        bs.tool_br      = tb;
        setup_branch_stat(bs, cat_dirs[cat]);

        if (bs.type == SKIP)
            std::cout << "[SKIP unsupported type] " << bname << "\n";
    }

    // Report branches only in toolkit
    TObjArray* tool_branches = ttree->GetListOfBranches();
    for (int bi = 0; bi < tool_branches->GetEntries(); ++bi) {
        TBranch* tb = (TBranch*)tool_branches->At(bi);
        std::string bname(tb->GetName());
        if (!ptree->GetBranch(bname.c_str()))
            std::cout << "[ONLY in toolkit] " << bname << "\n";
    }
}

// ============================================================
//  Write T_summary tree into category directory
// ============================================================
static void write_category_summary(const std::string& cat,
                                   const std::list<BranchStat>& stats,
                                   TDirectory* dir)
{
    if (!dir) return;
    dir->cd();
    TTree* ts = new TTree("T_summary", "per-branch comparison summary");
    char branch_name[256];
    long long n_compared, n_diff, n_sentinel;
    double max_abs_diff, mean_abs_diff, frac_diff;
    ts->Branch("branch_name",  branch_name,     "branch_name/C");
    ts->Branch("n_compared",   &n_compared,     "n_compared/L");
    ts->Branch("n_diff",       &n_diff,         "n_diff/L");
    ts->Branch("n_sentinel",   &n_sentinel,     "n_sentinel/L");
    ts->Branch("max_abs_diff", &max_abs_diff,   "max_abs_diff/D");
    ts->Branch("mean_abs_diff",&mean_abs_diff,  "mean_abs_diff/D");
    ts->Branch("frac_diff",    &frac_diff,      "frac_diff/D");

    for (const auto& bs : stats) {
        if (bs.category != cat || bs.type == SKIP) continue;
        std::strncpy(branch_name, bs.name.c_str(), 255);
        branch_name[255] = '\0';
        n_compared  = bs.n_compared;
        n_diff      = bs.n_diff;
        n_sentinel  = bs.n_sentinel;
        max_abs_diff  = bs.max_abs_diff;
        mean_abs_diff = (bs.n_diff > 0) ? bs.sum_abs_diff / bs.n_diff : 0.0;
        frac_diff     = (bs.n_compared > 0) ? (double)bs.n_diff / bs.n_compared : 0.0;
        ts->Fill();
    }
    ts->Write();
}

// ============================================================
//  main
// ============================================================
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: wire-cell-uboone-tagger-compare"
                     " -p<proto.root> -t<toolkit.root>"
                     " [-o<report.root>] [-v]\n";
        return 1;
    }

    TString proto_filename;
    TString toolkit_filename;
    TString out_filename = "tagger_compare_report.root";
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-' || argv[i][1] == '\0') continue;
        switch (argv[i][1]) {
        case 'p': proto_filename   = &argv[i][2]; break;
        case 't': toolkit_filename = &argv[i][2]; break;
        case 'o': out_filename     = &argv[i][2]; break;
        case 'v': verbose = true; break;
        }
    }

    if (proto_filename.IsNull() || toolkit_filename.IsNull()) {
        std::cerr << "Need -p<proto.root> and -t<toolkit.root>\n";
        return 1;
    }

    TFile* proto_file = TFile::Open(proto_filename, "READ");
    TFile* tool_file  = TFile::Open(toolkit_filename, "READ");
    if (!proto_file || proto_file->IsZombie()) {
        std::cerr << "Cannot open prototype file: " << proto_filename << "\n"; return 1;
    }
    if (!tool_file || tool_file->IsZombie()) {
        std::cerr << "Cannot open toolkit file: " << toolkit_filename << "\n"; return 1;
    }

    auto* proto_tagger = (TTree*)proto_file->Get("T_tagger");
    auto* proto_kine   = (TTree*)proto_file->Get("T_kine");
    auto* tool_tagger  = (TTree*)tool_file->Get("T_tagger");
    auto* tool_kine    = (TTree*)tool_file->Get("T_kine");

    for (auto* t : {proto_tagger, proto_kine, tool_tagger, tool_kine}) {
        if (!t) {
            std::cerr << "Missing T_tagger or T_kine in one of the files\n";
            return 1;
        }
    }

    long N_tagger = (long)std::min(proto_tagger->GetEntries(),
                                   tool_tagger->GetEntries());
    long N_kine   = (long)std::min(proto_kine->GetEntries(),
                                   tool_kine->GetEntries());
    if (N_tagger <= 0) {
        std::cerr << "T_tagger has no entries\n"; return 1;
    }
    std::cout << "Comparing " << N_tagger << " T_tagger events, "
              << N_kine << " T_kine events\n";

    auto cats = make_categories();

    TFile* out_file = TFile::Open(out_filename, "RECREATE");
    if (!out_file || out_file->IsZombie()) {
        std::cerr << "Cannot create report file: " << out_filename << "\n"; return 1;
    }

    std::map<std::string, TDirectory*> cat_dirs;
    std::list<BranchStat> all_stats;

    build_stats(proto_tagger, tool_tagger, false, cats, cat_dirs, out_file, all_stats);
    build_stats(proto_kine,   tool_kine,   true,  cats, cat_dirs, out_file, all_stats);

    // Event loop — T_tagger
    for (long ev = 0; ev < N_tagger; ++ev) {
        for (auto& bs : all_stats) {
            if (!bs.in_kine) compare_one(bs, ev, ev);
        }
    }
    // Event loop — T_kine
    for (long ev = 0; ev < N_kine; ++ev) {
        for (auto& bs : all_stats) {
            if (bs.in_kine) compare_one(bs, ev, ev);
        }
    }

    // Write histograms and summary trees
    for (const auto& kv : cat_dirs) {
        const std::string& cat = kv.first;
        TDirectory* dir = kv.second;
        dir->cd();
        // Write all histograms already owned by the directory
        for (auto& bs : all_stats) {
            if (bs.category == cat && bs.h_norm_diff) {
                bs.h_norm_diff->Write();
            }
        }
        write_category_summary(cat, all_stats, dir);
    }

    out_file->Write();
    out_file->Close();

    // ---- Terminal summary table ----
    // Collect per-category aggregate stats
    struct CatSummary {
        long   n_branches      = 0;
        long   n_diff_branches = 0;
        long   n_events        = 0;
        double max_abs_diff    = 0;
        std::string worst_branch;
    };
    std::map<std::string, CatSummary> cat_summary;
    for (const auto& bs : all_stats) {
        if (bs.type == SKIP) continue;
        CatSummary& cs = cat_summary[bs.category];
        cs.n_branches++;
        cs.n_events = std::max(cs.n_events, bs.n_compared);
        if (bs.n_diff > 0) {
            cs.n_diff_branches++;
            if (bs.max_abs_diff > cs.max_abs_diff) {
                cs.max_abs_diff = bs.max_abs_diff;
                cs.worst_branch = bs.name;
            }
        }
    }

    std::cout << "\n";
    std::cout << std::left
              << std::setw(16) << "category"
              << std::setw(11) << "n_branches"
              << std::setw(10) << "n_events"
              << std::setw(17) << "n_diff_branches"
              << std::setw(32) << "worst_branch"
              << std::setw(12) << "max|diff|"
              << "\n";
    std::cout << std::string(95, '-') << "\n";

    bool any_diff = false;
    // Print in category order defined by make_categories(), then uncategorized
    std::vector<std::string> ordered_cats;
    for (const auto& cat : cats) ordered_cats.push_back(cat.name);
    ordered_cats.push_back("_uncategorized");

    for (const auto& cat : ordered_cats) {
        auto it = cat_summary.find(cat);
        if (it == cat_summary.end()) continue;
        const CatSummary& cs = it->second;
        if (cs.n_diff_branches > 0) any_diff = true;
        std::cout << std::left
                  << std::setw(16) << cat
                  << std::setw(11) << cs.n_branches
                  << std::setw(10) << cs.n_events
                  << std::setw(17) << cs.n_diff_branches
                  << std::setw(32) << cs.worst_branch
                  << std::setw(12) << cs.max_abs_diff
                  << "\n";
    }
    std::cout << "\nReport written to: " << out_filename << "\n";

    if (verbose) {
        std::cout << "\n--- Per-branch details (differing only) ---\n";
        for (const auto& bs : all_stats) {
            if (bs.n_diff == 0 || bs.type == SKIP) continue;
            std::cout << "  " << std::left << std::setw(48) << bs.name
                      << "  n_diff=" << bs.n_diff
                      << "/" << bs.n_compared
                      << "  max|diff|=" << bs.max_abs_diff
                      << "\n";
        }
    }

    proto_file->Close();
    tool_file->Close();

    return any_diff ? 1 : 0;
}
