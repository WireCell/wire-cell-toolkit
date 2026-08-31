// Comparison app: prototype T_tagger/T_kine vs toolkit T_tagger/T_kine.
// Opens both ROOT files, walks all scalar and vector branches, groups diffs
// by originating tagger function, and prints results to stdout.
// Exit code 0 means all categories agree.
//
// Usage:
//   wire-cell-uboone-tagger-compare -p <proto.root> -t <toolkit.root> [-v]
//                                   [-n <NeutrinoTaggerInfo.h>]
//
// Matching is positional (entry 0 vs 0, etc.).

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "TBranch.h"
#include "TFile.h"
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
enum BrType { SCALAR_F, SCALAR_D, SCALAR_I, VEC_F, VEC_I, SKIP };

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
    if (tn == "Float_t"  || tn == "float")  return SCALAR_F;
    if (tn == "Double_t" || tn == "double") return SCALAR_D;
    if (tn == "Int_t"    || tn == "int")    return SCALAR_I;
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
    float  proto_f = 0, tool_f = 0;
    double proto_d = 0, tool_d = 0;
    int    proto_i = 0, tool_i = 0;
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

    // values at first differing event (printed in verbose mode)
    float  first_diff_proto_f{0}, first_diff_tool_f{0};
    double first_diff_proto_d{0}, first_diff_tool_d{0};
    int    first_diff_proto_i{0}, first_diff_tool_i{0};
    std::vector<float> first_diff_proto_vf, first_diff_tool_vf;
    std::vector<int>   first_diff_proto_vi, first_diff_tool_vi;
};

static void setup_branch_stat(BranchStat& bs)
{
    BrType pt = detect_type(bs.proto_br);
    BrType tt = detect_type(bs.tool_br);

    // Allow Double_t ↔ Float_t mixed pairs: read both via their native buffers,
    // compare_one promotes both to double when computing the diff.
    bool mixed_df = (pt == SCALAR_D && tt == SCALAR_F) ||
                    (pt == SCALAR_F && tt == SCALAR_D);
    if (!mixed_df && pt != tt) {
        bs.type = SKIP;
        return;
    }
    // Canonical type: prefer SCALAR_D if either side is double.
    bs.type = (pt == SCALAR_D || tt == SCALAR_D) ? SCALAR_D : pt;

    switch (pt) {
    case SCALAR_F: bs.proto_br->SetAddress(&bs.proto_f); break;
    case SCALAR_D: bs.proto_br->SetAddress(&bs.proto_d); break;
    case SCALAR_I: bs.proto_br->SetAddress(&bs.proto_i); break;
    case VEC_F:    bs.proto_br->SetAddress(&bs.proto_vf); break;
    case VEC_I:    bs.proto_br->SetAddress(&bs.proto_vi); break;
    default: return;
    }
    switch (tt) {
    case SCALAR_F: bs.tool_br->SetAddress(&bs.tool_f); break;
    case SCALAR_D: bs.tool_br->SetAddress(&bs.tool_d); break;
    case SCALAR_I: bs.tool_br->SetAddress(&bs.tool_i); break;
    case VEC_F:    bs.tool_br->SetAddress(&bs.tool_vf); break;
    case VEC_I:    bs.tool_br->SetAddress(&bs.tool_vi); break;
    default: return;
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
    constexpr double SENTINEL_D = -999.0;
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
            if (bs.n_diff == 0) {
                bs.first_diff_proto_f = bs.proto_f;
                bs.first_diff_tool_f  = bs.tool_f;
            }
            bs.n_diff++;
            bs.max_abs_diff  = std::max(bs.max_abs_diff, std::abs(diff));
            bs.sum_abs_diff += std::abs(diff);
        }
        break;
    }
    case SCALAR_D: {
        // Covers pure Double_t/Double_t and mixed Double_t/Float_t pairs.
        // Each side is read into its native buffer; promote to double here.
        bs.n_compared++;
        double pv = (bs.proto_br->GetLeaf(bs.name.c_str()) &&
                     std::string(bs.proto_br->GetLeaf(bs.name.c_str())->GetTypeName()).find("Double") != std::string::npos)
                    ? bs.proto_d : (double)bs.proto_f;
        double tv = (bs.tool_br->GetLeaf(bs.name.c_str()) &&
                     std::string(bs.tool_br->GetLeaf(bs.name.c_str())->GetTypeName()).find("Double") != std::string::npos)
                    ? bs.tool_d : (double)bs.tool_f;
        if (pv == SENTINEL_D || tv == SENTINEL_D) { bs.n_sentinel++; return; }
        double diff  = tv - pv;
        double scale = std::max({std::abs(pv), std::abs(tv), EPS});
        double nd    = diff / scale;
        if (std::abs(nd) > 1e-3) {
            if (bs.n_diff == 0) {
                bs.first_diff_proto_d = pv;
                bs.first_diff_tool_d  = tv;
            }
            bs.n_diff++;
            bs.max_abs_diff  = std::max(bs.max_abs_diff, std::abs(diff));
            bs.sum_abs_diff += std::abs(diff);
        }
        break;
    }
    case SCALAR_I: {
        bs.n_compared++;
        if (bs.proto_i != bs.tool_i) {
            if (bs.n_diff == 0) {
                bs.first_diff_proto_i = bs.proto_i;
                bs.first_diff_tool_i  = bs.tool_i;
            }
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
            if (bs.n_diff == 0) {
                bs.first_diff_proto_vf = *bs.proto_vf;
                bs.first_diff_tool_vf  = *bs.tool_vf;
            }
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
            }
        }
        if (event_differs) {
            if (bs.n_diff == 0) {
                bs.first_diff_proto_vf = *bs.proto_vf;
                bs.first_diff_tool_vf  = *bs.tool_vf;
            }
            bs.n_diff++;
        }
        break;
    }
    case VEC_I: {
        if (!bs.proto_vi || !bs.tool_vi) return;
        bs.n_compared++;
        if (bs.proto_vi->size() != bs.tool_vi->size()) {
            if (bs.n_diff == 0) {
                bs.first_diff_proto_vi = *bs.proto_vi;
                bs.first_diff_tool_vi  = *bs.tool_vi;
            }
            bs.n_diff++;
            return;
        }
        for (size_t k = 0; k < bs.proto_vi->size(); ++k) {
            if ((*bs.proto_vi)[k] != (*bs.tool_vi)[k]) {
                if (bs.n_diff == 0) {
                    bs.first_diff_proto_vi = *bs.proto_vi;
                    bs.first_diff_tool_vi  = *bs.tool_vi;
                }
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

        stat_list.push_back(BranchStat{});   // stable address in list
        BranchStat& bs  = stat_list.back();
        bs.name         = bname;
        bs.category     = classify_branch(bname, cats);
        bs.in_kine      = is_kine;
        bs.proto_br     = pb;
        bs.tool_br      = tb;
        setup_branch_stat(bs);

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
//  Auto-locate NeutrinoTaggerInfo.h
//
//  The installed header lives at <prefix>/include/WireCellClus/NeutrinoTaggerInfo.h
//  which is one level up from <prefix>/bin/ where this binary lives.
// ============================================================
static std::string find_tagger_info_header()
{
    char buf[4096] = {0};
    ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n <= 0) return "";
    std::string exe(buf, n);
    size_t slash = exe.rfind('/');
    if (slash == std::string::npos) return "";
    // <prefix>/bin → <prefix>/include/WireCellClus/NeutrinoTaggerInfo.h
    std::string candidate = exe.substr(0, slash) +
                            "/../include/WireCellClus/NeutrinoTaggerInfo.h";
    char resolved[4096];
    if (realpath(candidate.c_str(), resolved) && std::ifstream(resolved).good())
        return std::string(resolved);
    return "";
}

// ============================================================
//  NeutrinoTaggerInfo.h default-value parser
//
//  Reads the header and builds a map: branch_name → default_value_string.
//  Handles:  float/int name{val};   →  name → "val"  (omitted {} → "0")
//            std::vector<…> name;   →  name → "[]"
// ============================================================
static std::map<std::string, std::string>
parse_tagger_defaults(const std::string& path)
{
    std::map<std::string, std::string> defs;
    std::ifstream f(path);
    if (!f) {
        std::cerr << "[warn] cannot open " << path << " — no defaults shown\n";
        return defs;
    }
    std::string line;
    while (std::getline(f, line)) {
        // strip // comments
        auto cpos = line.find("//");
        if (cpos != std::string::npos) line = line.substr(0, cpos);
        // trim leading whitespace
        size_t s = line.find_first_not_of(" \t");
        if (s == std::string::npos) continue;
        line = line.substr(s);

        // vector fields → default []
        if (line.size() >= 12 && line.substr(0, 12) == "std::vector<") {
            size_t close = line.find('>');
            if (close == std::string::npos) continue;
            size_t ns = line.find_first_not_of(" \t", close + 1);
            if (ns == std::string::npos) continue;
            size_t ne = line.find_first_of("{; \t", ns);
            if (ne == std::string::npos) continue;
            defs[line.substr(ns, ne - ns)] = "[]";
            continue;
        }
        // scalar float/int fields
        bool is_float = (line.size() >= 6 && line.substr(0, 6) == "float ");
        bool is_int   = (line.size() >= 4 && line.substr(0, 4) == "int ");
        if (!is_float && !is_int) continue;

        size_t pfx = is_float ? 6 : 4;
        size_t ns  = line.find_first_not_of(" \t", pfx);
        if (ns == std::string::npos) continue;
        size_t brace = line.find('{', ns);
        size_t semi  = line.find(';', ns);
        if (semi == std::string::npos) continue;

        size_t ne = (brace != std::string::npos && brace < semi) ? brace : semi;
        std::string name = line.substr(ns, ne - ns);
        while (!name.empty() && std::isspace((unsigned char)name.back())) name.pop_back();
        if (name.empty()) continue;

        std::string defval = "0";
        if (brace != std::string::npos && brace < semi) {
            size_t cbrace = line.find('}', brace);
            if (cbrace != std::string::npos)
                defval = line.substr(brace + 1, cbrace - brace - 1);
        }
        defs[name] = defval;
    }
    return defs;
}

// ============================================================
//  Sentinel / uninitialized value detection
//
//  Returns a short annotation when one side of a diff looks like a known
//  sentinel and the other does not.  Only scalar types are checked; vector
//  diffs are already obvious from the size and content display.
//
//  Known sentinels:
//    |v| >= 1e7  — "no segment found" style (e.g. min_dis = 1e9 / units::cm = 1e8)
//    v == -1     — uninitialized counter idiom (int branches)
// ============================================================
static std::string sentinel_note(const BranchStat& bs)
{
    auto check_f = [](double a, double b) -> std::string {
        bool a_big = std::abs(a) >= 1e7;
        bool b_big = std::abs(b) >= 1e7;
        if (a_big && !b_big) return "[SENT proto]";
        if (!a_big && b_big) return "[SENT tool]";
        return "";
    };
    auto check_i = [](int a, int b) -> std::string {
        if (a == -1 && b != -1) return "[UNINIT proto]";
        if (b == -1 && a != -1) return "[UNINIT tool]";
        return "";
    };
    switch (bs.type) {
    case SCALAR_F: return check_f(bs.first_diff_proto_f, bs.first_diff_tool_f);
    case SCALAR_D: return check_f(bs.first_diff_proto_d, bs.first_diff_tool_d);
    case SCALAR_I: return check_i(bs.first_diff_proto_i, bs.first_diff_tool_i);
    default:       return "";
    }
}

// ============================================================
//  main
// ============================================================
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: wire-cell-uboone-tagger-compare"
                     " -p <proto.root> -t <toolkit.root>"
                     " [-v] [-n <NeutrinoTaggerInfo.h>]\n";
        return 1;
    }

    TString proto_filename;
    TString toolkit_filename;
    std::string info_header_path;
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-' || argv[i][1] == '\0') continue;
        char flag = argv[i][1];
        // Support both -p<file> and -p <file> forms
        const char* val = (argv[i][2] != '\0') ? argv[i] + 2
                          : (i + 1 < argc)     ? argv[++i]
                                               : nullptr;
        switch (flag) {
        case 'p': if (val) proto_filename   = val; break;
        case 't': if (val) toolkit_filename = val; break;
        case 'n': if (val) info_header_path = val; break;
        case 'v': verbose = true; break;
        }
    }

    if (info_header_path.empty())
        info_header_path = find_tagger_info_header();

    std::map<std::string, std::string> branch_defaults;
    if (!info_header_path.empty())
        branch_defaults = parse_tagger_defaults(info_header_path);

    if (proto_filename.IsNull() || toolkit_filename.IsNull()) {
        std::cerr << "Need -p <proto.root> and -t <toolkit.root>\n";
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

    std::list<BranchStat> all_stats;
    build_stats(proto_tagger, tool_tagger, false, cats, all_stats);
    build_stats(proto_kine,   tool_kine,   true,  cats, all_stats);

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

    // ---- Terminal summary table ----
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

    if (verbose) {
        // ---- Fingerprint cluster report ----------------------------------------
        // Many branches report the same (proto, toolkit) diff value because they
        // all read from the same upstream quantity (main shower energy, main shower
        // length, main vertex, …).  Group scalar diffs by their value fingerprint
        // and surface any cluster of ≥ 3 branches — these point to ONE root cause.
        //
        // Known sentinel patterns flagged by sentinel_note():
        //   |v| >= 1e7  →  "no-segment-found" style (e.g. min_dis = 1e9/units::cm)
        //   int == -1   →  uninitialized counter idiom
        {
            // fingerprint key = "proto_val|tool_val" in scientific(5) notation
            std::map<std::string, std::vector<std::string>> fp_map;
            for (const auto& bs : all_stats) {
                if (bs.n_diff == 0 || bs.type == SKIP) continue;
                std::ostringstream key;
                key << std::scientific << std::setprecision(5);
                switch (bs.type) {
                case SCALAR_F:
                    key << (double)bs.first_diff_proto_f << "|" << (double)bs.first_diff_tool_f;
                    break;
                case SCALAR_D:
                    key << bs.first_diff_proto_d << "|" << bs.first_diff_tool_d;
                    break;
                case SCALAR_I:
                    key.str(""); key << bs.first_diff_proto_i << "|" << bs.first_diff_tool_i;
                    break;
                default: continue;
                }
                fp_map[key.str()].push_back(bs.name);
            }

            // collect clusters of >= 3, sort largest first
            std::vector<std::pair<int, std::string>> clusters;
            for (const auto& kv : fp_map)
                if ((int)kv.second.size() >= 3)
                    clusters.push_back({(int)kv.second.size(), kv.first});
            std::sort(clusters.begin(), clusters.end(),
                      [](const auto& a, const auto& b){ return a.first > b.first; });

            if (!clusters.empty()) {
                std::cout << "\nFingerprint clusters"
                             " (≥3 scalar branches with identical first-diff values"
                             " → shared upstream input):\n"
                          << std::string(95, '-') << "\n";
                constexpr int MAX_SHOW = 6;
                for (const auto& kv : clusters) {
                    int sz = kv.first;
                    const std::string& key = kv.second;
                    size_t pipe = key.find('|');
                    std::cout << "  proto=" << std::setw(14) << key.substr(0, pipe)
                              << "  tool="  << std::setw(14) << key.substr(pipe + 1)
                              << "  (" << sz << " branches)  ";
                    const auto& names = fp_map[key];
                    for (int k = 0; k < std::min(sz, MAX_SHOW); ++k)
                        std::cout << (k ? ", " : "") << names[k];
                    if (sz > MAX_SHOW)
                        std::cout << ", …+" << (sz - MAX_SHOW) << " more";
                    std::cout << "\n";
                }
            }
        }

        // ---- Per-branch verbose table ------------------------------------------
        constexpr int W_NAME   = 46;
        constexpr int W_NDIFF  =  7;  // digits for n_diff and n_compared
        constexpr int W_VAL    = 14;  // width for each numeric value
        constexpr int W_NOTE   = 14;  // width for sentinel/uninit annotation

        auto fmt_val = [&](std::ostream& os, double v) -> std::ostream& {
            os << std::setw(W_VAL) << std::right << std::scientific
               << std::setprecision(5) << v;
            return os;
        };
        auto fmt_ival = [&](std::ostream& os, int v) -> std::ostream& {
            os << std::setw(W_VAL) << std::right << v;
            return os;
        };

        std::cout << "\n"
                  << std::left  << std::setw(W_NAME) << "branch"
                  << "  " << std::right << std::setw(W_NDIFF) << "n_diff"
                  << "/" << std::left  << std::setw(W_NDIFF) << "n_cmp"
                  << "  " << std::right << std::setw(W_VAL) << "max|diff|"
                  << "  " << std::right << std::setw(W_VAL) << "proto"
                  << "  " << std::right << std::setw(W_VAL) << "toolkit"
                  << "  " << std::right << std::setw(W_VAL) << "default"
                  << "  " << std::left  << std::setw(W_NOTE) << "note"
                  << "\n"
                  << std::string(W_NAME + 2 + 2*W_NDIFF + 1 + 4*(W_VAL+2) + 2 + W_NOTE, '-') << "\n";

        std::cout << std::defaultfloat;

        // Branch name → struct field name aliases (output visitor uses different names).
        static const std::map<std::string,std::string> kAliases = {
            {"nu_x", "kine_nu_x_corr"},
            {"nu_y", "kine_nu_y_corr"},
            {"nu_z", "kine_nu_z_corr"},
        };
        auto lookup = [&](const std::string& name) -> std::string {
            auto it = branch_defaults.find(name);
            if (it != branch_defaults.end()) return it->second;
            auto al = kAliases.find(name);
            if (al != kAliases.end()) {
                auto it2 = branch_defaults.find(al->second);
                if (it2 != branch_defaults.end()) return it2->second;
            }
            return "";
        };

        for (const auto& bs : all_stats) {
            if (bs.n_diff == 0 || bs.type == SKIP) continue;

            std::string defstr  = lookup(bs.name);
            std::string note    = sentinel_note(bs);

            switch (bs.type) {
            case SCALAR_F:
            case SCALAR_D:
            case SCALAR_I: {
                std::cout << std::left  << std::setw(W_NAME) << bs.name
                          << "  " << std::right << std::setw(W_NDIFF) << bs.n_diff
                          << "/" << std::left  << std::setw(W_NDIFF) << bs.n_compared
                          << "  ";
                fmt_val(std::cout, bs.max_abs_diff) << "  ";
                if (bs.type == SCALAR_F) {
                    fmt_val(std::cout, bs.first_diff_proto_f) << "  ";
                    fmt_val(std::cout, bs.first_diff_tool_f)  << "  ";
                } else if (bs.type == SCALAR_D) {
                    fmt_val(std::cout, bs.first_diff_proto_d) << "  ";
                    fmt_val(std::cout, bs.first_diff_tool_d)  << "  ";
                } else {
                    fmt_ival(std::cout, bs.first_diff_proto_i) << "  ";
                    fmt_ival(std::cout, bs.first_diff_tool_i)  << "  ";
                }
                std::cout << std::right << std::setw(W_VAL)
                          << (defstr.empty() ? "-" : defstr)
                          << "  " << std::left << note << "\n";
                break;
            }
            case VEC_F:
            case VEC_I: {
                // Header line with stats
                std::cout << std::left  << std::setw(W_NAME) << bs.name
                          << "  " << std::right << std::setw(W_NDIFF) << bs.n_diff
                          << "/" << std::left  << std::setw(W_NDIFF) << bs.n_compared
                          << "  ";
                fmt_val(std::cout, bs.max_abs_diff);
                if (!defstr.empty())
                    std::cout << "  default=" << defstr;
                std::cout << "\n";
                // proto / toolkit rows
                if (bs.type == VEC_F) {
                    std::cout << "    proto  : [";
                    for (size_t k = 0; k < bs.first_diff_proto_vf.size(); ++k)
                        std::cout << (k ? ", " : "") << bs.first_diff_proto_vf[k];
                    std::cout << "]\n    toolkit: [";
                    for (size_t k = 0; k < bs.first_diff_tool_vf.size(); ++k)
                        std::cout << (k ? ", " : "") << bs.first_diff_tool_vf[k];
                    std::cout << "]\n";
                } else {
                    std::cout << "    proto  : [";
                    for (size_t k = 0; k < bs.first_diff_proto_vi.size(); ++k)
                        std::cout << (k ? ", " : "") << bs.first_diff_proto_vi[k];
                    std::cout << "]\n    toolkit: [";
                    for (size_t k = 0; k < bs.first_diff_tool_vi.size(); ++k)
                        std::cout << (k ? ", " : "") << bs.first_diff_tool_vi[k];
                    std::cout << "]\n";
                }
                break;
            }
            default: break;
            }
        }
    }

    proto_file->Close();
    tool_file->Close();

    return any_diff ? 1 : 0;
}
