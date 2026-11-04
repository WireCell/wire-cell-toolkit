// ai-helper issue 39: see the header for the layout.
#include "WireCellClus/PerEventName.h"
#include "WireCellRoot/RootFlatTreeVisitor.h"

#include "WireCellClus/Facade_Ensemble.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include "TSystem.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TFile.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TTree.h"

#include <memory>
#include <set>
#include <stdexcept>

WIRECELL_FACTORY(RootFlatTreeVisitor, WireCell::Root::RootFlatTreeVisitor,
                 WireCell::Clus::IEnsembleVisitor, WireCell::IConfigurable)

using namespace WireCell;

namespace {

    const int kBufsize = 4096;

    // ROOT leaf-list type code of a C++ type.
    template <typename T> char code();
    template <> char code<Int_t>() { return 'I'; }
    template <> char code<UInt_t>() { return 'i'; }
    template <> char code<Float_t>() { return 'F'; }
    template <> char code<Double_t>() { return 'D'; }
    template <> char code<Long64_t>() { return 'L'; }
    template <> char code<ULong64_t>() { return 'l'; }
    template <> char code<Short_t>() { return 'S'; }
    template <> char code<UShort_t>() { return 's'; }
    template <> char code<Char_t>() { return 'B'; }
    template <> char code<UChar_t>() { return 'b'; }

    // An Int_t column of the output (a counter, a ..length or ..idx array).
    struct IntArray {
        std::vector<Int_t> v;
    };

    struct Column {
        virtual ~Column() = default;
        virtual void bind(TTree* in) = 0;    // set the input branch address
        virtual void take() = 0;             // after in->GetEntry(row)
        // create the output branches (counters first) under `base` with row
        // counter `rows`; called after all rows are taken.
        virtual void define(TTree* out, const std::string& base, const std::string& rows) = 0;
    };

    // Make a branch "name[count]/X" on `out` whose data is `data` (n >= 0
    // elements).  An empty vector still needs a valid address.
    template <typename T>
    void array_branch(TTree* out, const std::string& name, const std::string& count, std::vector<T>& data)
    {
        if (data.empty()) data.reserve(1);
        const std::string leaf = name + "[" + count + "]/" + code<T>();
        out->Branch(name.c_str(), (void*) data.data(), leaf.c_str(), kBufsize);
    }
    void scalar_int_branch(TTree* out, const std::string& name, Int_t* addr)
    {
        out->Branch(name.c_str(), (void*) addr, (name + "/I").c_str(), kBufsize);
    }

    template <typename T>
    struct ScalarColumn : Column {
        std::string name;
        T in{};
        std::vector<T> out;
        explicit ScalarColumn(const std::string& n) : name(n) {}
        void bind(TTree* t) override { t->SetBranchAddress(name.c_str(), &in); }
        void take() override { out.push_back(in); }
        void define(TTree* o, const std::string& base, const std::string& rows) override
        {
            array_branch(o, base + "." + name, rows, out);
        }
    };

    // vector<T>, or std::string as a vector of Char_t.
    template <typename T, typename Container = std::vector<T>>
    struct VectorColumn : Column {
        std::string name;
        Container* in{nullptr};
        std::vector<Int_t> len, idx;
        Int_t tot{0};
        std::vector<T> flat;
        explicit VectorColumn(const std::string& n) : name(n) {}
        void bind(TTree* t) override { t->SetBranchAddress(name.c_str(), &in); }
        void take() override
        {
            const size_t n = in ? in->size() : 0;
            idx.push_back((Int_t) flat.size());
            len.push_back((Int_t) n);
            if (n) flat.insert(flat.end(), in->begin(), in->end());
            tot = (Int_t) flat.size();
        }
        void define(TTree* o, const std::string& base, const std::string& rows) override
        {
            const std::string b = base + "." + name;
            scalar_int_branch(o, b + "..totarraysize", &tot);
            array_branch(o, b + "..length", rows, len);
            array_branch(o, b + "..idx", rows, idx);
            array_branch(o, b, b + "..totarraysize", flat);
        }
    };

    template <typename T>
    struct VectorVectorColumn : Column {
        std::string name;
        std::vector<std::vector<T>>* in{nullptr};
        std::vector<Int_t> len, idx;         // per row: outer elements
        Int_t tot{0};                        // outer elements in all rows
        std::vector<Int_t> vlen, vidx;       // per outer element: inner elements
        Int_t vtot{0};
        std::vector<T> flat;
        explicit VectorVectorColumn(const std::string& n) : name(n) {}
        void bind(TTree* t) override { t->SetBranchAddress(name.c_str(), &in); }
        void take() override
        {
            const size_t n = in ? in->size() : 0;
            idx.push_back(tot);
            len.push_back((Int_t) n);
            for (size_t k = 0; k < n; ++k) {
                const auto& inner = (*in)[k];
                vidx.push_back((Int_t) flat.size());
                vlen.push_back((Int_t) inner.size());
                flat.insert(flat.end(), inner.begin(), inner.end());
            }
            tot += (Int_t) n;
            vtot = (Int_t) flat.size();
        }
        void define(TTree* o, const std::string& base, const std::string& rows) override
        {
            const std::string b = base + "." + name;
            scalar_int_branch(o, b + "..totarraysize", &tot);
            array_branch(o, b + "..length", rows, len);
            array_branch(o, b + "..idx", rows, idx);
            scalar_int_branch(o, b + ".v..totarraysize", &vtot);
            array_branch(o, b + ".v..length", b + "..totarraysize", vlen);
            array_branch(o, b + ".v..idx", b + "..totarraysize", vidx);
            array_branch(o, b + ".v", b + ".v..totarraysize", flat);
        }
    };

    // Bool_t: std::vector<bool> has no data(), so keep the bytes and declare the
    // leaf as Bool_t ('O', one byte each).
    struct BoolColumn : Column {
        std::string name;
        Bool_t in{};
        std::vector<UChar_t> out;
        explicit BoolColumn(const std::string& n) : name(n) {}
        void bind(TTree* t) override { t->SetBranchAddress(name.c_str(), &in); }
        void take() override { out.push_back(in ? 1 : 0); }
        void define(TTree* o, const std::string& base, const std::string& rows) override
        {
            const std::string b = base + "." + name;
            if (out.empty()) out.reserve(1);
            o->Branch(b.c_str(), (void*) out.data(), (b + "[" + rows + "]/O").c_str(), kBufsize);
        }
    };

    template <typename T>
    std::unique_ptr<Column> make_scalar(const std::string& n) { return std::make_unique<ScalarColumn<T>>(n); }

    std::unique_ptr<Column> make_column(TBranch* br, std::string& why)
    {
        const std::string name = br->GetName();
        const std::string cls = br->GetClassName();
        if (!cls.empty()) {
            if (cls == "string") return std::make_unique<VectorColumn<Char_t, std::string>>(name);
            if (cls == "vector<int>") return std::make_unique<VectorColumn<Int_t>>(name);
            if (cls == "vector<unsigned int>") return std::make_unique<VectorColumn<UInt_t>>(name);
            if (cls == "vector<float>") return std::make_unique<VectorColumn<Float_t>>(name);
            if (cls == "vector<double>") return std::make_unique<VectorColumn<Double_t>>(name);
            if (cls == "vector<short>") return std::make_unique<VectorColumn<Short_t>>(name);
            if (cls == "vector<Long64_t>" || cls == "vector<long long>") return std::make_unique<VectorColumn<Long64_t>>(name);
            if (cls == "vector<vector<int> >") return std::make_unique<VectorVectorColumn<Int_t>>(name);
            if (cls == "vector<vector<float> >") return std::make_unique<VectorVectorColumn<Float_t>>(name);
            if (cls == "vector<vector<double> >") return std::make_unique<VectorVectorColumn<Double_t>>(name);
            why = "class " + cls;
            return nullptr;
        }
        auto* leaves = br->GetListOfLeaves();
        if (!leaves || leaves->GetEntries() != 1) {
            why = "branch with " + std::to_string(leaves ? leaves->GetEntries() : 0) + " leaves";
            return nullptr;
        }
        auto* leaf = (TLeaf*) leaves->At(0);
        if (leaf->GetLeafCount() || leaf->GetLen() != 1) {
            why = "array leaf " + std::string(leaf->GetTitle());
            return nullptr;
        }
        const std::string t = leaf->GetTypeName();
        if (t == "Int_t") return make_scalar<Int_t>(name);
        if (t == "UInt_t") return make_scalar<UInt_t>(name);
        if (t == "Float_t") return make_scalar<Float_t>(name);
        if (t == "Double_t") return make_scalar<Double_t>(name);
        if (t == "Long64_t") return make_scalar<Long64_t>(name);
        if (t == "ULong64_t") return make_scalar<ULong64_t>(name);
        if (t == "Short_t") return make_scalar<Short_t>(name);
        if (t == "UShort_t") return make_scalar<UShort_t>(name);
        if (t == "Char_t") return make_scalar<Char_t>(name);
        if (t == "UChar_t") return make_scalar<UChar_t>(name);
        if (t == "Bool_t") return std::make_unique<BoolColumn>(name);
        why = "leaf type " + t;
        return nullptr;
    }

}  // namespace

Root::RootFlatTreeVisitor::RootFlatTreeVisitor()
  : log(Log::logger("root"))
{
}

Root::RootFlatTreeVisitor::~RootFlatTreeVisitor() {}

void Root::RootFlatTreeVisitor::configure(const WireCell::Configuration& cfg)
{
    m_output_filename = get<std::string>(cfg, "output_filename", m_output_filename);
    m_tree_name = get<std::string>(cfg, "tree_name", m_tree_name);
    m_prefix = get<std::string>(cfg, "prefix", m_prefix);
    m_drop_row_trees = get<bool>(cfg, "drop_row_trees", m_drop_row_trees);
}

WireCell::Configuration Root::RootFlatTreeVisitor::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = m_output_filename;  // the SAME (per-event) name the tracking writers use
    cfg["tree_name"] = m_tree_name;
    cfg["prefix"] = m_prefix;
    cfg["drop_row_trees"] = m_drop_row_trees;    // true: remove the row trees after flattening
    return cfg;
}

bool Root::RootFlatTreeVisitor::flatten_file(const std::string& filename, const std::string& tree_name,
                                             const std::string& prefix, bool drop_row_trees,
                                             std::string& error)
{
    // "UPDATE" CREATES a missing file, which then got one entry with run = -1
    // and valid = 1 (PR 544 review, item 2).  The file must already be there.
    if (gSystem->AccessPathName(filename.c_str())) {   // kTRUE = does NOT exist
        error = "no such file " + filename;
        return false;
    }
    std::unique_ptr<TFile> tf(TFile::Open(filename.c_str(), "UPDATE"));
    if (!tf || tf->IsZombie()) {
        error = "cannot open " + filename;
        return false;
    }

    // The trees of the file, in key order, once each (highest cycle).
    std::vector<std::string> names;
    std::set<std::string> seen;
    for (auto* obj : *tf->GetListOfKeys()) {
        auto* key = (TKey*) obj;
        const std::string cls = key->GetClassName();
        const std::string nm = key->GetName();
        if (cls != "TTree" || nm == tree_name || seen.count(nm)) continue;
        seen.insert(nm);
        names.push_back(nm);
    }

    struct TreeCols {
        std::string name;
        Int_t nrows{0};
        std::vector<std::unique_ptr<Column>> cols;
    };
    std::vector<TreeCols> all;
    Int_t run = -1, subrun = -1, evt = -1;
    bool have_trun = false;   // valid = 1 only with the header RSE (PR 544 review, item 2)

    for (const auto& nm : names) {
        auto* in = dynamic_cast<TTree*>(tf->Get(nm.c_str()));
        if (!in) {
            error = "cannot read tree " + nm;
            return false;
        }
        TreeCols tc;
        tc.name = nm;
        tc.nrows = (Int_t) in->GetEntries();
        for (auto* bobj : *in->GetListOfBranches()) {
            auto* br = (TBranch*) bobj;
            std::string why;
            auto col = make_column(br, why);
            if (!col) {
                error = "unsupported branch " + nm + "." + br->GetName() + " (" + why + ")";
                return false;
            }
            col->bind(in);
            tc.cols.push_back(std::move(col));
        }
        for (Long64_t i = 0; i < in->GetEntries(); ++i) {
            in->GetEntry(i);
            for (auto& c : tc.cols) c->take();
        }
        if (nm == "Trun" && tc.nrows > 0) {
            // header RSE, exactly as written in Trun
            Int_t r = -1, s = -1, e = -1;
            in->ResetBranchAddresses();
            if (in->GetBranch("runNo")) in->SetBranchAddress("runNo", &r);
            if (in->GetBranch("subRunNo")) in->SetBranchAddress("subRunNo", &s);
            if (in->GetBranch("eventNo")) in->SetBranchAddress("eventNo", &e);
            in->GetEntry(0);
            run = r; subrun = s; evt = e;
            have_trun = true;
        }
        in->ResetBranchAddresses();
        all.push_back(std::move(tc));
    }

    tf->cd();
    auto* out = new TTree(tree_name.c_str(), "Wire-Cell per-event tree (ai-helper issue 39)");
    out->SetDirectory(tf.get());
    Int_t valid = have_trun ? 1 : 0;
    scalar_int_branch(out, prefix + ".run", &run);
    scalar_int_branch(out, prefix + ".subrun", &subrun);
    scalar_int_branch(out, prefix + ".evt", &evt);
    scalar_int_branch(out, prefix + ".valid", &valid);
    for (auto& tc : all) {
        const std::string base = prefix + "." + tc.name;
        const std::string rows = base + "..length";
        scalar_int_branch(out, rows, &tc.nrows);
        for (auto& c : tc.cols) c->define(out, base, rows);
    }
    out->Fill();
    out->Write("", TObject::kOverwrite);

    if (drop_row_trees) {
        for (const auto& nm : names) tf->Delete((nm + ";*").c_str());
    }
    tf->Close();
    return true;
}

void Root::RootFlatTreeVisitor::visit(Clus::Facade::Ensemble& ensemble) const
{
    const std::string fname = Clus::per_event_name(m_output_filename, ensemble);
    std::string error;
    if (!flatten_file(fname, m_tree_name, m_prefix, m_drop_row_trees, error)) {
        log->error("RootFlatTreeVisitor: {}: no {} written: {}", fname, m_tree_name, error);
        return;
    }
    log->debug("RootFlatTreeVisitor: {}: wrote {}", fname, m_tree_name);
}
