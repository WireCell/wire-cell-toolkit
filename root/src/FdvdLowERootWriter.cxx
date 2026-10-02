#include "WireCellRoot/FdvdLowERootWriter.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <limits>
#include <map>

WIRECELL_FACTORY(FdvdLowERootWriter, WireCell::Root::FdvdLowERootWriter,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Root;

namespace {
    constexpr int NCH = 184;
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    // a named f8 tensor as (rows, ncol); missing -> 0 rows
    struct Tab {
        const double* d{nullptr};
        size_t n{0}, nc{0};
        double at(size_t r, size_t c) const { return d[r * nc + c]; }
    };
    Tab table(const ITensorSet::pointer& in, const std::string& name)
    {
        for (const auto& t : *in->tensors()) {
            if (t->metadata()["name"].asString() != name) continue;
            const auto sh = t->shape();
            if (sh.size() != 2 || t->dtype().find("f8") == std::string::npos) {
                raise<ValueError>("FdvdLowERootWriter: tensor %s is not a 2-D f8 table", name);
            }
            return Tab{(const double*) t->data(), sh[0], sh[1]};
        }
        return Tab{};
    }
}  // namespace

struct FdvdLowERootWriter::Buf {
    int event{0};
    // T_meta
    std::string arms, specs, producer;
    // T_event
    int n_tree_clusters{0}, n_flash{0}, n_groups_stored{0}, n_window_groups{0}, n_clusters{0}, n_pairs{0}, n_drift{0};
    Long64_t n_tree_blobs{0}, n_unjoined{0}, n_table_blobs{0};
    // T_cluster
    int k{0}, tree_index{0}, rep{0}, nblob{0};
    double Q{0}, x_app{0}, t_lo{0}, t_hi{0}, mu{0}, sigma{0};
    std::vector<int> match_g;
    std::vector<double> match_t0, match_x;
    // T_flash
    int g{0}, group_index{0}, npd{0};
    double t{0}, tot{0};
    float pe[NCH];
    // T_pair
    float r{0}, ks{0}, dc{0}, npp{0};
    std::vector<int> pass;
    // T_match
    int arm{0}, spec{0};
    double t0{0}, x_corr{0};
};

FdvdLowERootWriter::FdvdLowERootWriter()
  : Aux::Logger("FdvdLowERootWriter", "root")
{
}

FdvdLowERootWriter::~FdvdLowERootWriter() { close(); }

WireCell::Configuration FdvdLowERootWriter::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = m_filename;
    cfg["vcm"] = m_vcm;
    cfg["event"] = m_event;
    cfg["pairs"] = m_pairs;
    cfg["pair_qmin"] = m_pair_qmin;
    cfg["pair_npd"] = m_pair_npd;
    return cfg;
}

void FdvdLowERootWriter::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "output_filename", m_filename);
    m_vcm = get(cfg, "vcm", m_vcm);
    m_event = get(cfg, "event", m_event);
    m_pairs = get(cfg, "pairs", m_pairs);
    m_pair_qmin = get(cfg, "pair_qmin", m_pair_qmin);
    m_pair_npd = get(cfg, "pair_npd", m_pair_npd);
    if (m_pairs != "qc" && m_pairs != "all" && m_pairs != "none") {
        raise<ValueError>("FdvdLowERootWriter: pairs must be qc, all or none, got '%s'", m_pairs);
    }
}

void FdvdLowERootWriter::open(const ITensorSet::pointer& in)
{
    const auto md = in->metadata();
    m_narm = md["arms"].size();
    m_nspec = md["specs"].size();
    m_buf = std::make_unique<Buf>();
    auto& b = *m_buf;
    for (Json::ArrayIndex i = 0; i < md["arms"].size(); ++i) b.arms += (i ? "," : "") + md["arms"][i].asString();
    for (Json::ArrayIndex i = 0; i < md["specs"].size(); ++i) b.specs += (i ? "," : "") + md["specs"][i].asString();
    b.producer = md["producer"].asString();
    m_file = TFile::Open(m_filename.c_str(), "RECREATE");
    if (!m_file || m_file->IsZombie()) raise<IOError>("FdvdLowERootWriter: cannot open %s", m_filename);

    m_tmeta = new TTree("T_meta", "FD-VD low-energy reconstruction: arms, specs");
    m_tmeta->Branch("arms", &b.arms);
    m_tmeta->Branch("specs", &b.specs);
    m_tmeta->Branch("producer", &b.producer);
    m_tmeta->Branch("vcm", &m_vcm);
    m_tmeta->Fill();

    m_tevent = new TTree("T_event", "per event counts");
    m_tevent->Branch("event", &b.event);
    m_tevent->Branch("n_tree_clusters", &b.n_tree_clusters);
    m_tevent->Branch("n_tree_blobs", &b.n_tree_blobs);
    m_tevent->Branch("n_unjoined", &b.n_unjoined);
    m_tevent->Branch("n_table_blobs", &b.n_table_blobs);
    m_tevent->Branch("n_flash", &b.n_flash);
    m_tevent->Branch("n_groups_stored", &b.n_groups_stored);
    m_tevent->Branch("n_window_groups", &b.n_window_groups);
    m_tevent->Branch("n_clusters", &b.n_clusters);
    m_tevent->Branch("n_pairs", &b.n_pairs);
    m_tevent->Branch("n_drift", &b.n_drift);

    m_tcluster = new TTree("T_cluster", "Q-L clusters");
    m_tcluster->Branch("event", &b.event);
    m_tcluster->Branch("k", &b.k);
    m_tcluster->Branch("tree_index", &b.tree_index);
    m_tcluster->Branch("rep", &b.rep);
    m_tcluster->Branch("nblob", &b.nblob);
    m_tcluster->Branch("Q", &b.Q);
    m_tcluster->Branch("x_app", &b.x_app);
    m_tcluster->Branch("t_lo", &b.t_lo);
    m_tcluster->Branch("t_hi", &b.t_hi);
    m_tcluster->Branch("mu", &b.mu);
    m_tcluster->Branch("sigma", &b.sigma);
    m_tcluster->Branch("match_g", &b.match_g);
    m_tcluster->Branch("match_t0", &b.match_t0);
    m_tcluster->Branch("match_x", &b.match_x);

    m_tflash = new TTree("T_flash", "flash groups in the time window");
    m_tflash->Branch("event", &b.event);
    m_tflash->Branch("g", &b.g);
    m_tflash->Branch("group_index", &b.group_index);
    m_tflash->Branch("t", &b.t);
    m_tflash->Branch("tot", &b.tot);
    m_tflash->Branch("npd", &b.npd);
    m_tflash->Branch("pe", b.pe, "pe[184]/F");

    m_tpair = new TTree("T_pair", "cluster x flash-group pairs in the drift window");
    m_tpair->Branch("event", &b.event);
    m_tpair->Branch("k", &b.k);
    m_tpair->Branch("g", &b.g);
    m_tpair->Branch("r", &b.r);
    m_tpair->Branch("ks", &b.ks);
    m_tpair->Branch("dc", &b.dc);
    m_tpair->Branch("npp", &b.npp);
    m_tpair->Branch("pass", &b.pass);

    m_tmatch = new TTree("T_match", "decisions: one row per matched cluster, arm and spec");
    m_tmatch->Branch("event", &b.event);
    m_tmatch->Branch("arm", &b.arm);
    m_tmatch->Branch("spec", &b.spec);
    m_tmatch->Branch("k", &b.k);
    m_tmatch->Branch("g", &b.g);
    m_tmatch->Branch("tree_index", &b.tree_index);
    m_tmatch->Branch("rep", &b.rep);
    m_tmatch->Branch("t0", &b.t0);
    m_tmatch->Branch("x_corr", &b.x_corr);

    m_tdrift = new TTree("T_drift", "drift regressor per regressed cluster");
    m_tdrift->Branch("event", &b.event);
    m_tdrift->Branch("rep", &b.rep);
    m_tdrift->Branch("mu", &b.mu);
    m_tdrift->Branch("sigma", &b.sigma);
}

void FdvdLowERootWriter::close()
{
    if (!m_file) return;
    m_file->cd();
    m_file->Write();
    m_file->Close();
    delete m_file;   // owns the trees
    m_file = nullptr;
}

bool FdvdLowERootWriter::operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out)
{
    out = in;   // pass-through: the set continues downstream unchanged
    if (!in) {
        log->debug("EOS at call={}, closing {}", m_count, m_filename);
        close();
        return true;
    }
    if (!m_file) open(in);
    auto& b = *m_buf;
    const auto md = in->metadata();
    b.event = m_event >= 0 ? m_event : in->ident();
    const auto C = table(in, "clusters"), G = table(in, "groups"), P = table(in, "pairs"), D = table(in, "decisions"),
               R = table(in, "drift");
    if (md["arms"].size() != m_narm || md["specs"].size() != m_nspec) {
        raise<ValueError>("FdvdLowERootWriter: arms/specs changed between events");
    }

    b.n_tree_clusters = md["n_tree_clusters"].asInt();
    b.n_tree_blobs = md["n_tree_blobs"].asInt64();
    b.n_unjoined = md["n_tree_blobs_unjoined"].asInt64();
    b.n_table_blobs = md["n_table_blobs"].asInt64();
    b.n_flash = md["n_flash"].asInt();
    b.n_groups_stored = md["n_groups_stored"].asInt();
    b.n_window_groups = (int) G.n;
    b.n_clusters = (int) C.n;
    b.n_pairs = (int) P.n;
    b.n_drift = (int) R.n;
    m_tevent->Fill();

    // decisions per cluster: index arm * nspec + spec
    const size_t nas = m_narm * m_nspec;
    std::map<int, std::vector<int>> mg;
    for (size_t r = 0; r < D.n; ++r) {
        auto& v = mg[(int) D.at(r, 2)];
        if (v.empty()) v.assign(nas, -1);
        v[(size_t) D.at(r, 0) * m_nspec + (size_t) D.at(r, 1)] = (int) D.at(r, 3);
    }
    auto gtime = [&](int g) { return (g >= 0 && (size_t) g < G.n) ? G.at(g, 2) : NaN; };
    for (size_t k = 0; k < C.n; ++k) {
        b.k = (int) k;
        b.tree_index = (int) C.at(k, 0);
        b.rep = (int) C.at(k, 1);
        b.nblob = (int) C.at(k, 2);
        b.Q = C.at(k, 3);
        b.x_app = C.at(k, 4);
        b.t_lo = C.at(k, 5);
        b.t_hi = C.at(k, 6);
        b.mu = C.at(k, 8);
        b.sigma = C.at(k, 9);
        const auto it = mg.find((int) k);
        b.match_g = it == mg.end() ? std::vector<int>(nas, -1) : it->second;
        b.match_t0.assign(nas, NaN);
        b.match_x.assign(nas, NaN);
        for (size_t a = 0; a < nas; ++a) {
            if (b.match_g[a] < 0) continue;
            b.match_t0[a] = gtime(b.match_g[a]);
            b.match_x[a] = b.x_app + m_vcm * b.match_t0[a];
        }
        m_tcluster->Fill();
    }
    for (size_t g = 0; g < G.n; ++g) {
        b.g = (int) G.at(g, 0);
        b.group_index = (int) G.at(g, 1);
        b.t = G.at(g, 2);
        b.tot = G.at(g, 3);
        b.npd = (int) G.at(g, 4);
        for (int ch = 0; ch < NCH && 5 + (size_t) ch < G.nc; ++ch) b.pe[ch] = (float) G.at(g, 5 + ch);
        m_tflash->Fill();
    }
    for (size_t p = 0; p < P.n && m_pairs != "none"; ++p) {
        b.k = (int) P.at(p, 0);
        b.g = (int) P.at(p, 1);
        if (m_pairs == "qc" && !(C.at(b.k, 3) >= m_pair_qmin && G.at(b.g, 4) >= m_pair_npd)) continue;
        b.r = (float) P.at(p, 2);
        b.ks = (float) P.at(p, 3);
        b.dc = (float) P.at(p, 4);
        b.npp = (float) P.at(p, 5);
        b.pass.assign(m_nspec, 0);
        for (size_t s = 0; s < m_nspec && 6 + s < P.nc; ++s) b.pass[s] = (int) P.at(p, 6 + s);
        m_tpair->Fill();
    }
    for (size_t r = 0; r < D.n; ++r) {
        b.arm = (int) D.at(r, 0);
        b.spec = (int) D.at(r, 1);
        b.k = (int) D.at(r, 2);
        b.g = (int) D.at(r, 3);
        b.tree_index = (int) C.at(b.k, 0);
        b.rep = (int) C.at(b.k, 1);
        b.t0 = gtime(b.g);
        b.x_corr = C.at(b.k, 4) + m_vcm * b.t0;
        m_tmatch->Fill();
    }
    for (size_t r = 0; r < R.n; ++r) {
        b.rep = (int) R.at(r, 0);
        b.mu = R.at(r, 1);
        b.sigma = R.at(r, 2);
        m_tdrift->Fill();
    }
    log->debug("call={} event={} clusters {} flashes {} pairs {} matches {} drift {}", m_count, b.event, C.n, G.n, P.n,
               D.n, R.n);
    ++m_count;
    return true;
}
