#include "WireCellMatch/FdvdLowEQLMatching.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellClus/Facade.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/String.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

WIRECELL_FACTORY(FdvdLowEQLMatching, WireCell::Match::FdvdLowEQLMatching,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Match;
namespace L = WireCell::Match::FdvdLowE;

FdvdLowEQLMatching::FdvdLowEQLMatching()
  : Aux::Logger("FdvdLowEQLMatching", "match")
{
}

FdvdLowEQLMatching::~FdvdLowEQLMatching() {}

std::vector<std::string> FdvdLowEQLMatching::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

WireCell::Configuration FdvdLowEQLMatching::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = (int) m_multiplicity;
    cfg["inpath"] = m_inpath;
    cfg["library"] = m_library;
    cfg["calibration"] = m_calibration;
    cfg["geom_file"] = m_geom_file;
    cfg["dump_tables"] = m_dump_tables;
    cfg["qmin_cluster"] = m_C.qmin_cluster;
    for (size_t i = 0; i < m_veto_k.size(); ++i) {
        cfg["specs"][(int) i]["name"] = m_spec_names[i];
        cfg["specs"][(int) i]["veto_k"] = m_veto_k[i];
    }
    return cfg;
}

void FdvdLowEQLMatching::configure(const WireCell::Configuration& cfg)
{
    const int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m != 3 && m != 4) raise<ValueError>("FdvdLowEQLMatching: multiplicity must be 3 or 4 (with drift)");
    m_multiplicity = m;
    m_inpath = get(cfg, "inpath", m_inpath);
    m_library = get(cfg, "library", m_library);
    m_calibration = get(cfg, "calibration", m_calibration);
    m_geom_file = get(cfg, "geom_file", m_geom_file);
    m_dump_tables = get(cfg, "dump_tables", m_dump_tables);
    m_C.qmin_cluster = get(cfg, "qmin_cluster", m_C.qmin_cluster);
    if (cfg["specs"].isArray()) {
        m_veto_k.clear();
        m_spec_names.clear();
        for (const auto& s : cfg["specs"]) {
            m_spec_names.push_back(s["name"].asString());
            m_veto_k.push_back(s["veto_k"].asDouble());
        }
    }
    if (m_library.empty() || m_calibration.empty() || m_geom_file.empty()) {
        raise<ValueError>("FdvdLowEQLMatching: library, calibration and geom_file are required");
    }
    m_lib = std::make_unique<PhotonLibraryModel>(m_library);
    m_drift.reset();
    if (cfg["drift"].isObject()) {
        if (m_multiplicity == 4) raise<ValueError>("FdvdLowEQLMatching: give a drift block or a drift port, not both");
        m_drift = std::make_unique<FdvdDriftRegressor>(cfg["drift"]);
    }
    if ((int) m_lib->nchan() != L::NCH) raise<ValueError>("FdvdLowEQLMatching: library has %d channels", (int) m_lib->nchan());
    m_cal = L::load_calibration(m_calibration);
    m_arms = L::default_arms();
    // ql_m2m_proto.py PD_POS (line 61): OpDet order, mm / 10
    const auto geo = Persist::load(Persist::resolve(m_geom_file));
    std::map<int, std::array<double, 3>> bypd;
    for (const auto& o : geo["opdets"]) {
        bypd[o["opdet"].asInt()] = {o["x"].asDouble() / 10.0, o["y"].asDouble() / 10.0, o["z"].asDouble() / 10.0};
    }
    m_pd_pos.clear();
    for (const auto& [od, p] : bypd) m_pd_pos.push_back(p);
    if ((int) m_pd_pos.size() != L::NCH) raise<ValueError>("FdvdLowEQLMatching: geometry has %d OpDets", (int) m_pd_pos.size());
}

namespace {
    ITensor::pointer named(const ITensorSet::pointer& ts, const std::string& name)
    {
        for (const auto& t : *ts->tensors()) {
            if (t->metadata()["name"].asString() == name) return t;
        }
        return nullptr;
    }

    template <typename T>
    ITensor::pointer tensor(const std::string& name, const std::vector<T>& data, size_t ncol)
    {
        Configuration md;
        md["name"] = name;
        const size_t nrow = ncol ? data.size() / ncol : 0;
        return std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{nrow, ncol}, data.data(), md);
    }
}  // namespace

bool FdvdLowEQLMatching::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    size_t neos = 0;
    for (const auto& in : invec) neos += !in;
    if (neos == invec.size()) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    if (neos) raise<ValueError>("FdvdLowEQLMatching: %d of %d inputs at EOS", (int) neos, (int) invec.size());

    // port 1: blob table
    const auto bt = named(invec[1], "blobs");
    if (!bt || bt->shape().size() != 2) raise<ValueError>("FdvdLowEQLMatching: no 'blobs' table on port 1");
    const size_t nbt = bt->shape()[0], ncolb = bt->shape()[1];
    const double* B = (const double*) bt->data();

    // port 0: clustering tree
    const int ident = invec[0]->ident();
    std::string inpath = m_inpath;
    if (inpath.find('%') != std::string::npos) inpath = String::format(inpath, ident);
    auto root = Aux::TensorDM::as_pctree(*invec[0]->tensors(), inpath + "/live");
    if (!root) raise<ValueError>("FdvdLowEQLMatching: no point-cloud tree at %s/live", inpath);
    auto grouping = root->value.facade<Clus::Facade::Grouping>();

    // drift: port 3 or the regressor
    std::map<int, std::pair<double, double>> drift;
    FdvdDriftRegressor::Result dres;
    if (m_drift) {
        dres = m_drift->compute(*grouping, (const double*) bt->data(), bt->shape()[0], bt->shape()[1], ident);
        for (size_t r = 0; r + 2 < dres.rows.size(); r += 3) drift[(int) dres.rows[r]] = {dres.rows[r + 1], dres.rows[r + 2]};
    }
    if (m_multiplicity == 4) {
        const auto dt = named(invec[3], "drift");
        if (dt && dt->shape().size() == 2 && dt->shape()[0] > 0) {
            const double* D = (const double*) dt->data();
            for (size_t r = 0; r < dt->shape()[0]; ++r) drift[(int) D[3 * r]] = {D[3 * r + 1], D[3 * r + 2]};
        }
    }

    // clusters (ql_m2m_proto.collect_one 141-154, ql08.cl_one 88), join tree blobs -> table rows
    std::vector<L::Cluster> cl;
    std::vector<double> mu, membership;
    size_t nmiss = 0, ntree = 0;
    const auto joined = L::join_clusters(*grouping, B, nbt, ncolb, ntree, nmiss);
    const size_t ntree_clusters = joined.size();
    for (const auto& jc : joined) {
        const auto ci = jc.index;
        const auto& rows = jc.rows;
        for (int r : rows) {
            membership.push_back((double) ci);
            membership.push_back(B[r * ncolb]);
        }
        if (rows.empty()) continue;
        double Q = 0;   // np.bincount: sequential from 0
        for (int r : rows) Q += B[r * ncolb + 14];
        if (!(Q >= m_C.qmin_cluster)) continue;
        std::vector<std::array<double, 3>> P;
        std::vector<double> pq, px;
        for (int r : rows) {
            const double* b = B + r * ncolb;
            P.push_back({b[11], b[12], b[13]});
            pq.push_back(std::max(b[14], 0.0));
            px.push_back(b[11] * pq.back());
        }
        const double spq = L::np_sum(pq.data(), pq.size());
        if (spq <= 0) continue;
        L::Cluster c;
        c.index = (int) ci;
        c.rep = (int) B[rows.front() * ncolb];
        c.Q = Q;
        c.nblob = (int) rows.size();
        c.x_app = L::np_sum(px.data(), px.size()) / spq;
        if (!L::shift_table(*m_lib, m_C, P, pq, c)) continue;
        const auto d = drift.find(c.rep);
        mu.push_back(d == drift.end() ? std::numeric_limits<double>::quiet_NaN() : d->second.first);
        cl.push_back(std::move(c));
    }
    if (nmiss) log->warn("call={} ident={}: {} of {} tree blobs not in the blob table", m_count, ident, nmiss, ntree);

    // port 2: flashes -> groups -> features -> decisions
    const auto ft = named(invec[2], "opflash");
    std::vector<double> ftime, fpe;
    if (ft && ft->shape().size() == 2 && ft->shape()[0] > 0) {
        if (ft->shape()[1] != (size_t) L::NCH + 1) raise<ValueError>("FdvdLowEQLMatching: opflash has %d columns", (int) ft->shape()[1]);
        const double* Fp = (const double*) ft->data();
        for (size_t r = 0; r < ft->shape()[0]; ++r) {
            ftime.push_back(Fp[r * (L::NCH + 1)]);
            fpe.insert(fpe.end(), Fp + r * (L::NCH + 1) + 1, Fp + (r + 1) * (L::NCH + 1));
        }
    }
    const auto G = L::build_groups(ftime, fpe, m_C);
    const auto Fe = L::features(cl, G, m_cal, m_C, m_pd_pos);
    std::vector<std::vector<bool>> masks;
    std::vector<double> dec;
    for (size_t s = 0; s < m_veto_k.size(); ++s) {
        masks.push_back(L::drift_mask(Fe, m_veto_k[s], mu, cl, m_cal, m_C));
        for (size_t a = 0; a < m_arms.size(); ++a) {
            for (const auto& [k, g] : L::decide(Fe, masks.back(), cl, m_arms[a], m_cal)) {
                dec.insert(dec.end(), {(double) a, (double) s, (double) k, (double) g});
            }
        }
    }

    // output tables
    std::vector<double> ctab, gtab, ptab;
    std::vector<float> tabs;
    std::vector<double> tidx;
    for (size_t k = 0; k < cl.size(); ++k) {
        const auto& c = cl[k];
        const auto d = drift.find(c.rep);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        ctab.insert(ctab.end(), {(double) c.index, (double) c.rep, (double) c.nblob, c.Q, c.x_app, c.t_lo, c.t_hi, c.dt,
                                 d == drift.end() ? nan : d->second.first, d == drift.end() ? nan : d->second.second});
        if (m_dump_tables) {
            tidx.insert(tidx.end(), {(double) (tabs.size() / L::NCH), (double) c.nstep()});
            tabs.insert(tabs.end(), c.tab.begin(), c.tab.end());
        }
    }
    for (size_t i = 0; i < Fe.gt.size(); ++i) {
        gtab.insert(gtab.end(), {(double) i, (double) Fe.gidx[i], Fe.gt[i], Fe.gtot[i], (double) Fe.gnpd[i]});
        for (int ch = 0; ch < L::NCH; ++ch) gtab.push_back((double) G.pe[Fe.gidx[i] * L::NCH + ch]);
    }
    const size_t ns = m_veto_k.size();
    for (size_t p = 0; p < Fe.k.size(); ++p) {
        ptab.insert(ptab.end(), {(double) Fe.k[p], (double) Fe.g[p], (double) Fe.r[p], (double) Fe.ks[p],
                                 (double) Fe.dc[p], (double) Fe.npp[p]});
        for (size_t s = 0; s < ns; ++s) ptab.push_back(masks[s][p] ? 1.0 : 0.0);
    }
    auto tv = std::make_shared<ITensor::vector>();
    tv->push_back(tensor("clusters", ctab, 10));
    tv->push_back(tensor("membership", membership, 2));
    tv->push_back(tensor("groups", gtab, 5 + L::NCH));
    tv->push_back(tensor("pairs", ptab, 6 + ns));
    tv->push_back(tensor("decisions", dec, 4));
    if (m_drift) {
        tv->push_back(tensor("drift", dres.rows, 3));
        if (m_drift->dump_crops()) {
            Configuration cmd;
            cmd["name"] = "crops";
            tv->push_back(std::make_shared<Aux::SimpleTensor>(
                ITensor::shape_t{dres.ncrop, (size_t) m_drift->nch(), (size_t) m_drift->ntk()}, dres.crops.data(), cmd));
        }
    }
    if (m_dump_tables) {
        tv->push_back(tensor("tables", tabs, L::NCH));
        tv->push_back(tensor("table_index", tidx, 2));
        tv->push_back(bt);   // the port-1 blob table, for the doc 16 gates
    }
    Configuration md;
    md["producer"] = "FdvdLowEQLMatching";
    md["ident"] = ident;
    for (size_t a = 0; a < m_arms.size(); ++a) md["arms"][(int) a] = m_arms[a].name;
    for (size_t s = 0; s < ns; ++s) md["specs"][(int) s] = m_spec_names[s];
    md["n_tree_clusters"] = (int) ntree_clusters;
    md["n_tree_blobs"] = (Json::UInt64) ntree;
    md["n_tree_blobs_unjoined"] = (Json::UInt64) nmiss;
    md["n_table_blobs"] = (Json::UInt64) nbt;
    md["n_flash"] = (int) ftime.size();
    md["n_groups_stored"] = (int) G.t.size();
    out = std::make_shared<Aux::SimpleTensorSet>(ident, md, ITensor::shared_vector(tv));
    log->debug("call={} ident={} clusters {} (tree {}), groups {} in window, pairs {}, decisions {}", m_count, ident,
               cl.size(), ntree_clusters, Fe.gt.size(), Fe.k.size(), dec.size() / 4);
    ++m_count;
    return true;
}
