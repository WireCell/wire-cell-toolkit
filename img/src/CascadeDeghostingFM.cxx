// CascadeDeghostingFM (wcfm doc 30): CascadeDeghosting (doc 14) + the stage-A head-score columns at the final level.
// FORKED BY DUPLICATION from CascadeDeghosting.cxx, which stays untouched.  See WireCellImg/CascadeDeghostingFM.h.

#include "WireCellImg/CascadeDeghostingFM.h"
#include "WireCellImg/CascadeGraph.h"
#include "WireCellImg/CellSteiner.h"
#include "WireCellImg/GeomClusteringUtil.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/MemUsage.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/cnpy.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <unordered_map>

WIRECELL_FACTORY(CascadeDeghostingFM, WireCell::Img::CascadeDeghostingFM,
                 WireCell::INamed,
                 WireCell::Img::IClusterFrameTensorJoin, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

Img::CascadeDeghostingFM::CascadeDeghostingFM()
  : Aux::Logger("CascadeDeghostingFM", "img")
{
}
Img::CascadeDeghostingFM::~CascadeDeghostingFM() {}

WireCell::Configuration Img::CascadeDeghostingFM::default_configuration() const
{
    Configuration cfg;
    cfg["levels"] = Json::arrayValue;
    cfg["cut_max_depth"] = m_cut_max_depth;
    cfg["cut_min_length"] = m_cut_min_length;
    cfg["cut_nudge"] = m_cut_nudge;
    cfg["guard"] = m_guard;
    cfg["max_level_nodes"] = m_max_level_nodes;
    cfg["repair"] = m_repair;
    cfg["repair_p_term"] = m_repair_p_term;
    cfg["repair_q_floor"] = m_repair_q_floor;
    cfg["repair_budget"] = m_repair_budget;
    cfg["policy"] = m_policy;
    cfg["charge_tag"] = m_charge_tag;
    cfg["charge_scale"] = m_charge_scale;
    cfg["uncer_cut"] = m_uncer_cut;
    cfg["ident_base"] = m_ident_base;
    cfg["dump_dir"] = m_dump_dir;
    cfg["nthreads"] = m_nthreads;
    cfg["iso_fallback"] = Json::nullValue;   // off; see the header for the members
    cfg["final_guard"] = m_final_guard;
    cfg["keep_slices"] = m_keep_slices;
    cfg["head"] = Json::nullValue;   // doc 30: off; {forward, fm_dim} turns the head-score columns on
    return cfg;
}

void Img::CascadeDeghostingFM::configure(const WireCell::Configuration& cfg)
{
    m_cut_max_depth = get(cfg, "cut_max_depth", m_cut_max_depth);
    m_cut_min_length = get(cfg, "cut_min_length", m_cut_min_length);
    m_cut_nudge = get(cfg, "cut_nudge", m_cut_nudge);
    m_guard = get(cfg, "guard", m_guard);
    m_max_level_nodes = get(cfg, "max_level_nodes", m_max_level_nodes);
    m_repair = get(cfg, "repair", m_repair);
    m_repair_p_term = get(cfg, "repair_p_term", m_repair_p_term);
    m_repair_q_floor = get(cfg, "repair_q_floor", m_repair_q_floor);
    m_repair_budget = get(cfg, "repair_budget", m_repair_budget);
    m_policy = get(cfg, "policy", m_policy);
    m_charge_tag = get(cfg, "charge_tag", m_charge_tag);
    m_charge_scale = get(cfg, "charge_scale", m_charge_scale);
    m_uncer_cut = get(cfg, "uncer_cut", m_uncer_cut);
    m_ident_base = get(cfg, "ident_base", m_ident_base);
    m_dump_dir = get(cfg, "dump_dir", m_dump_dir);
    m_nthreads = std::max(1, get(cfg, "nthreads", m_nthreads));
    m_final_guard = get(cfg, "final_guard", m_final_guard);
    m_keep_slices = get(cfg, "keep_slices", m_keep_slices);
    if (m_final_guard || m_keep_slices) {
        log->debug("wcfm doc 22 knobs: final_guard={} keep_slices={}", m_final_guard, m_keep_slices);
    }
    const auto& jiso = cfg["iso_fallback"];
    m_iso = jiso.isObject();
    if (m_iso) {
        m_iso_nmin = get(jiso, "nmin", m_iso_nmin);
        m_iso_mmin = get(jiso, "mmin", m_iso_mmin);
        m_iso_amin = get(jiso, "amin", m_iso_amin);
        m_iso_t = get(jiso, "t_keep", m_iso_t);
        m_iso_amb_lo = get(jiso, "amb_lo", m_iso_amb_lo);
        m_iso_amb_hi = get(jiso, "amb_hi", m_iso_amb_hi);
        log->debug("iso_fallback on: nmin={} mmin={} amin={} t_keep={} amb=({}, {})", m_iso_nmin, m_iso_mmin, m_iso_amin,
                   m_iso_t, m_iso_amb_lo, m_iso_amb_hi);
    }

    const auto& jh = cfg["head"];
    m_head = jh.isObject();
    if (m_head) {
        m_head_tn = get<std::string>(jh, "forward", "");
        m_fm_dim = get(jh, "fm_dim", m_fm_dim);
        m_head_chunk = get(jh, "chunk", m_head_chunk);
        if (m_head_tn.empty() || m_fm_dim <= 0) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFM: head needs a forward and fm_dim > 0"});
        }
        m_head_forward = Factory::find_tn<ITensorForward>(m_head_tn);
        log->debug("doc 30 head-score columns on: forward={} fm_dim={} chunk={}", m_head_tn, m_fm_dim, m_head_chunk);
    }

    m_levels.clear();
    for (const auto& jl : cfg["levels"]) {
        LevelCfg lc;
        lc.width = get<int>(jl, "width", 0);
        lc.forward_tn = get<std::string>(jl, "forward", "");
        lc.superwire = get<int>(jl, "superwire", 1);
        lc.threshold = get<double>(jl, "threshold", 0.0);
        if (lc.forward_tn.empty()) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFM: every level needs a forward"});
        }
        if (m_levels.empty() != (lc.width <= 0)) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFM: level 0 is uncut (width 0), every later level has a width"});
        }
        lc.forward = Factory::find_tn<ITensorForward>(lc.forward_tn);
        m_levels.push_back(lc);
    }
    if (m_levels.empty()) {
        THROW(ValueError() << errmsg{"CascadeDeghostingFM: no levels"});
    }
    if (m_repair_budget <= 0) {
        const double t = m_levels.back().threshold;
        m_repair_budget = 3.0 * (-std::log(0.2)) + std::log1p(std::exp(-t));  // -log sigmoid(t)
    }
    std::string desc;
    for (const auto& lc : m_levels) {
        desc += String::format(" [width=%d k=%d thr=%.4f %s]", lc.width, lc.superwire, lc.threshold,
                               lc.forward_tn.c_str());
    }
    log->debug("levels:{} cut max_depth={} min_length={} guard={} max_level_nodes={} repair={} budget={:.4f} policy={} charge_scale={} nthreads={}", desc,
               m_cut_max_depth, m_cut_min_length, m_guard, m_max_level_nodes, m_repair, m_repair_budget, m_policy, m_charge_scale,
               m_nthreads);
}

namespace {
    using Cascade::Level;

    // A non-owning ITensor over one of the level's arrays (wcfm doc 15 round 4).  The forward only reads its inputs
    // (the ITensorForward contract; TorchTensorSetService wraps them with torch::from_blob, no copy) and the Level
    // outlives the call, so the level graph is no longer duplicated into SimpleTensor stores during the forward
    // (0.22 GB of the 587/2 4-wire peak, jemalloc profile).
    template <typename T>
    class ViewTensor : public ITensor {
      public:
        ViewTensor(const std::vector<T>& v, const ITensor::shape_t& shape)
          : m_data(v.empty() ? nullptr : reinterpret_cast<const std::byte*>(v.data()))
          , m_nbytes(v.size() * sizeof(T))
          , m_shape(shape)
        {
        }
        virtual ~ViewTensor() {}
        virtual const std::type_info& element_type() const { return typeid(T); }
        virtual size_t element_size() const { return sizeof(T); }
        virtual shape_t shape() const { return m_shape; }
        virtual const std::byte* data() const { return m_data; }
        virtual size_t size() const { return m_nbytes; }

      private:
        const std::byte* m_data;
        size_t m_nbytes;
        shape_t m_shape;
    };

    template <typename T>
    ITensor::pointer tens(const std::vector<T>& v, const ITensor::shape_t& shape)
    {
        return std::make_shared<ViewTensor<T>>(v, shape);
    }

    // (logit, qhat) of every node of the level; xb = the node features (N x ncol: 15, or 18 with the head columns)
    void run_forward(const ITensorForward::pointer& fwd, const Level& lev, const std::vector<float>& xb, size_t ncol,
                     int ident, std::vector<float>& logit, std::vector<float>& qhat)
    {
        const size_t N = lev.nnodes(), M = lev.wq.size(), E = lev.bw_src.size();
        if (xb.size() != N * ncol) {
            THROW(RuntimeError() << errmsg{"CascadeDeghostingFM: node feature array size mismatch"});
        }
        ITensor::vector tv{
            tens(xb, {N, ncol}),          tens(lev.wq, {M}),        tens(lev.wplane, {M}),
            tens(lev.bw_src, {E}),        tens(lev.bw_dst, {E}),    tens(lev.bw_w, {E}),
            tens(lev.bb, {lev.bb.size() / 2, 2}), tens(lev.bb_in, {lev.bb_in.size() / 2, 2}),
            tens(lev.ww, {lev.ww.size() / 2, 2})};
        auto in = std::make_shared<Aux::SimpleTensorSet>(ident, Configuration(), std::make_shared<ITensor::vector>(tv));
        auto out = fwd->forward(in);
        if (!out) {
            THROW(RuntimeError() << errmsg{"CascadeDeghostingFM: forward failed"});
        }
        auto ot = out->tensors();
        if (ot->size() < 2) {
            THROW(RuntimeError() << errmsg{"CascadeDeghostingFM: forward must return (logit, qhat)"});
        }
        for (int k = 0; k < 2; ++k) {
            const auto& t = ot->at(k);
            if (t->dtype() != "f4" || t->size() != N * sizeof(float)) {
                THROW(RuntimeError() << errmsg{String::format("CascadeDeghostingFM: output %d is %s of %d bytes, want f4 [%d]",
                                                              k, t->dtype(), (int) t->size(), (int) N)});
            }
            const float* p = reinterpret_cast<const float*>(t->data());
            (k == 0 ? logit : qhat).assign(p, p + N);
        }
    }

    // ---- doc 30: the FM feature tensor set (FMFeatureExtract output) and the head-score columns

    // IEEE 754 binary16 bits -> float (the training graphs' b_fm came from half-stored sidecars; img does not link torch)
    float half_to_float(uint16_t h)
    {
        const uint32_t sign = (uint32_t) (h & 0x8000) << 16;
        const uint32_t exp = (h >> 10) & 0x1f;
        const uint32_t man = h & 0x3ff;
        uint32_t bits;
        if (exp == 0) {
            if (man == 0) {
                bits = sign;
            }
            else {   // subnormal: normalise
                uint32_t m = man;
                int e = -1;
                do {
                    m <<= 1;
                    ++e;
                } while ((m & 0x400) == 0);
                bits = sign | ((uint32_t) (127 - 15 - e) << 23) | ((m & 0x3ff) << 13);
            }
        }
        else if (exp == 31) {
            bits = sign | 0x7f800000u | (man << 13);
        }
        else {
            bits = sign | ((exp + 127 - 15) << 23) | (man << 13);
        }
        float f;
        std::memcpy(&f, &bits, sizeof f);
        return f;
    }
    // float -> binary16 bits, round to nearest even (c10::Half's rounding), so e_p matches the f16 b_fm of the graphs
    uint16_t float_to_half(float f)
    {
        uint32_t x;
        std::memcpy(&x, &f, sizeof x);
        const uint32_t sign = (x >> 16) & 0x8000;
        const int32_t exp = (int32_t) ((x >> 23) & 0xff) - 127 + 15;
        uint32_t man = x & 0x7fffff;
        if (((x >> 23) & 0xff) == 0xff) {   // inf / nan
            return (uint16_t) (sign | 0x7c00 | (man ? 0x200 : 0));
        }
        if (exp >= 31) return (uint16_t) (sign | 0x7c00);
        if (exp <= 0) {
            if (exp < -10) return (uint16_t) sign;
            man |= 0x800000;
            const int shift = 14 - exp;
            uint32_t half_man = man >> shift;
            const uint32_t rem = man & ((1u << shift) - 1);
            const uint32_t halfway = 1u << (shift - 1);
            if (rem > halfway || (rem == halfway && (half_man & 1))) ++half_man;
            return (uint16_t) (sign | half_man);
        }
        uint32_t half_man = man >> 13;
        const uint32_t rem = man & 0x1fff;
        uint32_t bits = sign | ((uint32_t) exp << 10) | half_man;
        if (rem > 0x1000 || (rem == 0x1000 && (half_man & 1))) ++bits;   // may carry into the exponent: correct
        return (uint16_t) bits;
    }
    float round_half(float f) { return half_to_float(float_to_half(f)); }

    // per plane: (channel ident, slice index) -> row, and the rows' features as float
    struct FMIndex {
        int dim{0};
        std::array<std::unordered_map<int64_t, int>, 3> rows;
        std::array<std::vector<float>, 3> feat;
        static int64_t key(int ch, int k) { return (int64_t) ch * 1000000LL + k; }
        bool has(int p, int ch, int k) const { return rows[p].count(key(ch, k)) > 0; }
        const float* row(int p, int ch, int k) const
        {
            auto it = rows[p].find(key(ch, k));
            return it == rows[p].end() ? nullptr : feat[p].data() + (size_t) it->second * dim;
        }
    };

    FMIndex parse_fm(const ITensorSet::pointer& fm, int dim, size_t& npix)
    {
        FMIndex idx;
        idx.dim = dim;
        npix = 0;
        if (!fm || !fm->tensors()) return idx;
        const auto& tv = *fm->tensors();
        // tensors come in (coords, feat) pairs per plane, each tagged with metadata name and plane
        std::array<ITensor::pointer, 3> coords{}, feats{};
        for (const auto& t : tv) {
            const auto md = t->metadata();
            const int p = md.isMember("plane") ? md["plane"].asInt() : -1;
            const std::string name = md.isMember("name") ? md["name"].asString() : "";
            if (p < 0 || p > 2) continue;
            if (name == "coords") coords[p] = t;
            else if (name == "feat" || name == "feat_half") feats[p] = t;
        }
        for (int p = 0; p < 3; ++p) {
            if (!coords[p] || !feats[p]) continue;
            const auto cs = coords[p]->shape(), fs = feats[p]->shape();
            if (cs.size() != 2 || cs[1] != 2 || coords[p]->dtype() != "i4" || fs.size() != 2 || fs[0] != cs[0] ||
                (int) fs[1] != dim) {
                THROW(ValueError() << errmsg{String::format("CascadeDeghostingFM: FM tensors of plane %d have an unexpected shape or dtype", p)});
            }
            const size_t N = cs[0];
            const int32_t* c = reinterpret_cast<const int32_t*>(coords[p]->data());
            idx.feat[p].resize(N * dim);
            if (feats[p]->dtype() == "f4") {
                const float* f = reinterpret_cast<const float*>(feats[p]->data());
                for (size_t i = 0; i < N * dim; ++i) idx.feat[p][i] = round_half(f[i]);   // as the half-stored sidecars
            }
            else if (feats[p]->dtype() == "u2") {
                const uint16_t* f = reinterpret_cast<const uint16_t*>(feats[p]->data());
                for (size_t i = 0; i < N * dim; ++i) idx.feat[p][i] = half_to_float(f[i]);
            }
            else {
                THROW(ValueError() << errmsg{"CascadeDeghostingFM: FM feature tensor is neither f4 nor u2"});
            }
            idx.rows[p].reserve(N);
            for (size_t i = 0; i < N; ++i) idx.rows[p][FMIndex::key(c[2 * i], c[2 * i + 1])] = (int) i;
            npix += N;
        }
        return idx;
    }

    struct HeadStats {
        size_t nhas[3]{0, 0, 0};
        size_t nfilled{0};
        size_t nactive_wires{0};
        float smin{0}, smax{0};
    };

    // The three columns [score, ls_mean, ls_min] per node of a k = 1 level (wcfm scripts/d27_head.py graph_inputs /
    // head_scores / fill_nan, and gnn_dataset.py b_fm).  The slice index of the FM coords is the slice ident.
    std::vector<float> head_columns(const ITensorForward::pointer& head, const Level& lev, const Cascade::SliceCharge& sc,
                                    const FMIndex& fm, int ident, HeadStats& st, int nthreads, size_t chunk,
                                    double tsplit[2])
    {
        using clock = std::chrono::steady_clock;
        const size_t N = lev.nnodes(), NW = lev.wq.size(), E = lev.bw_src.size();
        const int D = fm.dim;
        // wcfm doc 34: the head's inputs are built and forwarded `chunk` cells at a time (0 = the whole level in one
        // forward, as doc 30), and a chunk's rows are filled in nthreads blocks: a row depends on its own cell only
        // and is written to its own slot, so the arrays do not depend on the threading.
        const size_t C = chunk > 0 ? std::min(chunk, std::max<size_t>(N, 1)) : std::max<size_t>(N, 1);
        std::array<std::vector<float>, 3> e;
        std::vector<float> has;
        std::vector<float> score(N, 0.0f);
        tsplit[0] = tsplit[1] = 0.0;
        size_t a = 0;
        do {
            const size_t b = std::min(N, a + C), n = b - a;
            const auto tc = clock::now();
            for (int p = 0; p < 3; ++p) e[p].assign(n * D, 0.0f);
            has.assign(n * 3, 0.0f);
            Cascade::parallel_blocks(n, nthreads, [&](size_t j0, size_t j1) {
                std::vector<double> acc(D);
                for (size_t j = j0; j < j1; ++j) {
                    const size_t i = a + j;
                    const auto ch = Cascade::blob_channels(lev.blobs[i]);
                    const int k = sc.slice_of[lev.sidx[i]]->ident();
                    for (int p = 0; p < 3; ++p) {
                        std::fill(acc.begin(), acc.end(), 0.0);
                        int nrow = 0;
                        for (int c : ch[p]) {
                            const float* r = fm.row(p, c, k);
                            if (!r) continue;
                            for (int d = 0; d < D; ++d) acc[d] += r[d];
                            ++nrow;
                        }
                        if (nrow) {
                            has[j * 3 + p] = 1.0f;
                            float* out = e[p].data() + j * D;
                            for (int d = 0; d < D; ++d) out[d] = round_half((float) (acc[d] / nrow));   // b_fm is stored f16
                        }
                    }
                }
            });
            for (size_t j = 0; j < n; ++j) {
                for (int p = 0; p < 3; ++p) st.nhas[p] += has[j * 3 + p] > 0.0f;
            }
            const auto tf = clock::now();
            // the head: score = CrossHead(eU, eV, eW, has)
            {
                ITensor::vector tv{tens(e[0], {n, (size_t) D}), tens(e[1], {n, (size_t) D}), tens(e[2], {n, (size_t) D}),
                                   tens(has, {n, 3})};
                auto in = std::make_shared<Aux::SimpleTensorSet>(ident, Configuration(), std::make_shared<ITensor::vector>(tv));
                auto out = head->forward(in);
                if (!out || !out->tensors() || out->tensors()->empty()) {
                    THROW(RuntimeError() << errmsg{"CascadeDeghostingFM: head forward failed"});
                }
                const auto& t = out->tensors()->front();
                if (t->dtype() != "f4" || t->size() != n * sizeof(float)) {
                    THROW(RuntimeError() << errmsg{String::format("CascadeDeghostingFM: head output is %s of %d bytes, want f4 [%d]",
                                                                  t->dtype(), (int) t->size(), (int) n)});
                }
                const float* p = reinterpret_cast<const float*>(t->data());
                std::copy(p, p + n, score.begin() + a);
            }
            tsplit[0] += std::chrono::duration<double>(tf - tc).count();
            tsplit[1] += std::chrono::duration<double>(clock::now() - tf).count();
            a = b;
        } while (a < N);
        for (auto& v : e) std::vector<float>().swap(v);
        std::vector<float>().swap(has);
        st.smin = N ? *std::min_element(score.begin(), score.end()) : 0.0f;
        st.smax = N ? *std::max_element(score.begin(), score.end()) : 0.0f;
        // active wire nodes: an FM pixel at (plane, channel, slice ident)
        std::vector<char> wact(NW, 0);
        for (size_t w = 0; w < NW; ++w) {
            const int k = sc.slice_of[lev.wsidx[w]]->ident();
            wact[w] = fm.has((int) lev.wplane[w], lev.wchan[w], k) ? 1 : 0;
            st.nactive_wires += wact[w];
        }
        // per active wire: log-sum-exp of the scores of its nodes (double, max-shifted), over all bw edges to it
        std::vector<double> wmax(NW, -1e300), wsum(NW, 0.0);
        for (size_t ed = 0; ed < E; ++ed) {
            const int64_t w = lev.bw_dst[ed];
            if (!wact[w]) continue;
            wmax[w] = std::max(wmax[w], (double) score[lev.bw_src[ed]]);
        }
        for (size_t ed = 0; ed < E; ++ed) {
            const int64_t w = lev.bw_dst[ed];
            if (!wact[w]) continue;
            wsum[w] += std::exp((double) score[lev.bw_src[ed]] - wmax[w]);
        }
        std::vector<double> tot(N, 0.0), mn(N, 1e300);
        std::vector<int> cnt(N, 0);
        for (size_t ed = 0; ed < E; ++ed) {
            const int64_t w = lev.bw_dst[ed];
            if (!wact[w]) continue;
            const int64_t i = lev.bw_src[ed];
            const double ls = (double) score[i] - (wmax[w] + std::log(wsum[w]));
            tot[i] += ls;
            mn[i] = std::min(mn[i], ls);
            ++cnt[i];
        }
        std::vector<float> cols(N * 3);
        float fmin_mean = std::numeric_limits<float>::infinity(), fmin_min = std::numeric_limits<float>::infinity();
        for (size_t i = 0; i < N; ++i) {
            cols[i * 3] = score[i];
            if (cnt[i]) {
                cols[i * 3 + 1] = (float) (tot[i] / cnt[i]);
                cols[i * 3 + 2] = (float) mn[i];
                fmin_mean = std::min(fmin_mean, cols[i * 3 + 1]);
                fmin_min = std::min(fmin_min, cols[i * 3 + 2]);
            }
        }
        // NaN fill (doc 27 sec 1): a node with no active wire takes the level's column minimum over the finite nodes
        if (!std::isfinite(fmin_mean)) fmin_mean = 0.0f;
        if (!std::isfinite(fmin_min)) fmin_min = 0.0f;
        for (size_t i = 0; i < N; ++i) {
            if (cnt[i]) continue;
            cols[i * 3 + 1] = fmin_mean;
            cols[i * 3 + 2] = fmin_min;
            ++st.nfilled;
        }
        return cols;
    }

    std::vector<std::array<int64_t, 2>> pairs_of(const std::vector<int64_t>& flat)
    {
        std::vector<std::array<int64_t, 2>> out(flat.size() / 2);
        for (size_t i = 0; i < out.size(); ++i) out[i] = {flat[2 * i], flat[2 * i + 1]};
        return out;
    }

    // BlobClustering.cxx add_slice / add_blobs, duplicated (fork rule: the production file is untouched)
    void add_slice(cluster_indexed_graph_t& grind, const ISlice::pointer& islice)
    {
        if (grind.has(islice)) {
            return;
        }
        for (const auto& ichv : islice->activity()) {
            const IChannel::pointer ich = ichv.first;
            if (grind.has(ich)) {
                continue;
            }
            for (const auto& iwire : ich->wires()) {
                grind.edge(ich, iwire);
            }
        }
    }
    void add_blobs(cluster_indexed_graph_t& grind, const IBlob::vector& iblobs)
    {
        for (const auto& iblob : iblobs) {
            auto islice = iblob->slice();
            add_slice(grind, islice);
            grind.edge(islice, iblob);
            auto iface = iblob->face();
            auto wire_planes = iface->planes();
            const auto& shape = iblob->shape();
            for (const auto& strip : shape.strips()) {
                const int num_nonplane_layers = 2;
                int iplane = strip.layer - num_nonplane_layers;
                if (iplane < 0) {
                    continue;
                }
                const auto& wires = wire_planes[iplane]->wires();
                for (int wip = strip.bounds.first; wip < strip.bounds.second and wip < int(wires.size()); ++wip) {
                    grind.edge(iblob, wires[wip]);
                }
            }
        }
    }
}  // namespace

bool Img::CascadeDeghostingFM::operator()(const input_tuple_type& intup, output_pointer& out)
{
    out = nullptr;
    const auto& in = std::get<0>(intup);
    const auto& frame = std::get<1>(intup);
    const auto& fmset = std::get<2>(intup);
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    const auto& gr = in->graph();

    bool use_frame = false;
    if (frame) {
        const size_t ntr = m_charge_tag.empty() ? frame->traces()->size() : frame->tagged_traces(m_charge_tag).size();
        use_frame = ntr > 0;
        if (!use_frame) log->warn("call={} frame {} has no trace tagged \"{}\": charge from the slice activity", m_count,
                                  frame->ident(), m_charge_tag);
    }
    else {
        log->warn("call={} no frame: charge from the slice activity (not the training charge)", m_count);
    }
    FMIndex fmidx;
    if (m_head) {
        size_t npix = 0;
        fmidx = parse_fm(fmset, m_fm_dim, npix);
        if (!fmset) log->warn("call={} head on but no FM tensor set on port 2: no cell has an active FM pixel", m_count);
        log->debug("call={} FM pixels={} (planes {}/{}/{}) dim={}", m_count, npix, fmidx.rows[0].size(), fmidx.rows[1].size(),
                   fmidx.rows[2].size(), m_fm_dim);
    }
    auto sc = use_frame ? Cascade::make_slice_charge_frame(gr, frame, m_charge_tag, m_charge_scale)
                        : Cascade::make_slice_charge(gr, m_charge_scale, m_uncer_cut);
    std::vector<IBlob::pointer> cur;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        if (gr[vtx].code() == 'b') cur.push_back(std::get<IBlob::pointer>(gr[vtx].ptr));
    }
    const size_t nin = cur.size();
    auto uni = Cascade::make_universe(cur, sc);
    std::vector<int> parent(cur.size(), -1);   // index in the previous level
    std::vector<int> depth(cur.size(), 0);     // bisections from the uncut blob
    Cascade::CutParams cpar;
    cpar.nudge = m_cut_nudge;
    cpar.max_depth = m_cut_max_depth;
    cpar.min_length = m_cut_min_length;
    int next_ident = m_ident_base;
    std::vector<IBlob::pointer> kept;
    Cascade::SteinerResult srep;
    size_t nkeep_thr = 0;

    const size_t nlev = m_levels.size();
    for (size_t li = 0; li < nlev; ++li) {
        const auto& lc = m_levels[li];
        const auto tl = clock::now();
        if (cur.size() < 2) {   // d12_cascade: a graph with < 2 nodes is empty, nothing survives
            log->debug("call={} cluster={} level={} nodes={}: < 2, nothing kept", m_count, in->ident(), li, cur.size());
            cur.clear();
            break;
        }
        auto lev = Cascade::build_level(cur, sc, uni, lc.superwire, m_policy, m_nthreads);
        const auto tb = clock::now();
        const double rss_build = memusage_resident();
        std::vector<float> logit, qhat;
        // doc 30: at the final level (k = 1) with the head on, xb gains [score, ls_mean, ls_min]
        std::vector<float> headcols;
        const bool with_head = m_head && li + 1 == nlev;
        size_t ncol = 15;
        std::vector<float> xb18;
        if (with_head) {
            HeadStats hs;
            double tsplit[2];
            headcols = head_columns(m_head_forward, lev, sc, fmidx, in->ident(), hs, m_nthreads, (size_t) m_head_chunk, tsplit);
            fmidx = FMIndex();   // doc 34: the pixel table is not needed past the head columns
            const size_t Nh = lev.nnodes();
            xb18.resize(Nh * 18);
            for (size_t i = 0; i < Nh; ++i) {
                std::copy(lev.xb.begin() + i * 15, lev.xb.begin() + (i + 1) * 15, xb18.begin() + i * 18);
                std::copy(headcols.begin() + i * 3, headcols.begin() + (i + 1) * 3, xb18.begin() + i * 18 + 15);
            }
            ncol = 18;
            log->debug("call={} cluster={} level={} head columns: nodes={} has U/V/W={}/{}/{} active wires={}/{} "
                       "no-active-wire filled={} score=[{:.3f}, {:.3f}]",
                       m_count, in->ident(), li, Nh, hs.nhas[0], hs.nhas[1], hs.nhas[2], hs.nactive_wires, lev.wq.size(),
                       hs.nfilled, hs.smin, hs.smax);
            // wcfm doc 34: the final level's forward= phase is these two plus the level model's forward
            log->debug("call={} cluster={} level={} head split: columns={:.3f}s head_forward={:.3f}s chunk={}", m_count,
                       in->ident(), li, tsplit[0], tsplit[1], m_head_chunk);
        }
        run_forward(lc.forward, lev, with_head ? xb18 : lev.xb, ncol, in->ident(), logit, qhat);
        const auto tf = clock::now();
        const double rss_fwd = memusage_resident();
        const size_t N = lev.nnodes();
        std::vector<char> decision(N, 1);
        size_t nguard = 0, ncand = 0, ncap = 0;
        std::vector<IBlob::pointer> next;
        std::vector<int> next_parent, next_depth;
        if (li + 1 < nlev) {
            std::vector<int> cand;
            for (size_t i = 0; i < N; ++i) if (logit[i] < lc.threshold) cand.push_back((int) i);
            std::stable_sort(cand.begin(), cand.end(), [&](int a, int b) { return logit[a] < logit[b]; });
            ncand = cand.size();
            std::vector<bool> pruned(N, false);
            if (m_guard) {
                pruned = Cascade::guard_prune(lev, cand, nguard);
            }
            else {
                for (int b : cand) pruned[b] = true;
            }
            for (size_t i = 0; i < N; ++i) decision[i] = !pruned[i];
            // cut the survivors to the next width, one parent at a time (exact parent map)
            std::vector<int> order;
            for (size_t i = 0; i < N; ++i) if (decision[i]) order.push_back((int) i);
            if (m_max_level_nodes > 0) {
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return logit[a] > logit[b]; });
            }
            struct Kid {
                int parent;
                IBlob::vector blobs;
                std::vector<int> depth;
            };
            std::vector<Kid> kids;
            size_t nnext = 0;
            const int width = m_levels[li + 1].width;
            // the bisection of every survivor that needs it, first (a pure function of the shape; nthreads blocks),
            // then consumed in `order` exactly as before (doc 15 round 3)
            std::vector<std::vector<std::pair<RayGrid::Blob, int>>> cut(order.size());
            Cascade::parallel_blocks(order.size(), m_nthreads, [&](size_t k0, size_t k1) {
                for (size_t k = k0; k < k1; ++k) {
                    const auto& b = cur[order[k]];
                    if (Cascade::needs_cutting(b->shape(), width)) {
                        cut[k] = Cascade::cut_shape(b->face()->raygrid(), b->shape(), depth[order[k]], width, cpar);
                    }
                }
            });
            for (size_t ko = 0; ko < order.size(); ++ko) {
                const int i = order[ko];
                Kid k{i, {}, {}};
                const auto& b = cur[i];
                if (!Cascade::needs_cutting(b->shape(), width)) {   // passes through (as BlobCutting does)
                    k.blobs.push_back(b);
                    k.depth.push_back(depth[i]);
                }
                else {
                    auto leaves = std::move(cut[ko]);
                    const float value = b->value() / leaves.size();   // BlobCutting: parent charge shared equally
                    for (auto& [shape, d] : leaves) {
                        k.blobs.push_back(std::make_shared<Aux::SimpleBlob>(next_ident++, value, b->uncertainty(), shape,
                                                                            b->slice(), b->face()));
                        k.depth.push_back(d);
                    }
                }
                if (m_max_level_nodes > 0 && nnext + k.blobs.size() > (size_t) m_max_level_nodes) {
                    decision[i] = 0;
                    ++ncap;
                    continue;
                }
                nnext += k.blobs.size();
                kids.push_back(std::move(k));
            }
            // back to node order
            std::stable_sort(kids.begin(), kids.end(), [](const Kid& a, const Kid& b) { return a.parent < b.parent; });
            for (auto& k : kids) {
                for (size_t j = 0; j < k.blobs.size(); ++j) {
                    next.push_back(k.blobs[j]);
                    next_parent.push_back(k.parent);
                    next_depth.push_back(k.depth[j]);
                }
            }
        }
        else {
            std::vector<bool> keep(N);
            for (size_t i = 0; i < N; ++i) keep[i] = logit[i] >= lc.threshold;
            nkeep_thr = std::count(keep.begin(), keep.end(), true);
            if (m_repair) {
                auto edges = pairs_of(lev.bb);
                auto ein = pairs_of(lev.bb_in);
                edges.insert(edges.end(), ein.begin(), ein.end());
                Cascade::SteinerParams par;
                par.p_term = m_repair_p_term;
                par.q_floor = m_repair_q_floor;
                par.budget = m_repair_budget;
                srep = Cascade::steiner_repair(N, edges, logit, qhat, keep, par);
                keep = srep.keep;
            }
            if (m_iso) {   // wcfm doc 17: dense ambiguous slices keep their cells down to t_keep
                std::vector<int> fc(N);
                for (size_t i = 0; i < N; ++i) fc[i] = cur[i]->face()->which();
                Cascade::IsoParams ip;
                ip.nmin = (int) m_iso_nmin;
                ip.mmin = m_iso_mmin;
                ip.amin = m_iso_amin;
                ip.t_keep = m_iso_t;
                ip.amb_lo = m_iso_amb_lo;
                ip.amb_hi = m_iso_amb_hi;
                const auto ir = Cascade::iso_fallback(fc, lev.sidx, lev.wq, lev.wsidx, logit, keep, ip);
                log->debug("call={} cluster={} iso_fallback slices={} added={}", m_count, in->ident(), ir.nslices,
                           ir.nadded);
            }
            if (m_final_guard) {   // wcfm doc 22: every charged wire node keeps at least one 3-D explanation
                const size_t ng = Cascade::final_guard(lev, logit, keep);
                log->debug("call={} cluster={} final_guard added={}", m_count, in->ident(), ng);
            }
            for (size_t i = 0; i < N; ++i) {
                decision[i] = keep[i] ? 1 : 0;
                if (keep[i]) kept.push_back(cur[i]);
            }
        }
        const auto te = clock::now();
        const double dt = std::chrono::duration<double>(te - tl).count();
        size_t nsurv = std::count(decision.begin(), decision.end(), 1);
        log->debug("call={} cluster={} level={} k={} nodes={} wires={} bw={} bb={} bb_in={} ww={} cand={} guarded={} "
                   "capped={} survive={} next={} t={:.2f}s",
                   m_count, in->ident(), li, lev.k, N, lev.wq.size(), lev.bw_src.size(), lev.bb.size() / 2,
                   lev.bb_in.size() / 2, lev.ww.size() / 2, ncand, nguard, ncap, nsurv, next.size(), dt);
        // wcfm doc 15: wall per phase (build = wires+features / bb / bb_in / ww), resident memory after each (MB)
        log->debug("call={} cluster={} level={} phases build={:.2f}s ({:.2f} {:.2f} {:.2f} {:.2f}) forward={:.2f}s "
                   "select={:.2f}s rss after build={:.0f} forward={:.0f} select={:.0f} MB",
                   m_count, in->ident(), li, std::chrono::duration<double>(tb - tl).count(), lev.tpart[0],
                   lev.tpart[1], lev.tpart[2], lev.tpart[3], std::chrono::duration<double>(tf - tb).count(),
                   std::chrono::duration<double>(te - tf).count(), rss_build / 1024, rss_fwd / 1024,
                   memusage_resident() / 1024);

        if (!m_dump_dir.empty()) {
            const std::string fn = String::format("%s/cascade-%d-L%d.npz", m_dump_dir.c_str(), in->ident(), (int) li);
            std::vector<int32_t> wip(N * 6), sl(N), fc(N), an(N), bid(N);
            for (size_t i = 0; i < N; ++i) {
                for (const auto& strip : cur[i]->shape().strips()) {
                    const int ip = strip.layer - 2;
                    if (ip < 0 || ip > 2) continue;
                    wip[i * 6 + 2 * ip] = strip.bounds.first;
                    wip[i * 6 + 2 * ip + 1] = strip.bounds.second;
                }
                sl[i] = cur[i]->slice()->ident();
                fc[i] = cur[i]->face()->which();
                an[i] = cur[i]->face()->anode();
                bid[i] = cur[i]->ident();
            }
            std::vector<int32_t> par32(parent.begin(), parent.end());
            std::vector<int32_t> dec(decision.begin(), decision.end());
            std::vector<int32_t> sidx(lev.sidx.begin(), lev.sidx.end());
            cnpy::npz_save(fn, "xb", with_head ? xb18.data() : lev.xb.data(), {N, ncol}, "w");
            if (with_head) {
                cnpy::npz_save(fn, "head_cols", headcols.data(), {N, 3}, "a");
            }
            cnpy::npz_save(fn, "wq", lev.wq.data(), {lev.wq.size()}, "a");
            cnpy::npz_save(fn, "wplane", lev.wplane.data(), {lev.wplane.size()}, "a");
            cnpy::npz_save(fn, "bw_src", lev.bw_src.data(), {lev.bw_src.size()}, "a");
            cnpy::npz_save(fn, "bw_dst", lev.bw_dst.data(), {lev.bw_dst.size()}, "a");
            cnpy::npz_save(fn, "bw_w", lev.bw_w.data(), {lev.bw_w.size()}, "a");
            cnpy::npz_save(fn, "bb", lev.bb.data(), {lev.bb.size() / 2, 2}, "a");
            cnpy::npz_save(fn, "bb_in", lev.bb_in.data(), {lev.bb_in.size() / 2, 2}, "a");
            cnpy::npz_save(fn, "ww", lev.ww.data(), {lev.ww.size() / 2, 2}, "a");
            cnpy::npz_save(fn, "logit", logit.data(), {N}, "a");
            cnpy::npz_save(fn, "qhat", qhat.data(), {N}, "a");
            cnpy::npz_save(fn, "decision", dec.data(), {N}, "a");
            cnpy::npz_save(fn, "parent", par32.data(), {N}, "a");
            std::vector<int32_t> dep32(depth.begin(), depth.end());
            cnpy::npz_save(fn, "depth", dep32.data(), {N}, "a");
            cnpy::npz_save(fn, "wip", wip.data(), {N, 6}, "a");
            cnpy::npz_save(fn, "slice_ident", sl.data(), {N}, "a");
            cnpy::npz_save(fn, "sidx", sidx.data(), {N}, "a");
            cnpy::npz_save(fn, "face", fc.data(), {N}, "a");
            cnpy::npz_save(fn, "anode", an.data(), {N}, "a");
            cnpy::npz_save(fn, "ident", bid.data(), {N}, "a");
            std::vector<int32_t> sid(sc.slice_of.size());
            for (size_t s = 0; s < sid.size(); ++s) sid[s] = sc.slice_of[s]->ident();
            cnpy::npz_save(fn, "slice_ident_of_sidx", sid.data(), {sid.size()}, "a");
            cnpy::npz_save(fn, "wchan", lev.wchan.data(), {lev.wchan.size()}, "a");
            cnpy::npz_save(fn, "wsidx", lev.wsidx.data(), {lev.wsidx.size()}, "a");
        }
        if (li + 1 < nlev) {
            cur = std::move(next);
            parent = std::move(next_parent);
            depth = std::move(next_depth);
        }
    }

    // ---- the output cluster: kept blobs, renumbered, built as BlobClustering builds one
    IBlob::vector outblobs;
    outblobs.reserve(kept.size());
    int ident = m_ident_base;
    for (const auto& b : kept) {
        outblobs.push_back(std::make_shared<Aux::SimpleBlob>(ident++, b->value(), b->uncertainty(), b->shape(),
                                                             b->slice(), b->face()));
    }
    std::vector<IBlob::vector> per(sc.slice_of.size());
    for (const auto& b : outblobs) per[sc.index(b->slice())].push_back(b);
    IBlobSet::vector sets;
    for (size_t s = 0; s < per.size(); ++s) {
        if (per[s].empty()) continue;
        sets.push_back(std::make_shared<Aux::SimpleBlobSet>((int) s, sc.slice_of[s], per[s]));
    }
    cluster_indexed_graph_t grind;
    for (auto it = sets.begin(); it != sets.end(); ++it) {
        add_blobs(grind, (*it)->blobs());
        Img::geom_clustering(grind, it, sets.end(), m_policy);
    }
    if (m_keep_slices) {   // wcfm doc 22: every input slice node of a time whose blobs were all dropped is kept,
                           // with its activity (PointTreeBuilding's ctpc reads every s-node), in input vertex order
        size_t nkept_slices = 0;
        for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
            if (gr[vtx].code() != 's') continue;
            const auto& islice = std::get<ISlice::pointer>(gr[vtx].ptr);
            if (!per[sc.index(islice)].empty() || grind.has(islice)) continue;
            add_slice(grind, islice);
            grind.vertex(islice);
            ++nkept_slices;
        }
        log->debug("call={} cluster={} keep_slices: {} blob-less slice nodes kept", m_count, in->ident(), nkept_slices);
    }
    out = std::make_shared<Aux::SimpleCluster>(std::move(grind.graph()), in->ident());
    const double dt = std::chrono::duration<double>(clock::now() - t0).count();
    log->debug("call={} cluster={} blobs in={} kept at threshold={} weak dropped={} bridges={} added={} out={} "
               "vertices={} edges={} t={:.2f}s",
               m_count, in->ident(), nin, nkeep_thr, srep.nweak_cells, srep.nbridges, srep.nadded, outblobs.size(),
               boost::num_vertices(out->graph()), boost::num_edges(out->graph()), dt);
    ++m_count;
    return true;
}
