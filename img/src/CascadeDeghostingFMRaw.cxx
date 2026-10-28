// CascadeDeghostingFMRaw (wcfm doc 53 part 43 "D3 raw", pr_ml doc 30 sec 9): CascadeDeghostingFM with the head call
// replaced by a PCA projection of the per-view FM means at the final level.
// FORKED BY DUPLICATION from CascadeDeghostingFM.cxx, which stays untouched.  See WireCellImg/CascadeDeghostingFMRaw.h.

#include "WireCellImg/CascadeDeghostingFMRaw.h"
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
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/cnpy.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <unordered_map>

WIRECELL_FACTORY(CascadeDeghostingFMRaw, WireCell::Img::CascadeDeghostingFMRaw,
                 WireCell::INamed,
                 WireCell::Img::IClusterFrameTensorJoin, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

Img::CascadeDeghostingFMRaw::CascadeDeghostingFMRaw()
  : Aux::Logger("CascadeDeghostingFMRaw", "img")
{
}
Img::CascadeDeghostingFMRaw::~CascadeDeghostingFMRaw() {}

WireCell::Configuration Img::CascadeDeghostingFMRaw::default_configuration() const
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
    cfg["slice_start_relative"] = m_slice_start_relative;
    cfg["uncer_cut"] = m_uncer_cut;
    cfg["ident_base"] = m_ident_base;
    cfg["dump_dir"] = m_dump_dir;
    cfg["nthreads"] = m_nthreads;
    cfg["iso_fallback"] = Json::nullValue;   // off; see the header for the members
    cfg["final_guard"] = m_final_guard;
    cfg["keep_slices"] = m_keep_slices;
    cfg["fmraw"] = Json::nullValue;   // off; {pca: "<json on WIRECELL_PATH>", fm_dim: 128} turns the projected FM columns on
    return cfg;
}

void Img::CascadeDeghostingFMRaw::configure(const WireCell::Configuration& cfg)
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
    m_slice_start_relative = get(cfg, "slice_start_relative", m_slice_start_relative);
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

    const auto& jr = cfg["fmraw"];
    m_fmraw = jr.isObject();
    if (m_fmraw) {
        m_pca_file = get<std::string>(jr, "pca", "");
        m_fm_dim = get(jr, "fm_dim", m_fm_dim);
        if (m_pca_file.empty() || m_fm_dim <= 0) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: fmraw needs a pca file and fm_dim > 0"});
        }
        // {mean[dim], basis[dim][k], k, dim}: x_p = ((e_p - mean) @ basis) * has_p
        const auto jp = Persist::load(m_pca_file);
        m_pca_k = jp["k"].asInt();
        if (jp["dim"].asInt() != m_fm_dim || m_pca_k <= 0 || (int) jp["mean"].size() != m_fm_dim ||
            (int) jp["basis"].size() != m_fm_dim) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: pca file shape does not match fm_dim"});
        }
        m_pca_mean.assign(m_fm_dim, 0.0f);
        m_pca_basis.assign((size_t) m_fm_dim * m_pca_k, 0.0f);
        for (int d = 0; d < m_fm_dim; ++d) {
            m_pca_mean[d] = jp["mean"][d].asFloat();
            if ((int) jp["basis"][d].size() != m_pca_k) {
                THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: pca basis row length != k"});
            }
            for (int j = 0; j < m_pca_k; ++j) m_pca_basis[(size_t) d * m_pca_k + j] = jp["basis"][d][j].asFloat();
        }
        log->debug("fmraw columns on: pca={} fm_dim={} k={} (final level xb = 15 + 3 x {} + 3)", m_pca_file, m_fm_dim,
                   m_pca_k, m_pca_k);
    }

    m_levels.clear();
    for (const auto& jl : cfg["levels"]) {
        LevelCfg lc;
        lc.width = get<int>(jl, "width", 0);
        lc.forward_tn = get<std::string>(jl, "forward", "");
        lc.superwire = get<int>(jl, "superwire", 1);
        lc.threshold = get<double>(jl, "threshold", 0.0);
        if (lc.forward_tn.empty()) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: every level needs a forward"});
        }
        if (m_levels.empty() != (lc.width <= 0)) {
            THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: level 0 is uncut (width 0), every later level has a width"});
        }
        lc.forward = Factory::find_tn<ITensorForward>(lc.forward_tn);
        m_levels.push_back(lc);
    }
    if (m_levels.empty()) {
        THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: no levels"});
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
            THROW(RuntimeError() << errmsg{"CascadeDeghostingFMRaw: node feature array size mismatch"});
        }
        ITensor::vector tv{
            tens(xb, {N, ncol}),          tens(lev.wq, {M}),        tens(lev.wplane, {M}),
            tens(lev.bw_src, {E}),        tens(lev.bw_dst, {E}),    tens(lev.bw_w, {E}),
            tens(lev.bb, {lev.bb.size() / 2, 2}), tens(lev.bb_in, {lev.bb_in.size() / 2, 2}),
            tens(lev.ww, {lev.ww.size() / 2, 2})};
        auto in = std::make_shared<Aux::SimpleTensorSet>(ident, Configuration(), std::make_shared<ITensor::vector>(tv));
        auto out = fwd->forward(in);
        if (!out) {
            THROW(RuntimeError() << errmsg{"CascadeDeghostingFMRaw: forward failed"});
        }
        auto ot = out->tensors();
        if (ot->size() < 2) {
            THROW(RuntimeError() << errmsg{"CascadeDeghostingFMRaw: forward must return (logit, qhat)"});
        }
        for (int k = 0; k < 2; ++k) {
            const auto& t = ot->at(k);
            if (t->dtype() != "f4" || t->size() != N * sizeof(float)) {
                THROW(RuntimeError() << errmsg{String::format("CascadeDeghostingFMRaw: output %d is %s of %d bytes, want f4 [%d]",
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
                THROW(ValueError() << errmsg{String::format("CascadeDeghostingFMRaw: FM tensors of plane %d have an unexpected shape or dtype", p)});
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
                THROW(ValueError() << errmsg{"CascadeDeghostingFMRaw: FM feature tensor is neither f4 nor u2"});
            }
            idx.rows[p].reserve(N);
            for (size_t i = 0; i < N; ++i) idx.rows[p][FMIndex::key(c[2 * i], c[2 * i + 1])] = (int) i;
            npix += N;
        }
        return idx;
    }

    struct FMRawStats {
        size_t nhas[3]{0, 0, 0};
        size_t nnone{0};
    };

    // wcfm doc 53 part 43 (D3 raw): per node of a k = 1 level and view p, e_p = the mean of the FM feature rows of the
    // cell's channels with a pixel at (channel, slice ident) -- half-rounded, exactly the head stage's e_p of
    // CascadeDeghostingFM -- and has_p; the columns are [x_0 (k), x_1 (k), x_2 (k), has_0, has_1, has_2] with
    // x_p = ((e_p - mean) @ basis) * has_p (float, the mean subtracted per component, the dot product over d in order).
    // When `means` is non-null it receives e_p (N x 3 x dim) and `has` the flags (N x 3), for the parity dumps.
    std::vector<float> fmraw_columns(const Level& lev, const Cascade::SliceCharge& sc, const FMIndex& fm,
                                     const std::vector<float>& mean, const std::vector<float>& basis, int K,
                                     FMRawStats& st, int nthreads, std::vector<float>* means, std::vector<float>* hasout)
    {
        const size_t N = lev.nnodes();
        const int D = fm.dim;
        const size_t NC = 3 * (size_t) K + 3;
        std::vector<float> cols(N * NC, 0.0f);
        std::vector<float> has(N * 3, 0.0f);
        if (means) means->assign(N * 3 * (size_t) D, 0.0f);
        Cascade::parallel_blocks(N, nthreads, [&](size_t j0, size_t j1) {
            std::vector<double> acc(D);
            std::vector<float> e(D), c(D);
            for (size_t i = j0; i < j1; ++i) {
                const auto ch = Cascade::blob_channels(lev.blobs[i]);
                const int k = sc.slice_of[lev.sidx[i]]->ident();
                for (int p = 0; p < 3; ++p) {
                    std::fill(acc.begin(), acc.end(), 0.0);
                    int nrow = 0;
                    for (int cch : ch[p]) {
                        const float* r = fm.row(p, cch, k);
                        if (!r) continue;
                        for (int d = 0; d < D; ++d) acc[d] += r[d];
                        ++nrow;
                    }
                    if (!nrow) continue;
                    has[i * 3 + p] = 1.0f;
                    for (int d = 0; d < D; ++d) e[d] = round_half((float) (acc[d] / nrow));   // as the head stage
                    if (means) std::copy(e.begin(), e.end(), means->begin() + (i * 3 + p) * D);
                    for (int d = 0; d < D; ++d) c[d] = e[d] - mean[d];
                    float* out = cols.data() + i * NC + (size_t) p * K;
                    for (int j = 0; j < K; ++j) {
                        float s = 0.0f;
                        for (int d = 0; d < D; ++d) s += c[d] * basis[(size_t) d * K + j];
                        out[j] = s;
                    }
                }
                for (int p = 0; p < 3; ++p) cols[i * NC + 3 * (size_t) K + p] = has[i * 3 + p];
            }
        });
        for (size_t i = 0; i < N; ++i) {
            int nv = 0;
            for (int p = 0; p < 3; ++p) {
                st.nhas[p] += has[i * 3 + p] > 0.0f;
                nv += has[i * 3 + p] > 0.0f;
            }
            st.nnone += nv == 0;
        }
        if (hasout) *hasout = std::move(has);
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

bool Img::CascadeDeghostingFMRaw::operator()(const input_tuple_type& intup, output_pointer& out)
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
    if (m_fmraw) {
        size_t npix = 0;
        fmidx = parse_fm(fmset, m_fm_dim, npix);
        if (!fmset) log->warn("call={} fmraw on but no FM tensor set on port 2: no cell has an active FM pixel", m_count);
        log->debug("call={} FM pixels={} (planes {}/{}/{}) dim={}", m_count, npix, fmidx.rows[0].size(), fmidx.rows[1].size(),
                   fmidx.rows[2].size(), m_fm_dim);
    }
    auto sc = use_frame ? Cascade::make_slice_charge_frame(gr, frame, m_charge_tag, m_charge_scale, m_slice_start_relative)
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
        // doc 53 part 43: at the final level (k = 1) with fmraw on, xb gains [x_0, x_1, x_2, has_0, has_1, has_2]
        std::vector<float> rawcols, rawmeans, rawhas;
        const bool with_head = m_fmraw && li + 1 == nlev;
        size_t ncol = 15;
        std::vector<float> xb18;
        if (with_head) {
            FMRawStats rs;
            const auto t0 = clock::now();
            const bool keep_means = !m_dump_dir.empty();   // the parity dumps carry e_p and has_p
            rawcols = fmraw_columns(lev, sc, fmidx, m_pca_mean, m_pca_basis, m_pca_k, rs, m_nthreads,
                                    keep_means ? &rawmeans : nullptr, keep_means ? &rawhas : nullptr);
            fmidx = FMIndex();   // the pixel table is not needed past the columns
            const size_t Nh = lev.nnodes(), NC = 3 * (size_t) m_pca_k + 3;
            ncol = 15 + NC;
            xb18.resize(Nh * ncol);
            for (size_t i = 0; i < Nh; ++i) {
                std::copy(lev.xb.begin() + i * 15, lev.xb.begin() + (i + 1) * 15, xb18.begin() + i * ncol);
                std::copy(rawcols.begin() + i * NC, rawcols.begin() + (i + 1) * NC, xb18.begin() + i * ncol + 15);
            }
            log->debug("call={} cluster={} level={} fmraw columns: nodes={} has U/V/W={}/{}/{} no view={} ncol={} {:.3f}s",
                       m_count, in->ident(), li, Nh, rs.nhas[0], rs.nhas[1], rs.nhas[2], rs.nnone, ncol,
                       std::chrono::duration<double>(clock::now() - t0).count());
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
                cnpy::npz_save(fn, "fmraw_cols", rawcols.data(), {N, 3 * (size_t) m_pca_k + 3}, "a");
                cnpy::npz_save(fn, "fm_means", rawmeans.data(), {N, (size_t) 3, (size_t) m_fm_dim}, "a");
                cnpy::npz_save(fn, "fm_has", rawhas.data(), {N, (size_t) 3}, "a");
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
