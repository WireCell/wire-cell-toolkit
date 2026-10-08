// The cascade level loop of CascadeDeghosting (wcfm doc 14), extracted (pdvd doc 130 phase 2).  See
// WireCellAux/CascadeRun.h.  The body is CascadeDeghosting.cxx operator() at toolkit 02d90c1a, lines 278-487, with
// the node's members read from RunParams / LevelSpec; the messages are unchanged.

#include "WireCellAux/CascadeRun.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/MemUsage.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/cnpy.h"

#include <algorithm>
#include <chrono>
#include <cmath>

using namespace WireCell;
using namespace WireCell::Aux;

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

    // (logit, qhat) of every node of the level
    void run_forward(const ITensorForward::pointer& fwd, const Level& lev, int ident,
                     std::vector<float>& logit, std::vector<float>& qhat)
    {
        const size_t N = lev.nnodes(), M = lev.wq.size(), E = lev.bw_src.size();
        ITensor::vector tv{
            tens(lev.xb, {N, 15}),        tens(lev.wq, {M}),        tens(lev.wplane, {M}),
            tens(lev.bw_src, {E}),        tens(lev.bw_dst, {E}),    tens(lev.bw_w, {E}),
            tens(lev.bb, {lev.bb.size() / 2, 2}), tens(lev.bb_in, {lev.bb_in.size() / 2, 2}),
            tens(lev.ww, {lev.ww.size() / 2, 2})};
        auto in = std::make_shared<Aux::SimpleTensorSet>(ident, Configuration(), std::make_shared<ITensor::vector>(tv));
        auto out = fwd->forward(in);
        if (!out) {
            THROW(RuntimeError() << errmsg{"CascadeDeghosting: forward failed"});
        }
        auto ot = out->tensors();
        if (ot->size() < 2) {
            THROW(RuntimeError() << errmsg{"CascadeDeghosting: forward must return (logit, qhat)"});
        }
        for (int k = 0; k < 2; ++k) {
            const auto& t = ot->at(k);
            if (t->dtype() != "f4" || t->size() != N * sizeof(float)) {
                THROW(RuntimeError() << errmsg{String::format("CascadeDeghosting: output %d is %s of %d bytes, want f4 [%d]",
                                                              k, t->dtype(), (int) t->size(), (int) N)});
            }
            const float* p = reinterpret_cast<const float*>(t->data());
            (k == 0 ? logit : qhat).assign(p, p + N);
        }
    }

    std::vector<std::array<int64_t, 2>> pairs_of(const std::vector<int64_t>& flat)
    {
        std::vector<std::array<int64_t, 2>> out(flat.size() / 2);
        for (size_t i = 0; i < out.size(); ++i) out[i] = {flat[2 * i], flat[2 * i + 1]};
        return out;
    }

}  // namespace

double Cascade::default_repair_budget(double t)
{
    return 3.0 * (-std::log(0.2)) + std::log1p(std::exp(-t));  // -log sigmoid(t)
}

Cascade::RunResult Cascade::run_cascade(const std::vector<IBlob::pointer>& blobs, const SliceCharge& sc,
                                        const std::vector<LevelSpec>& levels, const RunParams& rp, int ident,
                                        Log::logptr_t log, int call)
{
    using clock = std::chrono::steady_clock;
    std::vector<IBlob::pointer> cur = blobs;
    const size_t nin = cur.size();
    auto uni = Cascade::make_universe(cur, sc);
    std::vector<int> parent(cur.size(), -1);   // index in the previous level
    std::vector<int> depth(cur.size(), 0);     // bisections from the uncut blob
    const Cascade::CutParams& cpar = rp.cut;
    int next_ident = rp.ident_base;
    std::vector<IBlob::pointer> kept;
    Cascade::SteinerResult srep;
    size_t nkeep_thr = 0;

    const size_t nlev = levels.size();
    for (size_t li = 0; li < nlev; ++li) {
        const auto& lc = levels[li];
        const auto tl = clock::now();
        if (cur.size() < 2) {   // d12_cascade: a graph with < 2 nodes is empty, nothing survives
            log->debug("call={} cluster={} level={} nodes={}: < 2, nothing kept", call, ident, li, cur.size());
            cur.clear();
            break;
        }
        auto lev = Cascade::build_level(cur, sc, uni, lc.superwire, rp.policy, rp.nthreads);
        const auto tb = clock::now();
        const double rss_build = memusage_resident();
        std::vector<float> logit, qhat;
        run_forward(lc.forward, lev, ident, logit, qhat);
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
            if (rp.guard) {
                pruned = Cascade::guard_prune(lev, cand, nguard);
            }
            else {
                for (int b : cand) pruned[b] = true;
            }
            for (size_t i = 0; i < N; ++i) decision[i] = !pruned[i];
            // cut the survivors to the next width, one parent at a time (exact parent map)
            std::vector<int> order;
            for (size_t i = 0; i < N; ++i) if (decision[i]) order.push_back((int) i);
            if (rp.max_level_nodes > 0) {
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return logit[a] > logit[b]; });
            }
            struct Kid {
                int parent;
                IBlob::vector blobs;
                std::vector<int> depth;
            };
            std::vector<Kid> kids;
            size_t nnext = 0;
            const int width = levels[li + 1].width;
            // the bisection of every survivor that needs it, first (a pure function of the shape; nthreads blocks),
            // then consumed in `order` exactly as before (doc 15 round 3)
            std::vector<std::vector<std::pair<RayGrid::Blob, int>>> cut(order.size());
            Cascade::parallel_blocks(order.size(), rp.nthreads, [&](size_t k0, size_t k1) {
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
                if (rp.max_level_nodes > 0 && nnext + k.blobs.size() > (size_t) rp.max_level_nodes) {
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
            if (rp.repair) {
                auto edges = pairs_of(lev.bb);
                auto ein = pairs_of(lev.bb_in);
                edges.insert(edges.end(), ein.begin(), ein.end());
                const Cascade::SteinerParams& par = rp.steiner;
                if (rp.edge_length) {   // pdvd doc 130: gap cells priced by length
                    std::vector<double> elen(edges.size());
                    for (size_t k = 0; k < edges.size(); ++k) elen[k] = rp.edge_length(cur[edges[k][0]], cur[edges[k][1]]);
                    srep = Cascade::steiner_repair(N, edges, logit, qhat, keep, par, &elen);
                }
                else {
                    srep = Cascade::steiner_repair(N, edges, logit, qhat, keep, par);
                }
                keep = srep.keep;
            }
            if (rp.iso) {   // wcfm doc 17: dense ambiguous slices keep their cells down to t_keep
                std::vector<int> fc(N);
                for (size_t i = 0; i < N; ++i) fc[i] = cur[i]->face()->which();
                const Cascade::IsoParams& ip = rp.iso_params;
                const auto ir = Cascade::iso_fallback(fc, lev.sidx, lev.wq, lev.wsidx, logit, keep, ip);
                log->debug("call={} cluster={} iso_fallback slices={} added={}", call, ident, ir.nslices,
                           ir.nadded);
            }
            if (rp.final_guard) {   // wcfm doc 22: every charged wire node keeps at least one 3-D explanation
                const size_t ng = Cascade::final_guard(lev, logit, keep);
                log->debug("call={} cluster={} final_guard added={}", call, ident, ng);
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
                   call, ident, li, lev.k, N, lev.wq.size(), lev.bw_src.size(), lev.bb.size() / 2,
                   lev.bb_in.size() / 2, lev.ww.size() / 2, ncand, nguard, ncap, nsurv, next.size(), dt);
        // wcfm doc 15: wall per phase (build = wires+features / bb / bb_in / ww), resident memory after each (MB)
        log->debug("call={} cluster={} level={} phases build={:.2f}s ({:.2f} {:.2f} {:.2f} {:.2f}) forward={:.2f}s "
                   "select={:.2f}s rss after build={:.0f} forward={:.0f} select={:.0f} MB",
                   call, ident, li, std::chrono::duration<double>(tb - tl).count(), lev.tpart[0],
                   lev.tpart[1], lev.tpart[2], lev.tpart[3], std::chrono::duration<double>(tf - tb).count(),
                   std::chrono::duration<double>(te - tf).count(), rss_build / 1024, rss_fwd / 1024,
                   memusage_resident() / 1024);

        if (!rp.dump_dir.empty()) {
            const std::string fn = String::format("%s/cascade-%d-L%d.npz", rp.dump_dir.c_str(), ident, (int) li);
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
            cnpy::npz_save(fn, "xb", lev.xb.data(), {N, 15}, "w");
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

    RunResult res;
    res.kept = std::move(kept);
    res.nin = nin;
    res.nkeep_thr = nkeep_thr;
    res.srep = std::move(srep);
    return res;
}
