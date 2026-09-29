// Level-graph builder of CascadeDeghosting (wcfm doc 14).  See WireCellImg/CascadeGraph.h.
// Python references (wcp-porting-img wcfm/scripts): gnn_dataset.py build_graph L152-276, inslice_adjacency L72-99,
// wire_adjacency L102-131; d12_cascade.py contract_level L223-299, guard_prune L340-357.

#include "WireCellImg/CascadeGraph.h"
#include "WireCellImg/GeomClusteringUtil.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

using namespace WireCell;
using namespace WireCell::Img;

int Cascade::SliceCharge::index(const ISlice::pointer& s) const
{
    auto it = sidx.find(s.get());
    if (it == sidx.end()) {
        THROW(ValueError() << errmsg{"CascadeGraph: blob slice not in the cluster graph"});
    }
    return it->second;
}

float Cascade::SliceCharge::charge(int s, int ch) const
{
    if (from_frame) {
        auto it = traces.find(ch);
        if (it == traces.end()) return 0.0f;
        const auto& tr = it->second;
        const auto& qs = tr->charge();
        const int tb = tr->tbin();
        double sum = 0;
        for (int t = t0[s]; t < t0[s] + nt[s]; ++t) {
            const int i = t - tb;
            if (i >= 0 && i < (int) qs.size()) sum += qs[i];
        }
        return (float) (sum * scale);
    }
    const auto& m = q[s];
    auto it = m.find(ch);
    return it == m.end() ? 0.0f : it->second;
}

Cascade::SliceCharge Cascade::make_slice_charge(const cluster_graph_t& gr, double scale, double uncer_cut)
{
    SliceCharge sc;
    std::vector<ISlice::pointer> slices;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& node = gr[vtx];
        if (node.code() != 's') continue;
        slices.push_back(std::get<ISlice::pointer>(node.ptr));
    }
    // by start time; vertex order breaks ties (stable)
    std::stable_sort(slices.begin(), slices.end(),
                     [](const ISlice::pointer& a, const ISlice::pointer& b) { return a->start() < b->start(); });
    for (const auto& s : slices) {
        if (sc.slice_of.empty() || s->start() != sc.slice_of.back()->start()) {
            sc.slice_of.push_back(s);
            sc.q.emplace_back();
        }
        const int si = (int) sc.slice_of.size() - 1;
        sc.sidx[s.get()] = si;
        auto& m = sc.q[si];
        // The activity map is keyed by channel pointer; only an order-independent reduction (max) is taken.
        for (const auto& [ich, val] : s->activity()) {
            if (val.uncertainty() >= uncer_cut) continue;   // dummy / masked plane entries
            const float v = (float) (val.value() * scale);
            auto [it, fresh] = m.emplace(ich->ident(), v);
            if (!fresh) it->second = std::max(it->second, v);
        }
    }
    return sc;
}

Cascade::SliceCharge Cascade::make_slice_charge_frame(const cluster_graph_t& gr, const IFrame::pointer& frame,
                                                     const std::string& tag, double scale)
{
    // the slice indexing of the activity mode (no activity read: uncer_cut = -1 skips every entry)
    SliceCharge sc = make_slice_charge(gr, scale, -1.0);
    sc.q.clear();
    sc.from_frame = true;
    sc.scale = scale;
    const double tick = frame->tick();
    for (const auto& s : sc.slice_of) {
        sc.t0.push_back((int) std::lround((s->start() - frame->time()) / tick));
        sc.nt.push_back((int) std::lround(s->span() / tick));
    }
    const auto& trs = tag.empty() ? *frame->traces() : ITrace::vector();
    ITrace::vector tagged;
    if (!tag.empty()) {
        auto all = frame->traces();
        for (size_t ind : frame->tagged_traces(tag)) tagged.push_back(all->at(ind));
    }
    for (const auto& tr : tag.empty() ? trs : tagged) {
        sc.traces.emplace(tr->channel(), tr);   // first trace of a channel wins (one per channel expected)
    }
    return sc;
}

namespace {
    // (plane, lo, hi) of the three wire-plane strips (RayGrid layers 2, 3, 4), half-open wire-in-plane ranges
    std::array<int, 6> wire_ranges(const IBlob::pointer& blob)
    {
        std::array<int, 6> r{0, 0, 0, 0, 0, 0};
        for (const auto& strip : blob->shape().strips()) {
            const int ip = strip.layer - 2;  // BlobClustering add_blobs: two non-plane layers first
            if (ip < 0 || ip > 2) continue;
            r[2 * ip] = strip.bounds.first;
            r[2 * ip + 1] = strip.bounds.second;
        }
        return r;
    }
}  // namespace

std::array<std::vector<int>, 3> Cascade::blob_channels(const IBlob::pointer& blob)
{
    std::array<std::vector<int>, 3> out;
    auto planes = blob->face()->planes();
    for (const auto& strip : blob->shape().strips()) {
        const int ip = strip.layer - 2;
        if (ip < 0 || ip > 2) continue;
        const auto& wires = planes[ip]->wires();
        auto& v = out[ip];
        for (int wip = strip.bounds.first; wip < strip.bounds.second && wip < (int) wires.size(); ++wip) {
            v.push_back(wires[wip]->channel());
        }
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return out;
}

int Cascade::Universe::index(int plane, int s, int ch) const
{
    auto it = lookup.find(key(plane, s, ch));
    return it == lookup.end() ? -1 : it->second;
}

Cascade::Universe Cascade::make_universe(const std::vector<IBlob::pointer>& blobs, const SliceCharge& sc)
{
    Universe uni;
    // coverage per (face ident, slice index), per plane: difference array over wire-in-plane indices
    struct Group {
        IAnodeFace::pointer face;
        std::array<std::vector<std::pair<int, int>>, 3> iv;
    };
    std::map<int64_t, Group> groups;
    std::vector<std::array<int, 3>> keys;
    for (const auto& b : blobs) {
        const int s = sc.index(b->slice());
        const auto ch = blob_channels(b);
        for (int p = 0; p < 3; ++p) {
            for (int c : ch[p]) keys.push_back({p, s, c});
        }
        const auto r = wire_ranges(b);
        auto& g = groups[(int64_t) b->face()->ident() * 10000000LL + s];
        g.face = b->face();
        for (int p = 0; p < 3; ++p) g.iv[p].emplace_back(r[2 * p], r[2 * p + 1]);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    uni.wires = keys;
    uni.q.resize(keys.size());
    uni.lookup.reserve(keys.size());
    for (size_t i = 0; i < keys.size(); ++i) {
        uni.lookup[Universe::key(keys[i][0], keys[i][1], keys[i][2])] = (int) i;
        uni.q[i] = sc.charge(keys[i][1], keys[i][2]);
    }
    // gnn_dataset.wire_adjacency: consecutive covered wips (w, w+1) of one (face, slice, plane) on two
    // different channels give an edge between their wire nodes; deduplicated over faces
    std::vector<std::array<int64_t, 2>> ww;
    for (const auto& [gk, g] : groups) {
        const int s = (int) (gk % 10000000LL);
        auto planes = g.face->planes();
        for (int p = 0; p < 3; ++p) {
            const auto& iv = g.iv[p];
            int wmin = iv[0].first, wmax = iv[0].second;
            for (const auto& [lo, hi] : iv) {
                wmin = std::min(wmin, lo);
                wmax = std::max(wmax, hi);
            }
            if (wmax - wmin < 2) continue;
            std::vector<int> cov(wmax - wmin + 1, 0);
            for (const auto& [lo, hi] : iv) {
                cov[lo - wmin] += 1;
                cov[hi - wmin] -= 1;
            }
            const auto& wires = planes[p]->wires();
            int run = 0;
            std::vector<char> covered(wmax - wmin, 0);
            for (int i = 0; i < wmax - wmin; ++i) {
                run += cov[i];
                covered[i] = run > 0;
            }
            for (int i = 0; i + 1 < wmax - wmin; ++i) {
                if (!(covered[i] && covered[i + 1])) continue;
                const int w = wmin + i;
                if (w + 1 >= (int) wires.size()) continue;
                const int c1 = wires[w]->channel(), c2 = wires[w + 1]->channel();
                if (c1 == c2) continue;
                const int i1 = uni.index(p, s, c1), i2 = uni.index(p, s, c2);
                if (i1 < 0 || i2 < 0) continue;
                ww.push_back({std::min(i1, i2), std::max(i1, i2)});
            }
        }
    }
    std::sort(ww.begin(), ww.end());
    ww.erase(std::unique(ww.begin(), ww.end()), ww.end());
    uni.ww = std::move(ww);
    return uni;
}

std::vector<std::array<int64_t, 2>> Cascade::inslice_pairs(const std::vector<int64_t>& group,
                                                           const std::vector<std::array<int, 6>>& rng)
{
    const size_t n = group.size();
    std::vector<int64_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    // np.lexsort((u0, keys)): by group, then u0, stable
    std::stable_sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
        if (group[a] != group[b]) return group[a] < group[b];
        return rng[a][0] < rng[b][0];
    });
    std::vector<std::array<int64_t, 2>> out;
    size_t g0 = 0;
    while (g0 < n) {
        size_t g1 = g0 + 1;
        while (g1 < n && group[order[g1]] == group[order[g0]]) ++g1;
        for (size_t a = g0; a + 1 < g1; ++a) {
            const auto& ra = rng[order[a]];
            for (size_t j = a + 1; j < g1; ++j) {
                const auto& rj = rng[order[j]];
                if (rj[0] > ra[1]) break;           // u0_j > u1_a: no later candidate (sorted by u0)
                if (std::max(ra[2], rj[2]) <= std::min(ra[3], rj[3]) &&
                    std::max(ra[4], rj[4]) <= std::min(ra[5], rj[5])) {
                    const int64_t x = order[a], y = order[j];
                    out.push_back({std::min(x, y), std::max(x, y)});
                }
            }
        }
        g0 = g1;
    }
    std::sort(out.begin(), out.end());
    return out;
}

Cascade::Level Cascade::build_level(const std::vector<IBlob::pointer>& blobs, const SliceCharge& sc,
                                    const Universe& uni, int k, const std::string& policy)
{
    Level lev;
    lev.k = std::max(k, 1);
    lev.blobs = blobs;
    const size_t N = blobs.size();
    lev.sidx.resize(N);
    for (size_t i = 0; i < N; ++i) lev.sidx[i] = sc.index(blobs[i]->slice());

    // ---- (super-)wire nodes: every universe wire, binned (plane, slice, channel // k)
    const size_t NU = uni.wires.size();
    std::vector<int> ubin(NU);
    size_t NB = 0;
    if (lev.k == 1) {
        std::iota(ubin.begin(), ubin.end(), 0);
        NB = NU;
    }
    else {
        std::vector<std::array<int, 3>> bk(NU);
        for (size_t i = 0; i < NU; ++i) bk[i] = {uni.wires[i][0], uni.wires[i][1], uni.wires[i][2] / lev.k};
        std::vector<std::array<int, 3>> ub = bk;
        std::sort(ub.begin(), ub.end());
        ub.erase(std::unique(ub.begin(), ub.end()), ub.end());
        NB = ub.size();
        for (size_t i = 0; i < NU; ++i) {
            ubin[i] = (int) (std::lower_bound(ub.begin(), ub.end(), bk[i]) - ub.begin());
        }
    }
    std::vector<double> wq(NB, 0.0);
    lev.wplane.assign(NB, 0);
    lev.wchan.assign(NB, 0);
    lev.wsidx.assign(NB, 0);
    for (size_t i = 0; i < NU; ++i) {
        wq[ubin[i]] += uni.q[i];
        lev.wplane[ubin[i]] = uni.wires[i][0];
        lev.wsidx[ubin[i]] = uni.wires[i][1];
        lev.wchan[ubin[i]] = uni.wires[i][2] / lev.k * lev.k;
    }
    lev.wq.resize(NB);
    for (size_t b = 0; b < NB; ++b) lev.wq[b] = (float) wq[b];

    // ---- node features and bw edges
    lev.xb.assign(N * 15, 0.0f);
    for (size_t i = 0; i < N; ++i) {
        const auto ch = blob_channels(blobs[i]);
        const int s = lev.sidx[i];
        double sums[3] = {0, 0, 0};
        for (int p = 0; p < 3; ++p) {
            double sq = 0;
            int npx = 0;
            int last_bin = -1;
            for (int c : ch[p]) {
                const int u = uni.index(p, s, c);
                if (u < 0) {
                    THROW(ValueError() << errmsg{"CascadeGraph: a blob channel is outside the wire universe"});
                }
                const float q = uni.q[u];
                sq += q;
                npx += (q > 0);
                const int b = ubin[u];
                if (b != last_bin) {
                    lev.bw_src.push_back((int64_t) i);
                    lev.bw_dst.push_back(b);
                    lev.bw_w.push_back(1.0f);
                    last_bin = b;
                }
                else {
                    lev.bw_w.back() += 1.0f;
                }
            }
            const size_t nch = ch[p].size();
            float* x = &lev.xb[i * 15 + 4 * p];
            x[0] = (float) nch;
            x[1] = (float) npx;
            x[2] = (float) std::log1p(sq);
            x[3] = nch > 0 ? (float) std::log1p(sq / nch) : 0.0f;
            sums[p] = std::log1p(sq);
        }
        lev.xb[i * 15 + 12] = (float) (sums[0] - sums[1]);
        lev.xb[i * 15 + 13] = (float) (sums[0] - sums[2]);
        lev.xb[i * 15 + 14] = (float) (sums[1] - sums[2]);
    }

    // ---- bb: Img::geom_clustering (the BlobClustering rule) over per-slice blob sets in time order
    {
        const size_t NS = sc.slice_of.size();
        std::vector<IBlob::vector> per(NS);
        for (size_t i = 0; i < N; ++i) per[lev.sidx[i]].push_back(blobs[i]);
        IBlobSet::vector sets;
        for (size_t s = 0; s < NS; ++s) {
            if (per[s].empty()) continue;
            sets.push_back(std::make_shared<Aux::SimpleBlobSet>((int) s, sc.slice_of[s], per[s]));
        }
        cluster_indexed_graph_t grind;
        for (auto it = sets.begin(); it != sets.end(); ++it) {
            Img::geom_clustering(grind, it, sets.end(), policy);
        }
        std::unordered_map<const IBlob*, int64_t> bidx;  // lookup only
        bidx.reserve(N);
        for (size_t i = 0; i < N; ++i) bidx[blobs[i].get()] = (int64_t) i;
        std::vector<std::array<int64_t, 2>> bb;
        const auto& g = grind.graph();
        for (auto e : boost::make_iterator_range(boost::edges(g))) {
            const auto& na = g[boost::source(e, g)];
            const auto& nb = g[boost::target(e, g)];
            if (na.code() != 'b' || nb.code() != 'b') continue;
            const int64_t a = bidx.at(std::get<IBlob::pointer>(na.ptr).get());
            const int64_t b = bidx.at(std::get<IBlob::pointer>(nb.ptr).get());
            if (a == b) continue;
            bb.push_back({std::min(a, b), std::max(a, b)});
        }
        std::sort(bb.begin(), bb.end());
        bb.erase(std::unique(bb.begin(), bb.end()), bb.end());
        for (const auto& p : bb) {
            lev.bb.push_back(p[0]);
            lev.bb.push_back(p[1]);
        }
    }

    // ---- bb_in: same face and slice, overlap or abut on all three planes
    {
        std::vector<int64_t> group(N);
        std::vector<std::array<int, 6>> rng(N);
        for (size_t i = 0; i < N; ++i) {
            group[i] = (int64_t) blobs[i]->face()->ident() * 10000000LL + lev.sidx[i];
            rng[i] = wire_ranges(blobs[i]);
        }
        for (const auto& p : inslice_pairs(group, rng)) {
            lev.bb_in.push_back(p[0]);
            lev.bb_in.push_back(p[1]);
        }
    }

    // ---- ww: the universe's wire-wire pairs, binned
    {
        std::vector<std::array<int64_t, 2>> ww;
        ww.reserve(uni.ww.size());
        for (const auto& p : uni.ww) {
            const int64_t a = ubin[p[0]], b = ubin[p[1]];
            if (a == b) continue;
            ww.push_back({std::min(a, b), std::max(a, b)});
        }
        std::sort(ww.begin(), ww.end());
        ww.erase(std::unique(ww.begin(), ww.end()), ww.end());
        for (const auto& p : ww) {
            lev.ww.push_back(p[0]);
            lev.ww.push_back(p[1]);
        }
    }
    return lev;
}

std::vector<bool> Cascade::guard_prune(const Level& lev, const std::vector<int>& cand_order, size_t& nguarded)
{
    const size_t N = lev.nnodes();
    const size_t E = lev.bw_src.size();
    // per node, its bw edges (bw_src is non-decreasing by construction)
    std::vector<size_t> ptr(N + 1, 0);
    for (size_t e = 0; e < E; ++e) ptr[lev.bw_src[e] + 1] += 1;
    for (size_t i = 0; i < N; ++i) ptr[i + 1] += ptr[i];
    std::vector<int> cov(lev.wq.size(), 0);
    for (size_t e = 0; e < E; ++e) cov[lev.bw_dst[e]] += 1;
    std::vector<bool> pruned(N, false);
    nguarded = 0;
    for (int b : cand_order) {
        bool guard = false;
        for (size_t e = ptr[b]; e < ptr[b + 1]; ++e) {
            if (cov[lev.bw_dst[e]] <= 1) {
                guard = true;
                break;
            }
        }
        if (guard) {
            ++nguarded;
            continue;
        }
        for (size_t e = ptr[b]; e < ptr[b + 1]; ++e) cov[lev.bw_dst[e]] -= 1;
        pruned[b] = true;
    }
    return pruned;
}
