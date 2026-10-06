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
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <atomic>
#include <thread>

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
                                                     const std::string& tag, double scale, bool slice_start_relative)
{
    // the slice indexing of the activity mode (no activity read: uncer_cut = -1 skips every entry)
    SliceCharge sc = make_slice_charge(gr, scale, -1.0);
    sc.q.clear();
    sc.from_frame = true;
    sc.scale = scale;
    const double tick = frame->tick();
    // legacy: slice start taken as absolute.  MaskSlice slice starts are frame-relative (pdvd doc 122 sec 3).
    const double origin = slice_start_relative ? 0.0 : frame->time();
    for (const auto& s : sc.slice_of) {
        sc.t0.push_back((int) std::lround((s->start() - origin) / tick));
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
                                    const Universe& uni, int k, const std::string& policy, int nthreads)
{
    using clock = std::chrono::steady_clock;
    auto t0 = clock::now();
    auto t1 = t0;
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

    // ---- node features and bw edges: nodes in contiguous blocks (nthreads), each block writes its xb rows and
    // appends its own bw edges; the blocks' edges are concatenated in node order (the arrays do not depend on it)
    lev.xb.assign(N * 15, 0.0f);
    // 8 blocks per thread (dynamic, parallel_blocks) for balance; 1 block when single-threaded
    const int nfb = nthreads <= 1 ? 1 : (int) std::max<size_t>(1, std::min<size_t>((size_t) nthreads * 8, N));
    struct BwPart {
        std::vector<int64_t> src, dst;
        std::vector<float> w;
    };
    std::vector<BwPart> bwp(nfb);
    parallel_blocks(nfb, nfb, [&](size_t fb0, size_t fb1) {
      for (size_t fb = fb0; fb < fb1; ++fb) {
        auto& part = bwp[fb];
        for (size_t i = N * fb / nfb; i < N * (fb + 1) / nfb; ++i) {
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
                        part.src.push_back((int64_t) i);
                        part.dst.push_back(b);
                        part.w.push_back(1.0f);
                        last_bin = b;
                    }
                    else {
                        part.w.back() += 1.0f;
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
      }
    });
    {
        size_t ne = 0;
        for (const auto& part : bwp) ne += part.src.size();
        lev.bw_src.reserve(ne);
        lev.bw_dst.reserve(ne);
        lev.bw_w.reserve(ne);
        for (auto& part : bwp) {
            lev.bw_src.insert(lev.bw_src.end(), part.src.begin(), part.src.end());
            lev.bw_dst.insert(lev.bw_dst.end(), part.dst.begin(), part.dst.end());
            lev.bw_w.insert(lev.bw_w.end(), part.w.begin(), part.w.end());
            part = BwPart();
        }
    }

    // ---- bb: the Img::geom_clustering (BlobClustering) pairs over per-slice blob sets in time order
    t1 = clock::now();
    lev.tpart[0] = std::chrono::duration<double>(t1 - t0).count();
    {
        const auto bb = geom_pairs(blobs, lev.sidx, sc.slice_of, policy, nthreads);
        lev.bb.reserve(2 * bb.size());
        for (const auto& p : bb) {
            lev.bb.push_back(p[0]);
            lev.bb.push_back(p[1]);
        }
    }
    t0 = clock::now();
    lev.tpart[1] = std::chrono::duration<double>(t0 - t1).count();

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
    t1 = clock::now();
    lev.tpart[2] = std::chrono::duration<double>(t1 - t0).count();

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
    lev.tpart[3] = std::chrono::duration<double>(clock::now() - t1).count();
    return lev;
}

namespace {
    // Img::geom_clustering's policy table (GeomClusteringUtil.cxx): max relative slice distance and the wire
    // tolerance per relative distance
    void geom_policy(const std::string& policy, int& max_rel_diff, std::map<int, int>& gap_tol)
    {
        static const std::set<std::string> known = {"simple", "uboone", "uboone_local", "dead_clus"};
        if (!known.count(policy)) {
            THROW(ValueError() << errmsg{"CascadeGraph: geom_clustering policy \"" + policy + "\" not implemented"});
        }
        max_rel_diff = 2;
        gap_tol = {{1, 2}, {2, 1}};
        if (policy == "uboone_local") {
            gap_tol = {{1, 2}, {2, 2}};
        }
        if (policy == "simple") {
            max_rel_diff = 1;
            gap_tol = {{1, 0}};
        }
    }

    // GeomClusteringUtil.cxx rel_time_diff
    int rel_time_diff(const ISlice::pointer& one, const ISlice::pointer& two)
    {
        return std::round((two->start() - one->start()) / one->span());
    }

    // per slice index, the node indices (ascending); the non-empty ones in slice-index (time) order
    std::vector<std::vector<int64_t>> slice_sets(const std::vector<int>& sidx, size_t nslices, std::vector<int>& set_slice)
    {
        std::vector<std::vector<int64_t>> per(nslices);
        for (size_t i = 0; i < sidx.size(); ++i) per[sidx[i]].push_back((int64_t) i);
        std::vector<std::vector<int64_t>> sets;
        set_slice.clear();
        for (size_t s = 0; s < nslices; ++s) {
            if (per[s].empty()) continue;
            sets.push_back(std::move(per[s]));
            set_slice.push_back((int) s);
        }
        return sets;
    }
}  // namespace

std::vector<std::array<int64_t, 2>> Cascade::geom_pairs_graph(const std::vector<IBlob::pointer>& blobs,
                                                              const std::vector<int>& sidx,
                                                              const std::vector<ISlice::pointer>& slice_of,
                                                              const std::string& policy)
{
    const size_t N = blobs.size();
    std::vector<int> set_slice;
    const auto per = slice_sets(sidx, slice_of.size(), set_slice);
    IBlobSet::vector sets;
    for (size_t k = 0; k < per.size(); ++k) {
        IBlob::vector v;
        for (auto i : per[k]) v.push_back(blobs[i]);
        sets.push_back(std::make_shared<Aux::SimpleBlobSet>(set_slice[k], slice_of[set_slice[k]], v));
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
    return bb;
}

// Img::geom_clustering (GeomClusteringUtil.cxx L40-104) without the cluster graph.  For a "beg" blob set and each
// later set ("test") within max_rel_diff slices, RayGrid::associate(test, beg, ...) pairs a test blob a with a beg
// blob b when RayGrid::overlap finds b through the projections of the wire layers nlayers-1 .. 2 (nlayers = the
// strips of beg's first blob).  At each layer overlap selects the projected blobs occupying a grid index g with
// max(0, a.lo - tol) <= g < min(proj.size(), a.hi + tol); b occupies g for b.lo <= g < b.hi, and b.hi <= proj.size().
// So b is found iff, on every one of those layers, max(a.lo - tol, b.lo) < min(a.hi + tol, b.hi) (b.lo >= 0).
// The edge is made only between blobs of one face.  The layer nlayers-1 test goes through a projection of beg by
// grid index (as RayGrid::projection), the lower layers through the interval test.
std::vector<std::array<int64_t, 2>> Cascade::geom_pairs(const std::vector<IBlob::pointer>& blobs,
                                                        const std::vector<int>& sidx,
                                                        const std::vector<ISlice::pointer>& slice_of,
                                                        const std::string& policy, int nthreads)
{
    int max_rel_diff = 2;
    std::map<int, int> gap_tol;
    geom_policy(policy, max_rel_diff, gap_tol);
    const size_t N = blobs.size();
    std::vector<int> set_slice;
    const auto sets = slice_sets(sidx, slice_of.size(), set_slice);

    // per node: face, and the bounds of every strip by strip index (as RayGrid::overlap reads
    // blob->strips()[layer]), flat: node i's strip L is bnd[off[i] + L]
    using pairii = std::pair<int, int>;
    std::vector<const IAnodeFace*> face(N);
    std::vector<size_t> off(N + 1, 0);
    for (size_t i = 0; i < N; ++i) {
        face[i] = blobs[i]->face().get();
        off[i + 1] = off[i] + blobs[i]->shape().strips().size();
    }
    std::vector<pairii> bnd(off[N]);
    for (size_t i = 0; i < N; ++i) {
        size_t k = off[i];
        for (const auto& st : blobs[i]->shape().strips()) bnd[k++] = {st.bounds.first, st.bounds.second};
    }

    // beg sets handed out one at a time to nthreads workers (the work per set is very uneven), one output per worker,
    // each sorted unique by its worker, then merged (the pair set does not depend on which worker took which set)
    const size_t nsets = sets.size();
    const int nblk = std::max(1, std::min<int>(nthreads, (int) nsets));
    std::vector<std::vector<std::array<int64_t, 2>>> outs(nblk);
    std::atomic<size_t> next_set{0};
    auto work = [&](int blk) {
        {
            auto& out = outs[blk];
            std::vector<std::vector<int64_t>> proj;
            for (size_t ib = next_set.fetch_add(1); ib < nsets; ib = next_set.fetch_add(1)) {
                const auto& beg = sets[ib];
                const auto& sbeg = slice_of[set_slice[ib]];
                const int top = (int) (off[beg[0] + 1] - off[beg[0]]) - 1;   // nlayers - 1 of beg's first blob
                bool projected = false;
                for (size_t it = ib + 1; it < nsets; ++it) {
                    const int rel = rel_time_diff(sbeg, slice_of[set_slice[it]]);
                    if (rel > max_rel_diff) break;
                    auto gt = gap_tol.find(rel);
                    if (gt == gap_tol.end()) continue;
                    const int tol = gt->second;
                    if (top < 2) continue;   // RayGrid::overlap never reaches a wire layer
                    if (!projected) {        // RayGrid::projection(references(beg), top)
                        proj.clear();
                        for (auto b : beg) {
                            const auto& bd = bnd[off[b] + top];
                            for (int g = bd.first; g < bd.second; ++g) {
                                if ((int) proj.size() <= g) proj.resize(g + 1);
                                proj[g].push_back(b);
                            }
                        }
                        projected = true;
                    }
                    for (auto a : sets[it]) {
                        const pairii* ab = &bnd[off[a]];
                        const int lo = std::max(0, ab[top].first - tol);
                        const int hi = std::min((int) proj.size(), ab[top].second + tol);
                        for (int g = lo; g < hi; ++g) {
                            for (auto b : proj[g]) {
                                const pairii* bb = &bnd[off[b]];
                                // b sits in the buckets [b.lo, b.hi): take it once, at its first bucket in [lo, hi)
                                if (g != std::max(lo, bb[top].first)) continue;
                                if (face[a] != face[b]) continue;
                                bool ok = true;
                                for (int L = top - 1; L >= 2 && ok; --L) {
                                    ok = std::max(ab[L].first - tol, bb[L].first) < std::min(ab[L].second + tol, bb[L].second);
                                }
                                if (ok) out.push_back({std::min(a, b), std::max(a, b)});
                            }
                        }
                    }
                }
            }
            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
        }
    };
    parallel_workers(nblk, work);
    std::vector<std::array<int64_t, 2>> out;
    size_t tot = 0;
    for (const auto& o : outs) tot += o.size();
    out.reserve(tot);
    for (auto& o : outs) {
        const size_t mid = out.size();
        out.insert(out.end(), o.begin(), o.end());
        std::vector<std::array<int64_t, 2>>().swap(o);
        std::inplace_merge(out.begin(), out.begin() + mid, out.end());
    }
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
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

size_t Cascade::final_guard(const Level& lev, const std::vector<float>& logit, std::vector<bool>& keep)
{
    const size_t N = lev.nnodes();
    const size_t NW = lev.wq.size();
    const size_t E = lev.bw_src.size();
    // per wire node, its nodes (CSR over bw_dst, nodes in edge order = node order)
    std::vector<size_t> ptr(NW + 1, 0);
    for (size_t e = 0; e < E; ++e) ptr[lev.bw_dst[e] + 1] += 1;
    for (size_t w = 0; w < NW; ++w) ptr[w + 1] += ptr[w];
    std::vector<int64_t> nodes(E);
    {
        std::vector<size_t> fill(ptr.begin(), ptr.end() - 1);
        for (size_t e = 0; e < E; ++e) nodes[fill[lev.bw_dst[e]]++] = lev.bw_src[e];
    }
    // per node, its wire nodes (bw_src is non-decreasing by construction)
    std::vector<size_t> nptr(N + 1, 0);
    for (size_t e = 0; e < E; ++e) nptr[lev.bw_src[e] + 1] += 1;
    for (size_t i = 0; i < N; ++i) nptr[i + 1] += nptr[i];
    std::vector<int> cov(NW, 0);
    for (size_t e = 0; e < E; ++e) {
        if (keep[lev.bw_src[e]]) cov[lev.bw_dst[e]] += 1;
    }
    size_t nadded = 0;
    for (size_t w = 0; w < NW; ++w) {
        if (!(lev.wq[w] > 0) || cov[w] > 0 || ptr[w] == ptr[w + 1]) continue;
        int64_t best = -1;
        for (size_t k = ptr[w]; k < ptr[w + 1]; ++k) {
            const int64_t i = nodes[k];
            if (best < 0 || logit[i] > logit[best] || (logit[i] == logit[best] && i < best)) best = i;
        }
        keep[best] = true;
        ++nadded;
        for (size_t e = nptr[best]; e < nptr[best + 1]; ++e) cov[lev.bw_dst[e]] += 1;
    }
    return nadded;
}

void Cascade::parallel_workers(int nthreads, const std::function<void(int)>& f)
{
    const int nw = std::max(1, nthreads);
    if (nw == 1) {
        f(0);
        return;
    }
    std::vector<std::thread> th;
    std::vector<std::exception_ptr> err(nw);
    for (int w = 0; w < nw; ++w) {
        th.emplace_back([&, w]() {
            try {
                f(w);
            }
            catch (...) {
                err[w] = std::current_exception();
            }
        });
    }
    for (auto& t : th) t.join();
    for (auto& e : err) {
        if (e) std::rethrow_exception(e);
    }
}

void Cascade::parallel_blocks(size_t n, int nthreads, const std::function<void(size_t, size_t)>& f)
{
    const size_t nw = std::min<size_t>(n, (size_t) std::max(1, nthreads));
    if (nw <= 1) {
        f(0, n);
        return;
    }
    // dynamic: ~32 chunks per thread handed out by an atomic counter (the work per item is uneven: an event's
    // activity sits in a narrow range of slices); callers write per-item or per-block slots, so the result does not
    // depend on which thread ran which chunk
    const size_t chunk = std::max<size_t>(1, n / (nw * 32));
    std::atomic<size_t> next{0};
    parallel_workers((int) nw, [&](int) {
        for (size_t c0 = next.fetch_add(chunk); c0 < n; c0 = next.fetch_add(chunk)) {
            f(c0, std::min(n, c0 + chunk));
        }
    });
}

Cascade::IsoResult Cascade::iso_fallback(const std::vector<int>& face, const std::vector<int>& sidx,
                                         const std::vector<float>& wq, const std::vector<int>& wsidx,
                                         const std::vector<float>& logit, std::vector<bool>& keep,
                                         const IsoParams& par)
{
    IsoResult res;
    const size_t N = logit.size();
    if (face.size() != N || sidx.size() != N || keep.size() != N || wsidx.size() != wq.size()) {
        THROW(ValueError() << errmsg{"Cascade::iso_fallback: array sizes differ"});
    }
    // charged wire nodes per slice index
    std::map<int, size_t> nch;
    for (size_t w = 0; w < wq.size(); ++w) {
        if (wq[w] > 0) ++nch[wsidx[w]];
    }
    // (face, slice) groups: node count and ambiguous count
    std::map<std::pair<int, int>, std::pair<size_t, size_t>> grp;
    for (size_t i = 0; i < N; ++i) {
        auto& g = grp[{face[i], sidx[i]}];
        ++g.first;
        const double l = logit[i];
        if (l > par.amb_lo && l < par.amb_hi) ++g.second;
    }
    std::map<std::pair<int, int>, bool> trig;
    for (const auto& [key, g] : grp) {
        auto it = nch.find(key.second);
        const double nw = std::max<double>(1.0, it == nch.end() ? 0.0 : (double) it->second);
        const double n = (double) g.first;
        const bool t = g.first >= (size_t) std::max(0, par.nmin) && n >= par.mmin * nw && (double) g.second >= par.amin * n;
        trig[key] = t;
        if (t) ++res.nslices;
    }
    for (size_t i = 0; i < N; ++i) {
        if (keep[i] || (double) logit[i] < par.t_keep) continue;
        if (trig[{face[i], sidx[i]}]) {
            keep[i] = true;
            ++res.nadded;
        }
    }
    return res;
}
