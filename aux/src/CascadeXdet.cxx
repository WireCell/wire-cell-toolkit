// The G + D columns of the pooled cascade model from a file (wcfm doc 53 W7).  See WireCellAux/CascadeXdet.h.

#include "WireCellAux/CascadeXdet.h"

#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/ISlice.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/cnpy.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Aux;

namespace {

    const cnpy::NpyArray& get_array(const cnpy::npz_t& npz, const std::string& path, const std::string& name,
                                    size_t rows, size_t cols, size_t level)
    {
        auto it = npz.find(name);
        if (it == npz.end()) {
            THROW(ValueError() << errmsg{String::format("xdet columns %s: level %d has no %s", path, (int) level, name)});
        }
        const auto& a = it->second;
        const bool shape_ok = cols == 0 ? (a.shape.size() == 1 && a.shape[0] == rows)
                                        : (a.shape.size() == 2 && a.shape[0] == rows && a.shape[1] == cols);
        if (a.word_size != 4 || a.fortran_order || !shape_ok) {
            THROW(ValueError() << errmsg{String::format("xdet columns %s: %s must be a C-order 4-byte array [%d%s]",
                                                        path, name, (int) rows,
                                                        cols ? String::format(", %d", (int) cols) : std::string(""))});
        }
        return a;
    }

    template <typename K, typename V>
    void insert_unique(std::map<K, V>& m, const K& k, const V& v, const std::string& path, const std::string& what,
                       size_t level)
    {
        auto [it, fresh] = m.emplace(k, v);
        if (!fresh && it->second != v) {
            THROW(ValueError() << errmsg{String::format("xdet columns %s: level %d %s key given twice with different values",
                                                        path, (int) level, what)});
        }
    }

}  // namespace

Cascade::XdetColumns Cascade::load_xdet_columns(const std::string& path)
{
    XdetColumns xc;
    xc.path = path;
    const cnpy::npz_t npz = cnpy::npz_load(path);
    for (size_t li = 0;; ++li) {
        const std::string pre = String::format("L%d_", (int) li);
        auto ib = npz.find(pre + "bkey");
        if (ib == npz.end()) break;
        if (ib->second.shape.size() != 2) {
            THROW(ValueError() << errmsg{String::format("xdet columns %s: %sbkey must be 2-D", path, pre)});
        }
        const size_t nb = ib->second.shape[0];
        const auto& bkey = get_array(npz, path, pre + "bkey", nb, 8, li);
        const auto& bcols = get_array(npz, path, pre + "bcols", nb, xdet_blob_cols, li);
        auto iw = npz.find(pre + "wkey");
        if (iw == npz.end() || iw->second.shape.size() != 2) {
            THROW(ValueError() << errmsg{String::format("xdet columns %s: level %d needs a 2-D %swkey", path, (int) li, pre)});
        }
        const size_t nw = iw->second.shape[0];
        const auto& wkey = get_array(npz, path, pre + "wkey", nw, 2, li);
        const auto& wcols = get_array(npz, path, pre + "wcols", nw, xdet_wire_cols, li);
        const auto& live = get_array(npz, path, pre + "live", nw, 0, li);

        XdetLevelTable tab;
        const int32_t* pk = nb ? bkey.data<int32_t>() : nullptr;
        const float* pc = nb ? bcols.data<float>() : nullptr;
        for (size_t i = 0; i < nb; ++i) {
            XdetBlobKey k;
            std::array<float, xdet_blob_cols> v;
            std::copy(pk + 8 * i, pk + 8 * (i + 1), k.begin());
            std::copy(pc + xdet_blob_cols * i, pc + xdet_blob_cols * (i + 1), v.begin());
            insert_unique(tab.blob, k, v, path, "blob", li);
        }
        const int32_t* qk = nw ? wkey.data<int32_t>() : nullptr;
        const float* qc = nw ? wcols.data<float>() : nullptr;
        const float* ql = nw ? live.data<float>() : nullptr;
        for (size_t m = 0; m < nw; ++m) {
            XdetWireKey k{qk[2 * m], qk[2 * m + 1]};
            std::array<float, xdet_wire_cols + 1> v;
            std::copy(qc + xdet_wire_cols * m, qc + xdet_wire_cols * (m + 1), v.begin());
            v[xdet_wire_cols] = ql[m];
            insert_unique(tab.wire, k, v, path, "wire", li);
        }
        xc.levels.push_back(std::move(tab));
    }
    if (xc.levels.empty()) {
        THROW(ValueError() << errmsg{String::format("xdet columns %s: no level (L0_bkey missing)", path)});
    }
    return xc;
}

Cascade::XdetBlobKey Cascade::xdet_blob_key(const IBlob::pointer& blob)
{
    // the dump's wip row (CascadeRun.cxx): strip layers 2-4 are U, V, W; a missing layer stays 0
    XdetBlobKey k{};
    k[0] = blob->face()->which();
    k[1] = blob->slice()->ident();
    for (const auto& strip : blob->shape().strips()) {
        const int ip = strip.layer - 2;
        if (ip < 0 || ip > 2) continue;
        k[2 + 2 * ip] = strip.bounds.first;
        k[3 + 2 * ip] = strip.bounds.second;
    }
    return k;
}

Cascade::XdetLevelCols Cascade::xdet_level_columns(const XdetColumns& xc, size_t level,
                                                   const std::vector<IBlob::pointer>& blobs, const Level& lev)
{
    if (level >= xc.levels.size()) {
        THROW(ValueError() << errmsg{String::format("xdet columns %s: no level %d (the file has %d)", xc.path,
                                                    (int) level, (int) xc.levels.size())});
    }
    const auto& tab = xc.levels[level];
    const size_t N = blobs.size(), M = lev.wq.size();
    XdetLevelCols out;
    out.xb_cols.assign(N * xdet_blob_cols, 0.0f);
    out.xw_cols.assign(M * xdet_wire_cols, 0.0f);
    out.live.assign(M, 1.0f);
    size_t nmiss_b = 0, nmiss_w = 0;
    std::string first_miss;
    for (size_t i = 0; i < N; ++i) {
        const auto k = xdet_blob_key(blobs[i]);
        auto it = tab.blob.find(k);
        if (it == tab.blob.end()) {
            if (nmiss_b++ == 0 && first_miss.empty()) {
                first_miss = String::format("blob (face %d, slice %d, %d %d %d %d %d %d)", k[0], k[1], k[2], k[3], k[4],
                                            k[5], k[6], k[7]);
            }
            continue;
        }
        std::copy(it->second.begin(), it->second.end(), out.xb_cols.begin() + xdet_blob_cols * i);
    }
    for (size_t m = 0; m < M; ++m) {
        const XdetWireKey k{(int32_t) lev.wplane[m], (int32_t) lev.wchan[m]};
        auto it = tab.wire.find(k);
        if (it == tab.wire.end()) {
            if (nmiss_w++ == 0 && first_miss.empty()) {
                first_miss = String::format("wire (plane %d, channel %d)", k[0], k[1]);
            }
            continue;
        }
        std::copy(it->second.begin(), it->second.begin() + xdet_wire_cols, out.xw_cols.begin() + xdet_wire_cols * m);
        out.live[m] = it->second[xdet_wire_cols];
    }
    if (nmiss_b || nmiss_w) {
        THROW(ValueError() << errmsg{String::format("xdet columns %s: level %d misses %d of %d blob and %d of %d wire keys; first %s",
                                                    xc.path, (int) level, (int) nmiss_b, (int) N, (int) nmiss_w,
                                                    (int) M, first_miss)});
    }
    return out;
}
