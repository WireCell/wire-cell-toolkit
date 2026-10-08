#include "RetileDump.h"

#include "WireCellAux/ClusterArrays.h"
#include "WireCellIface/ICluster.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/custard/custard_boost.hpp"

#include <boost/iostreams/filtering_stream.hpp>

#include <cmath>

using namespace WireCell;

namespace WireCell::Clus::RetileDump {

    std::set<cell_t> sentinel_cells(const activity_map_t& msm)
    {
        std::set<cell_t> out;
        for (const auto& [key, measures] : msm) {
            for (size_t layer = 2; layer < measures.size(); ++layer) {
                const auto& m = measures[layer];
                for (size_t w = 0; w < m.size(); ++w) {
                    if (is_sentinel(m[w])) out.insert({(int) layer - 2, key.first, (int) w});
                }
            }
        }
        return out;
    }

    std::vector<std::array<int, 4>> classify(const std::set<cell_t>& before, const activity_map_t& msm)
    {
        std::vector<std::array<int, 4>> out;
        for (const auto& c : sentinel_cells(msm)) {
            out.push_back({c[0], c[1], c[2], before.count(c) ? 1 : 2});
        }
        return out;
    }

    // BlobClustering.cxx add_slice / add_blobs, duplicated (fork rule: the production file is untouched)
    static void add_slice(cluster_indexed_graph_t& grind, const ISlice::pointer& islice)
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
    static void add_blobs(cluster_indexed_graph_t& grind, const IBlob::vector& iblobs)
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

    template <typename Rows> static Json::Value rows_json(const Rows& rows)
    {
        Json::Value out = Json::arrayValue;
        for (const auto& r : rows) {
            Json::Value jr = Json::arrayValue;
            for (int v : r) jr.append(v);
            out.append(jr);
        }
        return out;
    }

    void write(const std::string& dir, int cluster_ident, int serial, int apa, const std::vector<Face>& faces,
               const std::string& prefix)
    {
        const std::string base = String::format("%s/%s-c%d-k%d-apa%d", dir, prefix, cluster_ident, serial, apa);

        cluster_indexed_graph_t grind;
        Json::Value jfaces = Json::arrayValue;
        size_t nblobs = 0;
        for (const auto& f : faces) {
            add_blobs(grind, f.iblobs);
            nblobs += f.iblobs.size();
            Json::Value jf;
            jf["face"] = f.face;
            jf["tick_span"] = f.tick_span;
            // start tick, u0,u1,v0,v1,w0,w1, sampled
            Json::Value jb = Json::arrayValue;
            for (size_t i = 0; i < f.iblobs.size(); ++i) {
                const auto& ib = f.iblobs[i];
                Json::Value r = Json::arrayValue;
                r.append((int) std::lround(ib->slice()->start() / ib->slice()->span() * f.tick_span));
                int b[6] = {0, 0, 0, 0, 0, 0};
                for (const auto& strip : ib->shape().strips()) {
                    const int p = strip.layer - 2;
                    if (p < 0 || p > 2) continue;
                    b[2 * p] = strip.bounds.first;
                    b[2 * p + 1] = strip.bounds.second;
                }
                for (int v : b) r.append(v);
                r.append(i < f.sampled.size() && f.sampled[i] ? 1 : 0);
                jb.append(r);
            }
            jf["blobs"] = jb;
            jf["removed"] = rows_json(f.removed);
            jf["orig"] = rows_json(f.orig);
            jf["cells"] = rows_json(f.cells);
            jfaces.append(jf);
        }
        Json::Value top;
        top["cluster"] = cluster_ident;
        top["serial"] = serial;
        top["apa"] = apa;
        top["faces"] = jfaces;
        Persist::assuredir(base + ".json");
        Persist::dump(base + ".json", top, false);
        if (nblobs == 0) {
            return;
        }

        boost::iostreams::filtering_ostream out;
        custard::output_filters(out, base + ".tar.gz");
        if (out.empty()) {
            raise<ValueError>("RetileDump: unsupported outname: %s.tar.gz", base.c_str());
        }
        Aux::ClusterArrays::node_array_set_t nas;
        Aux::ClusterArrays::edge_array_set_t eas;
        Aux::ClusterArrays::to_arrays(grind.graph(), nas, eas);
        for (const auto& [nc, na] : nas) {
            std::string name = "_nodes";
            name[0] = nc;
            Stream::write(out, String::format("cluster_%d_%s.npy", cluster_ident, name), na);
        }
        for (const auto& [ec, ea] : eas) {
            Stream::write(out, String::format("cluster_%d_%sedges.npy", cluster_ident, Aux::ClusterArrays::to_string(ec)), ea);
        }
        out.flush();
        out.pop();
    }

}  // namespace WireCell::Clus::RetileDump
