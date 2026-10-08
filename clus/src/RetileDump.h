// doc pdvd/129 phase 1.  Study-only dump of what the Steiner-stage retiler
// (ImproveCluster_2, full mode, final pass) tiles, BEFORE its remove_bad_blobs
// filter.  Called only when the retiler's "dump_dir" knob is set (C++ default
// "" = off: no call, no file, output byte-identical).  Per cluster and anode:
//   <dir>/retile-c<ident>-k<serial>-apa<apa>.tar.gz   the tiled blobs as an
//       imaging cluster graph in the ClusterFileSink numpy layout (prefix
//       "cluster"), readable by ClusterFileSource;
//   <dir>/retile-c<ident>-k<serial>-apa<apa>.json     per-blob status and the
//       provenance of the cells that carry no measured charge.
#ifndef WIRECELLCLUS_RETILEDUMP
#define WIRECELLCLUS_RETILEDUMP

#include "WireCellIface/IBlob.h"
#include "WireCellUtil/RayHelpers.h"

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace WireCell::Clus::RetileDump {

    // the retiler's activity map: (start tick, end tick) -> 5 layers (2 dummy, U, V, W) of per-wire charge
    using activity_map_t = std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t>>;
    // (plane, start tick, wire)
    using cell_t = std::array<int, 3>;

    // The retiler marks every admitted cell without measured charge with this
    // value (improvecluster_1.cxx get_activity_improved / hack_activity_improved).
    inline bool is_sentinel(double v) { return v == 1e-3; }

    /// All cells of the map holding the sentinel.
    std::set<cell_t> sentinel_cells(const activity_map_t& msm);

    /// Rows (plane, start tick, wire, class) for the sentinel cells of `msm`:
    /// class 1 = already a sentinel in `before` (dead channel, or a footprint
    /// cell with no measured charge), class 2 = new since (painted).
    std::vector<std::array<int, 4>> classify(const std::set<cell_t>& before, const activity_map_t& msm);

    struct Face {
        int face{0};
        int tick_span{0};
        WireCell::IBlob::vector iblobs;                  // all tiled blobs, tiling order
        std::vector<char> sampled;                       // per iblob: got points (entered the filter)
        std::vector<std::array<int, 7>> removed;         // start tick, u0,u1,v0,v1,w0,w1 of the blobs the filter removed
        std::vector<std::array<int, 7>> orig;            // the same for the cluster's own blobs on this face
        std::vector<std::array<int, 4>> cells;           // classify()
    };

    void write(const std::string& dir, int cluster_ident, int serial, int apa, const std::vector<Face>& faces);

}  // namespace WireCell::Clus::RetileDump

#endif
