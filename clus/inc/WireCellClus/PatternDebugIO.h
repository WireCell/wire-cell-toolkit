#ifndef WIRECELLCLUS_PATTERNDEBUG_IO
#define WIRECELLCLUS_PATTERNDEBUG_IO

#include "WireCellClus/Facade.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellUtil/PointTree.h"

#include <string>
#include <memory>

namespace WireCell::Clus::PR::DebugIO {

    /// Dump the inputs to init_first_segment() as a JSON file.
    /// Call this from TaggerCheckNeutrino::visit() before init_first_segment.
    void dump_init_first_segment_inputs(
        const std::string& output_path,
        const Facade::Cluster& cluster,
        const Facade::Cluster* main_cluster,
        bool flag_back_search,
        const TrackFitting& track_fitter);

    /// Data reconstructed from a dump file for test replay.
    struct LoadedTestData {
        // Owns the tree; must outlive cluster/main_cluster pointers.
        std::unique_ptr<PointCloud::Tree::Points::node_t> grouping_node;

        Facade::Cluster* cluster{nullptr};
        Facade::Cluster* main_cluster{nullptr};  // may equal cluster
        bool flag_back_search{true};
        TrackFitting::Parameters trackfitting_params;

        // Precomputed boundary steiner graph indices (avoids needing anode)
        std::pair<size_t, size_t> boundary_steiner_indices{0, 0};
    };

    /// Load a dump file and reconstruct minimal Facade::Cluster objects
    /// suitable for calling init_first_segment().
    LoadedTestData load_init_first_segment_inputs(const std::string& input_path);

}

#endif
