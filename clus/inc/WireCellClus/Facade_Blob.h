/** A facade over a PC tree giving semantics to otherwise nodes.
 *
 */

#ifndef WIRECELL_CLUS_FACADEBLOB
#define WIRECELL_CLUS_FACADEBLOB

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Spdlog.h"
#include "WireCellUtil/Graph.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"

#include "WireCellClus/Facade_Util.h"


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {

    class Cluster;

    struct BlobCache {
        float_t charge{0};
        float_t center_x{0};
        float_t center_y{0};
        float_t center_z{0};
        int_t npoints{0};

        int_t slice_index_min{0};  // unit: tick
        int_t slice_index_max{0};

        int_t u_wire_index_min{0};
        int_t u_wire_index_max{0};
        int_t v_wire_index_min{0};
        int_t v_wire_index_max{0};
        int_t w_wire_index_min{0};
        int_t w_wire_index_max{0};

        int_t max_wire_interval{-1}; 
        int_t min_wire_interval{-1};
        int_t max_wire_type{-1}; // 0: u, 1: v, 2: w
        int_t min_wire_type{-1}; // 0: u, 1: v, 2: w

        std::vector<geo_point_t> corners_;
    };


    /// Give a node "Blob" semantics
    class Blob : public NaryTree::Facade<points_t>, public Mixin<Blob, BlobCache> {
       public:
        Blob() : Mixin<Blob, BlobCache>(*this, "scalar") {}
        virtual ~Blob() {}

        // Return the cluster to which this blob is a child.  May be nullptr.
        Cluster* cluster();
        const Cluster* cluster() const;

        geo_point_t center_pos() const;

        bool overlap_fast(const Blob& b, const int offset) const;

        float_t charge() const { return cache().charge; }
        float_t center_x() const { return cache().center_x; }
        float_t center_y() const { return cache().center_y; }
        float_t center_z() const { return cache().center_z; }
        int_t npoints() const { return cache().npoints; }

        // units are number of ticks
        /// FIXME: change min, max to begin end
        int_t slice_index_min() const { return cache().slice_index_min; }
        int_t slice_index_max() const { return cache().slice_index_max; }

        int_t u_wire_index_min() const { return cache().u_wire_index_min; }
        int_t u_wire_index_max() const { return cache().u_wire_index_max; }
        int_t v_wire_index_min() const { return cache().v_wire_index_min; }
        int_t v_wire_index_max() const { return cache().v_wire_index_max; }
        int_t w_wire_index_min() const { return cache().w_wire_index_min; }
        int_t w_wire_index_max() const { return cache().w_wire_index_max; }

        int_t get_max_wire_interval() const { return cache().max_wire_interval;}
        int_t get_min_wire_interval() const { return cache().min_wire_interval;}
        int_t get_max_wire_type() const { return cache().max_wire_type;}
        int_t get_min_wire_type() const { return cache().min_wire_type;}

        std::vector<geo_point_t> corners() const { return corners_; }

        // Return a value representing the content of this blob.
        size_t hash() const;

        // Return the scope points.
        std::vector<geo_point_t> points() const;
        size_t nbpoints() const;

        // Check facade consistency
        bool sanity(Log::logptr_t log = nullptr) const;

       private:
        // moved to cache

        std::vector<geo_point_t> corners_;

       protected:
        virtual void fill_cache(BlobCache& cache) const;
    };
    std::ostream& operator<<(std::ostream& os, const Blob& blob);

    // Return true if a is less than b.  May be used as 3rd arg in std::sort to
    // get ascending order.  For descending, pass to sort() rbegin()/rend()
    // instead of begin()/end()..
    bool blob_less(const Blob* a, const Blob* b);
    // Apply standard sort to put blobs in descending order.

    void sort_blobs(std::vector<const Blob*>& blobs);
    void sort_blobs(std::vector<Blob*>& blobs);

    struct blob_less_functor {
        bool operator()(const Facade::Blob* a, const Facade::Blob* b) const {
        return Facade::blob_less(a, b);  
        }
    };

}  // namespace WireCell::PointCloud::Facade

template <> struct fmt::formatter<WireCell::PointCloud::Facade::Blob> : fmt::ostream_formatter {};

#endif
