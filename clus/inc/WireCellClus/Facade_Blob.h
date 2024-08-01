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

    /// Give a node "Blob" semantics
    class Blob : public NaryTree::Facade<points_t> {
       public:
        Blob() = default;
        virtual ~Blob() {}

        // Return the cluster to which this blob is a child.  May be nullptr.
        Cluster* cluster();
        const Cluster* cluster() const;

        geo_point_t center_pos() const;

        bool overlap_fast(const Blob& b, const int offset) const;

        float_t charge() const { return charge_; }
        float_t center_x() const { return center_x_; }
        float_t center_y() const { return center_y_; }
        float_t center_z() const { return center_z_; }
        int_t npoints() const { return npoints_; }

        // units are number of ticks
        /// FIXME: change min, max to begin end
        int_t slice_index_min() const { return slice_index_min_; }
        int_t slice_index_max() const { return slice_index_max_; }

        int_t u_wire_index_min() const { return u_wire_index_min_; }
        int_t u_wire_index_max() const { return u_wire_index_max_; }
        int_t v_wire_index_min() const { return v_wire_index_min_; }
        int_t v_wire_index_max() const { return v_wire_index_max_; }
        int_t w_wire_index_min() const { return w_wire_index_min_; }
        int_t w_wire_index_max() const { return w_wire_index_max_; }

        int_t get_max_wire_interval() const { return max_wire_interval_;}
        int_t get_min_wire_interval() const { return min_wire_interval_;}
        int_t get_max_wire_type() const { return max_wire_type_;}
        int_t get_min_wire_type() const { return min_wire_type_;}

        // Return a value representing the content of this blob.
        size_t hash() const;

        // Return the scope points.
        std::vector<geo_point_t> points() const;
        size_t nbpoints() const;

        // Check facade consistency
        bool sanity(Log::logptr_t log = nullptr) const;

       private:
        float_t charge_{0};
        float_t center_x_{0};
        float_t center_y_{0};
        float_t center_z_{0};
        int_t npoints_{0};

        int_t slice_index_min_{0};  // unit: tick
        int_t slice_index_max_{0};

        int_t u_wire_index_min_{0};
        int_t u_wire_index_max_{0};
        int_t v_wire_index_min_{0};
        int_t v_wire_index_max_{0};
        int_t w_wire_index_min_{0};
        int_t w_wire_index_max_{0};

        int_t max_wire_interval_{-1}; 
        int_t min_wire_interval_{-1};
        int_t max_wire_type_{-1}; // 0: u, 1: v, 2: w
        int_t min_wire_type_{-1}; // 0: u, 1: v, 2: w

       protected:
        // Receive notification when this facade is created on a node.
        virtual void on_construct(node_type* node);
    };
    std::ostream& operator<<(std::ostream& os, const Blob& blob);

    // Return true if a is less than b.  May be used as 3rd arg in std::sort to
    // get ascending order.  For descending, pass to sort() rbegin()/rend()
    // instead of begin()/end()..
    bool blob_less(const Blob* a, const Blob* b);
    // Apply standard sort to put blobs in descending order.
    void sort_blobs(std::vector<const Blob*>& blobs);
    void sort_blobs(std::vector<Blob*>& blobs);

}  // namespace WireCell::PointCloud::Facade

template <> struct fmt::formatter<WireCell::PointCloud::Facade::Blob> : fmt::ostream_formatter {};

#endif
