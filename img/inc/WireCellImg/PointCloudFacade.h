/** A facade over a PC tree giving semantics to otherwise nodes.
 *
 */

#ifndef WIRECELLIMG_POINTCLOUDFACADE
#define WIRECELLIMG_POINTCLOUDFACADE

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {
    using node_t = WireCell::PointCloud::Tree::Points::node_t;
    using node_ptr = std::unique_ptr<node_t>;
    using geo_point_t = WireCell::Point;
    using geo_vector_t = WireCell::Vector;
    using float_t = double;
    using int_t = int;

    struct TPCParams {
        float_t pitch_u {3*units::mm};
        float_t pitch_v {3*units::mm};
        float_t pitch_w {3*units::mm};
        float_t angle_u {1.0472};  // 60 degrees    uboone geometry ...
        float_t angle_v {-1.0472};  //-60 degrees   uboone geometry ...
        float_t angle_w {0};        // 0 degrees    uboone geometry ...
        float_t tick_width {0.5*1.101*units::mm}; // width corresponding to one tick time
    };

    // Provide common types for an object to be shared via pointer.
    template<typename T>
    struct Shared {

        // The wrapped sub class type.
        using shared_type = T;

        // For holding a facade by a shared pointer.
        using pointer = std::shared_ptr<T>;
        using const_pointer = std::shared_ptr<const T>;

        // Simple collection of shared facades.
        using vector = std::vector<pointer>;
        using const_vector = std::vector<const_pointer>;
    };
    

    class Blob : public Shared<Blob> {
       public:
        Blob(node_t* n);

        geo_point_t center_pos() const;
	int_t num_points() const;
        bool overlap_fast(const Blob& b, const int offset) const;

        /// FIXME: cache all scalers?
        float_t charge {0};
        float_t center_x {0};
        float_t center_y {0};
        float_t center_z {0};
	int_t npoints {0};
	
        int_t slice_index_min {0}; // unit: tick
        int_t slice_index_max {0};

        int_t u_wire_index_min {0};
        int_t u_wire_index_max {0};
        int_t v_wire_index_min {0};
        int_t v_wire_index_max {0};
        int_t w_wire_index_min {0};
        int_t w_wire_index_max {0};

        node_t* node() { return m_node; }
        const node_t* node() const { return m_node; }

      private:
        node_t* m_node;  /// do not own
    };

    // A cluster facade adds to a PC tree node semantics of a set of blobs
    // likely due to connected activity.
    class Cluster : public Shared<Cluster> {
        Blob::vector m_blobs;
        node_t* m_node;  /// do not own

        // The expected scope.
        const WireCell::PointCloud::Tree::Scope scope = { "3d", {"x","y","z"} };

    public:
        Cluster(node_t* n);

        node_t* node() { return m_node; }
        const node_t* node() const { return m_node; }

        // Access the collection of blobs.
        Blob::const_vector blobs() const {
            Blob::const_vector ret(m_blobs.size());
            std::transform(m_blobs.begin(), m_blobs.end(), ret.begin(),
                           [](auto& bptr) { return std::const_pointer_cast<const Blob>(bptr); });
            return ret;
        }
        Blob::vector blobs() { return m_blobs; }

        geo_point_t calc_ave_pos(const geo_point_t& origin, const double dis, const int alg = 0) const;

        // Return blob containing the returned point that is closest to the given point.
	std::pair<geo_point_t, Blob::const_pointer > get_closest_point_mcell(const geo_point_t& point) const;

        // Return set of blobs each with an a characteristic point.  The set
        // includes blobs with at least one point within the given radius of the
        // given point.  The characteristic point is the point in the blob that
        // is closest to to the given point.
        //
        // Note: radius must provide a LINEAR distance measure.
	std::map<Blob::const_pointer, geo_point_t> get_closest_mcell(const geo_point_t& point, double radius) const;

	std::pair<geo_point_t, double> get_closest_point_along_vec(geo_point_t& p_test, geo_point_t dir, double test_dis, double dis_step, double angle_cut, double dis_cut) const;

        // Return the number of points in the k-d tree
	int get_num_points() const;

        // Return the number of points within radius of the point.  Note, radius
        // is a LINEAR distance through the L2 metric is used internally.
	int get_num_points(const geo_point_t& point, double radius) const;

        // Return the number of points in the k-d tree partitioned into pair
        // (#forward,#backward) based on given direction of view from the given
        // point.
	std::pair<int, int> get_num_points(const geo_point_t& point, const geo_point_t& dir) const;

        // Return the number of points with in the radius of the given point in
        // the k-d tree partitioned into pair (#forward,#backward) based on
        // given direction of view from the given point.
        //
        // Note: the radius is a LINEAR distance measure.
	std::pair<int, int> get_num_points(const geo_point_t& point, const geo_point_t& dir, double radius) const;

        // Return the points at the extremes of the X-axis.
        //
        // Note: the two points are in ASCENDING order!
	std::pair<geo_point_t, geo_point_t> get_earliest_latest_points() const;

        // Return the points at the extremes of the given Cartesian axis.  Default is Y-axis.
        //
        // Note: the two points are in DESCENDING order!
	std::pair<geo_point_t, geo_point_t> get_highest_lowest_points(size_t axis=1) const;
	
	
        Blob::const_vector is_connected(const Cluster& c, const int offset) const;

        // Return the angles characterizing the points within radius of given point.
        // The angles are pair (cos(theta), phi) if alg is 0 else (theta, phi).
        //
        // Note: radius must provide a LINEAR distance measure.
        std::pair<double, double> hough_transform(const geo_point_t& point, const double radius, const int alg = 1) const;

        // Call hough_transform() and transform result as to a directional vector representation.
        //
        // FIXME: this function should be removed.  Caller should do the
        // directional vector transform themselves.  We may add that function to
        // eg Point.h.
        geo_point_t vhough_transform(const geo_point_t& point, const double radius, const int alg = 1) const;

        // get the number of unique uvwt bins
        std::tuple<int, int, int, int> get_uvwt_range() const;
        double get_length(const TPCParams& tp) const;

	// Return blob at the front of the time blob map.
	Blob::const_pointer get_first_blob() const;

	// Return blob at the back of the time blob map.
	Blob::const_pointer get_last_blob() const;
	
       private:
       // start slice index (tick number) to blob facade pointer
       // can be duplicated, example usage: https://github.com/HaiwangYu/learn-cpp/blob/main/test-multimap.cxx
	std::multimap<int, Blob::const_pointer> m_time_blob_map;
    };


    inline double cal_proj_angle_diff(const geo_vector_t& dir1, const geo_vector_t& dir2, double plane_angle) {
        geo_vector_t temp_dir1;
        geo_vector_t temp_dir2;

        temp_dir1.set(dir1.x(), 0, -sin(plane_angle) * dir1.y() + cos(plane_angle) * dir1.z());
        temp_dir2.set(dir2.x(), 0, -sin(plane_angle) * dir2.y() + cos(plane_angle) * dir2.z());

        return temp_dir1.angle(temp_dir2);
    }

    inline bool is_angle_consistent(const geo_vector_t& dir1, const geo_vector_t& dir2, bool same_direction, double angle_cut, double uplane_angle, double vplane_angle, double wplane_angle, int num_cut=2) {
        double angle_u = cal_proj_angle_diff(dir1, dir2, uplane_angle);
        double angle_v = cal_proj_angle_diff(dir1, dir2, vplane_angle);
        double angle_w = cal_proj_angle_diff(dir1, dir2, wplane_angle);
        int num = 0;
        // input is degrees ...
        angle_cut *= 3.1415926 / 180.;

        if (same_direction) {
            if (angle_u <= angle_cut) num++;
            if (angle_v <= angle_cut) num++;
            if (angle_w <= angle_cut) num++;
        } else {
            if ((3.1415926 - angle_u) <= angle_cut) num++;
            if ((3.1415926 - angle_v) <= angle_cut) num++;
            if ((3.1415926 - angle_w) <= angle_cut) num++;
        }

        if (num >= num_cut) return true;
        return false;
    }

}  // namespace WireCell::PointCloud

#endif
