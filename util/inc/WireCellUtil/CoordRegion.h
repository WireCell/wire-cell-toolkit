#ifndef WIRECELL_COORDREGION
#define WIRECELL_COORDREGION

#include "WireCellUtil/Point.h"
#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    /** Base for all coordinate bounds.

        A bounds always considers a preferred axis and the less than comparison
        assumes that axis.  That is a point is "less than" the bounds if the
        point is on the "left hand side" of the bounds, aka, the side with
        decreasing coordinate values.
     */
    class CoordBounds {
    public:
        virtual ~CoordBounds() {}

        // The distance FROM bounds TO point.  A negative distance means the
        // point is "less than" the bounds.
        virtual double distance(const Point& pt) const = 0;

        // Return a nominal location along the axis.
        virtual double location() const = 0;
    };

    inline bool operator<(const Point& pt, const CoordBounds& cb) {
        return cb.distance(pt) < 0;
    }
    inline bool operator<(const CoordBounds& cb, const Point& pt) {
        return cb.distance(pt) > 0;
    }

    /** A boundary defined on a coordinate axis.

        The axis coordinate of a point is compared to the given bounds coord.
     */
    class CoordBoundsScalar : public CoordBounds{
    private:
        double m_coord{0};
        int m_axis{0};
    public:
        CoordBoundsScalar() {}
        explicit CoordBoundsScalar(double coord, int axis=0)
            : m_coord(coord), m_axis(axis) {}
        virtual ~CoordBoundsScalar() {};

        virtual double distance(const Point& pt) const {
            return pt[m_axis] - m_coord;
        }
        virtual double location() const { return m_coord; }
    };
    /** A boundary defined by a plane.

        A plane is defined by a ray.  The ray tail is interpreted as a point in
        the plane and the ray vector (from tail to head) is interpreted as being
        along the normal.  The absolute direction and the magnitude of the ray
        is not relevant.
     */
    class CoordBoundsRay : public CoordBounds {
    private:
        Point m_origin{};
        Point m_norm{1,0,0};
        int m_axis{0};
    public:
        CoordBoundsRay() {}
        explicit CoordBoundsRay(const Ray& ray, int axis=0) : m_axis(axis){
            m_origin = ray.first;
            m_norm = ray_unit(ray);
            if (m_norm[m_axis] < 0) {
                m_norm = -1.0 * m_norm;
            }
        }
        virtual ~CoordBoundsRay() {};

        virtual double distance(const Point& pt) const {
            return (pt - m_origin).dot(m_norm);
        }
        virtual double location() const { return m_origin[m_axis]; }
    };

    /** A boundary defined by a sampled surface.

        The tree closest points to a test point are used to form a plane.  The
        test point is then compared to be "below" the plane along the axis.
     */

    class CoordBoundsSampled : public CoordBounds {
        using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexStatic>;
        nfkd_t m_kd;
        int m_axes[3];
        double m_location{0};

    public:
        CoordBoundsSampled();

        // Construct with samples as SoA.
        explicit CoordBoundsSampled(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    const std::vector<double>& z, int axis=0)
            : m_kd(3) {

            m_axes[0] = axis;
            m_axes[1] = (axis+1)%3;
            m_axes[2] = (axis+2)%3;

            nfkd_t::points_type kdpts = {x, y, z};
            m_kd.append(kdpts);

            m_location = 0;
            for (const auto& val : kdpts[axis]) {
                m_location += val;
            }
            m_location /= kdpts[axis].size();

        }

        virtual double distance(const Point& pt) const {
            const auto res = m_kd.knn(3, pt);

            const Point a = m_kd.point3d(res[0].first);
            const Point b = m_kd.point3d(res[1].first);
            const Point c = m_kd.point3d(res[2].first);

            const auto v0 = pt - a; // relative vector from a to test point
            const auto v1 = b - a;  // relative vector in plane
            const auto v2 = c - a;  // relative vector in plane
            const auto norm = v1.cross(v2);
            return norm.dot(v0) * (norm[m_axes[0]] < 0 ? -1.0 : 1.0);
        }
        virtual double location() const { return m_location; }

    };

    // Accept a standard schema to make one of the coords.
    std::unique_ptr<CoordBounds> make_coordbounds(const Configuration& cfg, int axis=0) {
        if (cfg.isNumeric()) {
            return std::make_unique<CoordBoundsScalar>(cfg.asDouble(), axis);
        }
        if (cfg.isMember("tail")) {
            return std::make_unique<CoordBoundsRay>(convert<Ray>(cfg), axis);
        }
        if (cfg.isMember("x")) {
            return std::make_unique<CoordBoundsSampled>(
                get<std::vector<double>>(cfg, "x"),
                get<std::vector<double>>(cfg, "y"),
                get<std::vector<double>>(cfg, "z"),
                axis);
        }
        return nullptr;
    }
    

    /** CoordRegion - partitions the axis with a high and low bounds.
     */
    class CoordRegion {
        
    public:
        // Construct a region with two bounds.
        explicit CoordRegion(const CoordBounds& lo, const CoordBounds& hi)
            :m_lo(lo), m_hi(hi) {}

        // Return true if pt is in the open region (lo, hi) or (hi, lo).  Use
        // this when ordering between lo and hi is unknown.
        bool inside(const Point& pt) const {
            return (m_lo < pt && pt < m_hi) || (m_hi < pt && pt < m_lo);
        }

        // Treat lo < hi and check if point is strictly between the two.
        bool between(const Point& pt) const {
            return m_lo < pt && pt < m_hi;
        }

    private:
        const CoordBounds&  m_lo;
        const CoordBounds&  m_hi;
    };

}

#endif
