#ifndef WIRECELLUTIL_POINT_H
#define WIRECELLUTIL_POINT_H

#include "WireCellUtil/D3Vector.h"
#include "WireCellUtil/Configuration.h"

#include <set>
#include <memory>  // auto_ptr
#include <vector>

namespace WireCell {

    /// A 3D Cartesian point in double precision.

    typedef D3Vector<double> Point;

    /// An alias for Point.
    typedef Point Vector;

    /// A line segment running from a first (tail) to a second (head) point.
    typedef std::pair<Point, Point> Ray;

    // Pair of rays.
    typedef std::pair<Ray, Ray> ray_pair_t;

    // A collection of ray pairs;
    typedef std::vector<ray_pair_t> ray_pair_vector_t;

    /// PointVector - a collection of Points
    typedef std::vector<Point> PointVector;

    /// A scalar + vector, eg charge at a point.
    typedef std::pair<double, Vector> ScalarPoint;
    typedef std::vector<ScalarPoint> ScalarField;

    /// PointValue - an association of a point and a value
    typedef std::pair<Point, float> PointValue;

    /// PointValueVector - a collection of point-value associations
    typedef std::vector<PointValue> PointValueVector;

    /// Return true if lhs<rhs w/in tolerance.
    struct ComparePoints {
        bool operator()(const Point& lhs, const Point& rhs) const;
    };
    typedef std::set<Point, ComparePoints> PointSet;

    /// PointF - a 3D Cartesian point in single precision for when
    /// memory is constrained and double precision is not required.
    typedef D3Vector<float> PointF;

    /** Return true if point is contained in a rectangular solid
     * described by the ray bounds running between diagonally opposed
     * corners.*/
    bool point_contained(const Point& point, const Ray& bounds);

    /** Return true if point is contained by the bounding box along
     * the given axis (x=0, y=1, z=2) of the bounding box. */
    bool point_contained(const Point& point, const Ray& bounds, int axis);

    /** Return the angle from axis vector to vector. This is just
     * acos(dot).*/
    double point_angle(const Vector& axis, const Vector& vector);

    /** Return a ray representing the points of closest approach
     * between the two lines colinear with the two rays. */
    Ray ray_pitch(const Ray& ray1, const Ray& ray2);

    /** Return the distance from the tail to the head of the ray. */
    double ray_length(const Ray& ray);

    /** Return a vector going from ray's tail to ray's head. */
    Vector ray_vector(const Ray& ray);

    /** Return a unit vector pointing in the direction from the tail
     * to the head of the ray. */
    Vector ray_unit(const Ray& ray);

    /** Return the distance from the tail of the ray to the point
     * projected onto the ray's direction. */
    double ray_dist(const Ray& ray, const Point& point);

    /** Return the volume of a box aligned with axes and with the ray
     * at opposite corners. */
    double ray_volume(const Ray& ray);

    /** Return the ray representing the intersection of two BoundingBox. */
    Ray box_intersect(const Ray& r1, const Ray& r2);

    template <>
    inline  // fixme: ignores default
        WireCell::Point
        convert<WireCell::Point>(const Configuration& cfg, const WireCell::Point& def)
    {
        return Point(get<double>(cfg, "x"), get<double>(cfg, "y"), get<double>(cfg, "z"));
    }
    template <>
    inline  // fixme: ignores default
        WireCell::Ray
        convert<WireCell::Ray>(const Configuration& cfg, const WireCell::Ray& def)
    {
        return Ray(get<WireCell::Point>(cfg, "tail"), get<WireCell::Point>(cfg, "head"));
    }

    // std::ostream& operator<<(std::ostream& os, const WireCell::Ray& ray);
    inline std::ostream& operator<<(std::ostream& os, const WireCell::Ray& ray)
    {
        os << "[" << ray.first << " --> " << ray.second << "]";
        return os;
    }

    inline WireCell::Ray operator/(WireCell::Ray ray, const double& scale)
    {
        return WireCell::Ray(ray.first / scale, ray.second / scale);
    }

}  // namespace WireCell

// inline
// WireCell::Ray operator*(WireCell::Ray ray, const double& scale) {
//     return WireCell::Ray(ray.first*scale, ray.second*scale);
// }
// inline
// WireCell::Ray operator*(double scale, const WireCell::Ray& ray) {
//     return WireCell::Ray(ray.first*scale, ray.second*scale);
// }

#include <fmt/core.h>
#if FMT_VERSION >= 90000
#include <fmt/ostream.h>
template <> struct fmt::formatter<WireCell::Point> : fmt::ostream_formatter {};
template <> struct fmt::formatter<WireCell::Ray> : fmt::ostream_formatter {};
#endif


#endif
