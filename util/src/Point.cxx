#include "WireCellUtil/Point.h"

#include <algorithm>  // minmax

namespace instantiations {

    WireCell::Point a_point;
    WireCell::Ray a_ray;
    WireCell::PointVector a_point_vector;
    WireCell::PointValue a_point_value;
    WireCell::PointValueVector a_point_value_vector;
    WireCell::PointF a_sad_little_point;

}  // namespace instantiations

// std::ostream& operator<<(std::ostream& os, const WireCell::Ray& ray)
// {
//     os << "["  << ray.first << " --> " << ray.second << "]";
//     return os;
// }

bool WireCell::ComparePoints::operator()(const WireCell::Point& lhs, const WireCell::Point& rhs) const
{
    double mag = (lhs - rhs).magnitude();
    if (mag < 1e-10) {
        return false;  // consider them equal
    }
    // otherwise, order them by x,y,z
    for (int ind = 0; ind < 3; ++ind) {
        if (lhs[ind] < rhs[ind]) {
            return true;
        }
    }
    return false;
}

bool WireCell::point_contained(const WireCell::Point& point, const WireCell::Ray& bounds)
{
    for (int axis = 0; axis < 3; ++axis) {
        if (!point_contained(point, bounds, axis)) {
            return false;
        }
    }
    return true;
}

bool WireCell::point_contained(const WireCell::Point& point, const WireCell::Ray& bounds, int axis)
{
    std::pair<double, double> mm = std::minmax(bounds.first[axis], bounds.second[axis]);
    return mm.first <= point[axis] && point[axis] <= mm.second;
}

double WireCell::point_angle(const WireCell::Vector& axis, const WireCell::Vector& vector)
{
    return acos(axis.dot(vector));
}

double WireCell::ray_length(const WireCell::Ray& ray) { return (ray.second - ray.first).magnitude(); }

WireCell::Vector WireCell::ray_vector(const WireCell::Ray& ray) { return ray.second - ray.first; }
WireCell::Vector WireCell::ray_unit(const WireCell::Ray& ray) { return ray_vector(ray).norm(); }

WireCell::Ray WireCell::ray_pitch(const WireCell::Ray& pu, const WireCell::Ray& qv)
{
    // http://geomalgorithms.com/a07-_distance.html
    const WireCell::Vector w0 = pu.first - qv.first;
    const WireCell::Vector u = ray_unit(pu);
    const WireCell::Vector v = ray_unit(qv);
    const double a = u.dot(u), b = u.dot(v), c = v.dot(v);
    const double d = u.dot(w0), e = v.dot(w0);

    const double denom = a * c - b * b;
    if (denom < 1e-6) {  // parallel
        double t = e / c;
        return Ray(pu.first, qv.first + t * v);
    }
    const double s = (b * e - c * d) / denom;
    const double t = (a * e - b * d) / denom;
    return Ray(pu.first + s * u, qv.first + t * v);
}

double WireCell::ray_dist(const WireCell::Ray& ray, const WireCell::Point& point)
{
    return ray_unit(ray).dot(point - ray.first);
}

double WireCell::ray_volume(const WireCell::Ray& ray)
{
    auto diff = ray_vector(ray);
    return diff.x() * diff.y() * diff.z();
}

// This returns the intersection between two unrotated(!)
// box volumes, represented by a Ray.
//       ________
//      |    S2  |
//  ____|___     |
// |    |   |    |
// |    |___|____|
// | S1     |
// |________|
WireCell::Ray WireCell::box_intersect(const Ray& s1, const Ray& s2)
{
    Ray bb_ray;
    for (size_t ind = 0; ind < 3; ind++) {
        auto lb1 = s1.first[ind];            // left bound
        auto rb1 = s1.second[ind];           // right bound
        if (lb1 > rb1) std::swap(lb1, rb1);  // let left < right
        auto lb2 = s2.first[ind];
        auto rb2 = s2.second[ind];
        if (lb2 > rb2) std::swap(lb2, rb2);

        bb_ray.first[ind] = std::max(lb1, lb2);
        bb_ray.second[ind] = std::min(rb1, rb2);
    }
    return bb_ray;
}

bool WireCell::plane_split(const WireCell::Point& point, const WireCell::Vector& normal, const WireCell::Ray& ray)
{
    const auto& [v1,v2] = ray;

    // If projection onto normal are different signs, then the ray is split by
    // the plane.
    return normal.dot(v1-point) * normal.dot(v2-point) < 0;
}

WireCell::Point WireCell::plane_intersection(const WireCell::Point& point, const WireCell::Vector& normal, const WireCell::Ray& segment)
{
    const double denom = normal.dot(segment.second - segment.first);
    if (denom == 0) {           // line is parallel to the plane
        return point + normal;
    }

    double t = normal.dot(point - segment.first) / denom;
    return segment.first + t*(segment.second - segment.first);
}
