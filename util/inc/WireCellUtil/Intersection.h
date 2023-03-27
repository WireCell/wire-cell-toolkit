#ifndef WIRECELLUTIL_INTERSECTION
#define WIRECELLUTIL_INTERSECTION

#include "WireCellUtil/Point.h"

namespace WireCell {

    /** Find intersection of a ray and the sides of an axis-aligned
     * box which are orthogonal to the given axis.
     *
     * \param axis0 is the axis number (0,1,2).
     *
     * \param bounds is a ray from opposite corners of a box defining a CLOSED interval in 3-space.
     *
     * \param point is the origin point of a ray
     *
     * \param dir is the unit vector along the direction of a ray
     *
     * \param hits provides a pair of intersection(s).
     *
     * \return a "hit mask" categorizes the intersections:
     *
     * - 0 :: the ray does not intersect, hits is not valid.
     * 
     * - 1 :: the ray intersects at a single point (box corner), hits.first is valid.
     *
     * - 2 :: the ray intersects at a single point (box corner), hits.second is valid.
     *
     * - 3 :: the ray intersects at two points, hits is fully valid.
     *
     * Invalid points in hits are set to the origin point.
     *
     * The hits ray is made parallel to the ray "dir".
     *
     */
    int box_intersection(int axis0, const Ray& bounds, const Vector& point, const Vector& dir, Ray& hits);

    /** Determine if a ray hits a rectangular box aligned with the
     * Cartesian axes.  See above variant of this function for
     * arguments and return value description.
     */
    int box_intersection(const Ray& bounds, const Vector& point, const Vector& dir, Ray& hits);

}  // namespace WireCell

#endif
