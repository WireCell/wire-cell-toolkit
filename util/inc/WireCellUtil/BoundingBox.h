#ifndef WIRECELL_BOUNDINGBOX
#define WIRECELL_BOUNDINGBOX

#include "WireCellUtil/Point.h"

namespace WireCell {

    /** A bounding box parallel to the Cartesian axes.
     */
    class BoundingBox {
        Ray m_bounds;
        bool m_initialized = false;

       public:
        /// Create a bounding box without an initial point or ray
        BoundingBox() { m_bounds = Ray(); }

        /// Create a bounding box bounding an initial point.
        BoundingBox(const Point& initial);

        /// Create a bounding box bounding an initial ray.
        BoundingBox(const Ray& initial);

        /// Create a bounding box from an iterator pair.
        template <typename RayOrPointIterator>
        BoundingBox(const RayOrPointIterator& begin, const RayOrPointIterator& end)
        {
            for (auto it = begin; it != end; ++it) {
                (*this)(*it);
            }
        }

        // Create a union of this bounding box with another one
        BoundingBox unite(const BoundingBox& other) const;

        // Create an intersection of this bounding box with another one
        BoundingBox intersect(const BoundingBox& other) const;

        /// Return the intersection points of an infinite line going through the
        /// raw with the bounding box.
        ///
        /// The returned ray points in the same direction as the line.
        /// 
        /// If the line does not intersect the default Ray is returned (both
        /// points are at origin).
        ///
        /// In the unusual case that the line clips a corner, both ray endpoints
        /// are identical.  Note: this can be ambiguous with the returning the
        /// default Ray!
        ///
        Ray intersect(const Ray& line) const;

        /// Just like intersect(line) but consider a finite line segment.
        ///
        /// If an endpoint is inside the bounding box, it is retained in its
        /// place in the returned ray.
        Ray crop(const Ray& segment) const;


        // Create a bounding box that contains all provided bounding boxes
        static BoundingBox combine(const std::vector<BoundingBox>& boxes);

        // Check if this bounding box overlaps with another one
        bool overlaps(const BoundingBox& other) const;

                
        // Calculate the minimum distance from a point to the bounding box (0 if inside)
        double distance(const Point& point) const;

        /// Return distances FROM a point TO walls ALONG an axis.
        ///
        /// When the direction along the axis from the point intercepts the
        /// bounding box, two numbers are returned.  Otherwise the vector is
        /// empty.
        ///
        /// The first number gives the distance to the wall with the smaller
        /// axis coordinate.
        ///
        /// Distances are SIGNED.  A point inside the BB will have the first
        /// number negative, second positive.  Two positives means the point is
        /// outside the BB and at a lower coordinate.  The BB is "in front" of
        /// the point, along the axis.  Two negatives means the point is outside
        /// of the BB and the BB is "behind" the point.
        std::vector<double> axis_distances(const Point& point, int axis) const;

        /// Returns the closest point on the box to the given point
        Point closest_point(const Point& point) const;
       
        /// Get the volume of the box (0 if any dimension is 0)
        double volume() const;

        /// Get the center point of the box
        Point center() const;

        /// Get the dimensions of the box as a vector (width, height, depth)
        Vector dimensions() const;


        /// Return true if point is inside bounding box
        bool inside(const Point& point) const;

        /// Return the ray representing the bounds.
        const Ray& bounds() const { return m_bounds; }

        /// Enlarge bounds to hold point.
        void operator()(const Point& p);

        /// Enlarge bounds to hold ray.
        void operator()(const Ray& r);

        template <typename RayOrPointIterator>
        void operator()(const RayOrPointIterator& begin, const RayOrPointIterator& end)
        {
            for (auto it = begin; it != end; ++it) {
                (*this)(*it);
            }
        }

        // Pad the BB by adding the vector relative * diagonal_vector
        // to max BB point and subtracting it from min.  Negative
        // "relative" value causes the BB to shrink.
        void pad_rel(double relative);
        // Pad the BB by the vector distance * diagonal_unit_vector
        // otherwise as pad_rel().
        void pad_abs(double distance);

        bool empty() const { return !m_initialized; }
    };


}  // namespace WireCell

#endif
