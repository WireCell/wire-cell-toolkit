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
