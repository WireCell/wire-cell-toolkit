/** Common types and functions for pattern recognition code.

    This file must not depend on any (other) types in the WireCell::Clus::PR
    namespace.
 */
#ifndef WIRECELL_CLUS_PR_COMMON
#define WIRECELL_CLUS_PR_COMMON

#include "WireCellUtil/Point.h"

namespace WireCell::Clus::PR {

    /// A WCPoint is a 3D point and corresponding wire indices and an index.
    //
    // FIXME: does this need any change given we now support wrapped wires and
    // multi APA/face detectors?
    struct WCPoint {
        WireCell::Point point;   // 3D point 
        int uvw[3] = {-1,-1,-1}; // wire indices
        int index{-1};           // point index in some container

        // FIXME: WCP had this, does WCT need it?
        // blob* b;

        
        // Return true if the point information has been filled.
        bool valid() const {
            if (index < 0) return false;
            return true;
        }
    };
    using WCPointVector = std::vector<WCPoint>;

    /** A Fit holds information predicted about a point by some "fit".
     *
     * A PR::Vertex has a scalar Fit object and PR::Segment has a vector<Fit>
     *
     * Note, WCP's ProtoSegment had struct-of-array instead of vector<Fit>
     */
    struct Fit {
        WireCell::Point point;
        double dQ{-1}, dx{0}, pu{-1}, pv{-1}, pw{-1}, pt{0}, reduced_chi2{-1};

        int index{-1};
        double range{-1};
        bool flag_fix{false};

        // Explicitly NOT defined:

        // bool flag_fit.  This seems never actually used in WCP.  If needed,
        // can we simply test on index or range?

        // Restore values to invalid
        void reset() {
            index = -1;
            flag_fix = false;
            range = -1;
        }

        double distance(const Point& p) {
            return (p - point).magnitude();
        }

        /** Return true if fit information has been filled */
        bool valid() {
            if (index < 0 || range < 0) return false;
            return true;
        }
    };
    using FitVector = std::vector<Fit>;
    
    /** Some mixin classes, eg used by Vertex and Segment.

        See also "Graphed" from PRGraph.h and "Flagged" from util.
     */

    /// Transform an object-with-point to a point.
    template<typename OWP>
    Point owp_to_point(const OWP& owp) { return owp.point; };

    /// This type describes a `transform` function from some type to type Point.
    ///
    /// Functions that operate on Point and templated types and that take this
    /// `transform` function will apply it to non-Point types in order to
    /// produce a Point.
    ///
    /// A likely use will be to pass one of these:
    ///
    /// @code{.cpp}
    /// transform = owp_to_point<Fit>
    /// transform = owp_to_point<WCPoint>
    /// @endcode
    template<typename Val>
    using to_point_f = std::function<Point(const Val&)>;

    /// Return the closest to point from a collection of points.
    ///
    /// An iterator into the collection is returned.
    ///
    /// See `to_point_f` type for information about the `transform` argument.
    template<typename Vec>
    typename Vec::const_iterator closest_point(
        const Vec& points, const Point& point,
        to_point_f<typename Vec::value_type> transform = [](const typename Vec::value_type& a) { return a; })
    {
        return std::min_element(points.begin(), points.end(),
                                [&](const auto& a, const auto& b) {
                                    return (transform(a)-point).magnitude() < (transform(b)-point).magnitude();});
    }

    /// Return the closest to point from an iteration range of points.
    ///
    /// An iterator into the collection is returned.
    ///
    /// See `to_point_f` type for information about the `transform` argument.
    template<typename It>
    It closest_point(
        It begin, It end,
        const Point& point,
        to_point_f<typename std::iterator_traits<It>::value_type> transform =
        [](const typename std::iterator_traits<It>::value_type& a) { return a; })
    {
        return std::min_element(begin, end,
                                [&](const auto& a, const auto& b) {
                                    return (transform(a)-point).magnitude() < (transform(b)-point).magnitude();});
    }



    /// Return the "walk length" over a path of points in a vector-like collection.
    ///
    /// This returns:
    ///
    ///  length(v[0],v[1]) + length(v[1],v[2]) ...
    ///
    /// See `to_point_f` type for information about the `transform` argument.
    template<typename VP>
    double walk_length(const VP& points,
                       to_point_f<VP> transform = [](const typename VP::value_type& a) { return a; }) {
        const auto siz = points.size();
        if (siz < 2) { return 0.0; }
        double total_dist = 0.0;
        Point last_point = transform(points[0]);
        for (size_t i = 1; i < siz - 1; ++i) {
            Point next_point = transform(points[i+1]);
            total_dist += (last_point - next_point).magnitude();
            last_point = next_point;
        }
        return total_dist;
    }

    /// Return the "walk length" over a path of points in an integrator range.
    /// 
    /// This returns:
    ///
    ///  length(v[0],v[1]) + length(v[1],v[2]) ...
    ///
    /// See `to_point_f` type for information about the `transform` argument.
    template<typename It>
    double walk_length(It begin, It end, 
                       to_point_f<typename std::iterator_traits<It>::value_type> transform =
                       [](const typename std::iterator_traits<It>::value_type& a) { return a; })
    {
        const auto siz = std::distance(begin, end);
        if (siz < 2) { return 0.0; }
        double total_dist = 0.0;
        Point last_point = transform(*begin);
        for (++begin; begin != end; ++begin) {
            Point next_point = transform(*begin);
            total_dist += (last_point - next_point).magnitude();
            last_point = next_point;
        }
        return total_dist;
    }

    

}


#endif
