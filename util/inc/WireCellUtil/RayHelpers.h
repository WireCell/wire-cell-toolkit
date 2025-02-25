#ifndef WIRECELLUTIL_RAYHELPERS
#define WIRECELLUTIL_RAYHELPERS

#include "WireCellUtil/Units.h"
#include "WireCellUtil/RayTiling.h"

namespace WireCell::RayGrid {

    // Make a set of ray pairs suitable for construting a
    // RayGrid::Coordinates on the assumption of induction wire
    // planes symmetric at +/- angle relative to the third collection
    // plane.
    WireCell::ray_pair_vector_t symmetric_raypairs(double width = 100*units::cm,
                                                   double height = 100*units::cm,
                                                   double pitch_mag = 3*units::mm,
                                                   double angle = 60.0 * units::pi / 180.0);

    // Trivial, aka bogus, "simulation".  Points "hits" closest wire.
    using measure_t = std::vector<Activity::value_t>;
    std::vector<measure_t> make_measures(const Coordinates& coords, const std::vector<Point>& points);

    // Turn "measures" into activities
    activities_t make_activities(const Coordinates& coords, std::vector<measure_t>& measures);

}
#endif
