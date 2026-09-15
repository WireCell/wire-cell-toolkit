// doc pdvd/101 -- the default-OFF lattice knobs of the trajectory fit.
//
// fit_weight_pow (default 2) is the power of the position least-squares weight
// |q/err * factors|^pow; assoc_cont_center (default 0) centres
// form_point_association's wire window on the continuous wire coordinate.  The
// window centre comes from Facade::point2wind_cont, which must round to exactly
// the wire point2wind returns, or the knob-on window would sit on a different
// wire from the one every other part of the fit uses.  These cases pin:
//   * the defaults are the legacy path (2 and 0);
//   * the set/get round trip through the string dispatch (a typo throws);
//   * std::round(point2wind_cont) == point2wind over random points on the three
//     detectors' pitches and wire angles, and a point on wire k's centre has
//     continuous coordinate k.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Util.h"

#include <cmath>
#include <random>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("doc pdvd/101 lattice knobs: defaults are the legacy path")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("fit_weight_pow") == 2.0);
    CHECK(tf.get_parameter("assoc_cont_center") == 0.0);
    TrackFitting::Parameters p;
    CHECK(p.fit_weight_pow == 2.0);
    CHECK(p.assoc_cont_center == 0.0);
}

TEST_CASE("doc pdvd/101 lattice knobs: set_parameter / get_parameter round trip")
{
    TrackFitting tf;
    tf.set_parameter("fit_weight_pow", 1.5);
    tf.set_parameter("assoc_cont_center", 1.0);
    CHECK(tf.get_parameter("fit_weight_pow") == 1.5);
    CHECK(tf.get_parameter("assoc_cont_center") == 1.0);
    CHECK_THROWS_AS(tf.set_parameter("fit_weight_power", 1.5), ValueError);
}

TEST_CASE("doc pdvd/101 point2wind_cont rounds to point2wind")
{
    // (angle, pitch) of SBND, PDHD and PDVD induction/collection planes, mm.
    const double deg = M_PI / 180.0;
    const double angles[] = {60 * deg, -60 * deg, 35.71 * deg, -35.71 * deg, 0.0};
    const double pitches[] = {3.0, 4.669, 4.792, 7.65, 5.10};
    std::mt19937 rng(101);
    std::uniform_real_distribution<double> pos(-3000.0, 3000.0);
    std::uniform_real_distribution<double> cen(-2000.0, 2000.0);
    int n = 0, bad = 0;
    for (double a : angles) {
        for (double pitch : pitches) {
            const double center = cen(rng);
            for (int i = 0; i < 2000; ++i) {
                const Facade::geo_point_t p(pos(rng), pos(rng), pos(rng));
                const double c = Facade::point2wind_cont(p, a, pitch, center);
                if (static_cast<int>(std::round(c)) != Facade::point2wind(p, a, pitch, center)) ++bad;
                ++n;
            }
            // a point on wire k's centre: cos(a)*z - sin(a)*y = pitch*(k+0.5) + center
            const int k = 17;
            const double s = pitch * (k + 0.5) + center;
            const Facade::geo_point_t pc(0.0, -std::sin(a) * s, std::cos(a) * s);
            CHECK(Facade::point2wind_cont(pc, a, pitch, center) == doctest::Approx(k).epsilon(1e-9));
        }
    }
    CHECK(n == 50000);
    CHECK(bad == 0);
}
