// doc pdvd/44 -- TrackFitting::Parameters::gaus_nsigma, the acceptance window
// of cal_gaus_integral in sigmas.
//
// Before this round the window was the bare literal 4 at all nine call sites
// of cal_gaus_integral_seg.  It is now a runtime parameter (set_parameter
// "gaus_nsigma") whose default 4.0 reproduces the literal bit for bit.  These
// cases pin:
//   * the default is 4.0 -- changing it changes every production fit, so a
//     reverted or retyped default breaks case 1;
//   * the set/get round trip through the string dispatch -- a typo in either
//     branch breaks case 2 (an unknown name throws ValueError);
//   * the integral's dependence on nsigma: a wire bin farther than nsigma*sigma
//     from the centre receives nothing, a wider window admits it, and inside
//     the window the value does not depend on nsigma (no renormalisation).

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellClus/TrackFitting.h"

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("doc pdvd/44 gaus_nsigma: default is the legacy literal 4")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("gaus_nsigma") == 4.0);
    TrackFitting::Parameters p;
    CHECK(p.gaus_nsigma == 4.0);
}

TEST_CASE("doc pdvd/44 gaus_nsigma: set_parameter / get_parameter round trip")
{
    TrackFitting tf;
    tf.set_parameter("gaus_nsigma", 6.0);
    CHECK(tf.get_parameter("gaus_nsigma") == 6.0);
    CHECK_THROWS_AS(tf.set_parameter("gaus_nsigmas", 6.0), ValueError);
}

TEST_CASE("doc pdvd/44 gaus_nsigma: the window gates the wire bin, no renormalisation")
{
    TrackFitting tf;
    // centre at wire 0, sigma_T = 0.2 wire (PDVD induction under the pre-doc-44
    // constants); the first neighbour (wire 1) is 5 sigma away.
    const double t_center = 0.0, t_sigma = 1.0, w_center = 0.0, w_sigma = 0.2;
    double far4 = tf.cal_gaus_integral(0, 1, t_center, t_sigma, w_center, w_sigma, 0, 4.0, 1);
    double far6 = tf.cal_gaus_integral(0, 1, t_center, t_sigma, w_center, w_sigma, 0, 6.0, 1);
    CHECK(far4 == 0.0);
    CHECK(far6 > 0.0);
    // the centre bin is inside both windows and its value is independent of nsigma
    double c4 = tf.cal_gaus_integral(0, 0, t_center, t_sigma, w_center, w_sigma, 0, 4.0, 1);
    double c6 = tf.cal_gaus_integral(0, 0, t_center, t_sigma, w_center, w_sigma, 0, 6.0, 1);
    CHECK(c4 > 0.0);
    CHECK(c4 == c6);
    // and the literal-4 legacy equals the double default
    double lit = tf.cal_gaus_integral(0, 1, t_center, t_sigma, w_center, 0.3, 0, 4, 1);
    double dbl = tf.cal_gaus_integral(0, 1, t_center, t_sigma, w_center, 0.3, 0, tf.get_parameter("gaus_nsigma"), 1);
    CHECK(lit == dbl);
    CHECK(lit > 0.0);
}
