#include "WireCellUtil/CoordRegion.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"

#include <random>
#include <chrono>

using namespace WireCell;
using spdlog::debug;
using String::format;

TEST_CASE("coord region scalar")
{
    for (int axis=0; axis<3; ++axis) {

        const CoordBoundsScalar lo(-1, axis);
        const CoordBoundsScalar hi(+1, axis);

        CHECK(lo.location() == -1);
        CHECK(hi.location() == +1);

        Point origin, plo, phi, elo;

        plo[axis] = -2;
        phi[axis] = +2;
        elo[axis] = -1;

        CHECK(lo < origin);
        CHECK(origin < hi);
        CHECK(!(origin < lo));
        CHECK(!(hi < origin));

        CHECK(plo < lo);
        CHECK(plo < hi);
        CHECK(lo < phi);
        CHECK(hi < phi);

        // double dist = lo.distance(elo);
        // debug("axis {} point {} is distance {} from lo.  negative? {}",
        //       axis, elo, dist, dist<0);

        CHECK(!(elo<lo));
        CHECK(!(lo<elo));

        CoordRegion cr(lo, hi);

        CHECK(cr.inside(origin));
        CHECK(!cr.inside(plo));
        CHECK(!cr.inside(phi));
    }
}
static std::vector<double> uniform_values(size_t seed, size_t number = 3, 
                                          double vmin=-1.0, double vmax=1.0)
{
    std::default_random_engine re{seed};
    std::uniform_real_distribution<double> dist(vmin, vmax);

    std::vector<double> vals(number);
    for (size_t ind=0; ind<number; ++ind) {
        vals[ind] = dist(re);
    }
    return vals;
}

TEST_CASE("coord region sampled")
{
    const size_t npts = 100;

    for (int axis=0; axis<3; ++axis) {
        const size_t seed = npts * (axis+1);
        std::vector<double> x = uniform_values(seed+1, npts);
        std::vector<double> y = uniform_values(seed+2, npts);
        std::vector<double> z = uniform_values(seed+3, npts);

        // This makes an "ugly" "surface" which violates the assumptions of the
        // interpolation but we will avoid problems by not testing points in the
        // interior of the [-1,1]^3 cube of random points.
        const CoordBoundsSampled lo(x,y,z, axis);
        const CoordBoundsScalar hi(+2, axis);

        Point origin, pin, pout;
        pin[axis] = 1.5;
        pout[axis] = 2.5;

        CHECK(lo < pin);
        CHECK(lo < pout);
        CHECK(pin < hi);
        CHECK(hi < pout);
        
        CoordRegion cr(lo, hi);

        CHECK(cr.inside(pin));
        CHECK(!cr.inside(pout));
    }
}


static void dump_time(const std::chrono::steady_clock::time_point& t1,
                      const std::chrono::steady_clock::time_point& t2,
                      const std::string& msg="")
{
    double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    std::vector<std::string> units = {"ns","us","ms","s"};
    std::string found;
    
    for (const auto& unit : units) {
        if (dt < 1000) {
            found = unit;
            break;
        }
        dt /= 1000;
    }
    debug("{:.3f} {}\t{}", dt, found, msg);
}

TEST_CASE("coord region sampled speed")
{
    using clock = std::chrono::steady_clock;
    std::vector<size_t> npts_list = {100, 1000, 10000};
    std::vector<size_t> ndepos_list = {1000, 1000000};
    for (const auto& npts : npts_list) {

        const size_t seed = npts;
        std::vector<double> x = uniform_values(seed+1, npts);
        std::vector<double> y = uniform_values(seed+2, npts);
        std::vector<double> z = uniform_values(seed+3, npts);

        auto t0 = clock::now();
        const CoordBoundsSampled lo(x,y,z);
        const CoordBoundsScalar hi(+2);
        CoordRegion cr(lo, hi);
        auto t1 = clock::now();
        dump_time(t0,t1,format("make coord region for %d samples", npts));

        for (const auto& ndepos : ndepos_list) {
            const size_t dseed = seed + ndepos;
            std::vector<double> dx = uniform_values(dseed+1, ndepos);
            std::vector<double> dy = uniform_values(dseed+2, ndepos);
            std::vector<double> dz = uniform_values(dseed+3, ndepos);

            auto t2 = clock::now();
            for (size_t ind=0; ind<ndepos; ++ind) {
                cr.inside(Point(dx[ind], dy[ind], dz[ind]));
            }
            auto t3 = clock::now();
            dump_time(t2,t3,format("check %d depos against %d samples", ndepos, npts));
        }
    }
}
