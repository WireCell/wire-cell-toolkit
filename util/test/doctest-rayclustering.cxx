#include "WireCellUtil/RayClustering.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/RayHelpers.h"

#include <math.h>

#include <random>
#include <fstream>
#include <string>

using namespace WireCell;
using namespace WireCell::Waveform;
using namespace WireCell::RayGrid;
using namespace std;
using spdlog::debug;
using spdlog::info;
using spdlog::warn;

const int ndepos = 10;
const int neles = 10;
const double pitch_magnitude = 5;
const double gaussian_spread = 3;
const double border = 10;
const double width = 100;
const double height = 100;

#include "raygrid_dump.h"

static std::vector<Point> make_points(std::default_random_engine& generator, double x)
{
    std::vector<Point> points;
    std::uniform_real_distribution<double> position(0, std::max(width, height));
    std::normal_distribution<double> spread(0.0, gaussian_spread);
    for (int idepo = 0; idepo < ndepos; ++idepo) {
        Point cp(0, position(generator), position(generator));
        for (int iele = 0; iele < neles; ++iele) {
            const Point delta(x, spread(generator), spread(generator));
            const Point pt = cp + delta;
            if (pt.y() < -border or pt.y() > height + border or pt.z() < -border or pt.z() > width + border) {
                warn("Rejecting far away point: {} + {}", cp, delta);
                continue;
            }
            points.push_back(cp + delta);
        }
    }
    return points;
}

typedef std::vector<Activity::value_t> measure_t;

struct Chirp {
    const blobs_t& one;
    const blobs_t& two;
    Coordinates& coords;

    typedef typename std::unordered_set<std::size_t> indices_t;
    indices_t* sel1;
    indices_t* sel2;

    Chirp(const blobs_t& one, const blobs_t& two, Coordinates& coords, JsonEvent& dumper)
      : one(one)
      , two(two)
      , coords(coords)
      , sel1(new indices_t)
      , sel2(new indices_t)
    {
    }

    bool in(const blobref_t& a, const blobref_t& b)
    {
        if (surrounding(a, b)) {
            return true;
        }

        const auto& astrips = a->strips();
        const int nlayers = astrips.size();

        for (const auto& c : b->corners()) {
            int found = 0;
            for (layer_index_t layer = 0; layer < nlayers; ++layer) {
                const auto& astrip = astrips[layer];
                if (layer == c.first.layer) {
                    info("L{} A: {} {},{}", layer, astrip, c.first, c.second);
                    if (astrip.on(c.first.grid)) {
                        info("\ton with found={} nlayers={}", found, nlayers);
                        ++found;
                        continue;
                    }
                    info("\toff with found={} nlayers={}", found, nlayers);
                    break;
                }
                if (layer == c.second.layer) {
                    info("L{} A: {} {},{}", layer, astrip, c.first, c.second);
                    if (astrip.on(c.second.grid)) {
                        info("\ton with found={} nlayers={}", found, nlayers);
                        ++found;
                        continue;
                    }
                    info("\toff with found={} nlayers={}", found, nlayers);
                    break;
                }
                const double ploc = coords.pitch_location(c.first, c.second, layer);
                const int pind = coords.pitch_index(ploc, layer);

                info("L{} A: {} pind={} ploc={} {},{}", layer, astrip, pind, ploc, c.first, c.second);

                if (astrip.in(pind)) {
                    info("\tin with found={} nlayers={}", found, nlayers);
                    ++found;
                }
                else {
                    info("\tout with found={} nlayers={}", found, nlayers);
                    break;
                }
            }
            if (found == nlayers) {
                return true;
            }
        }
        return false;
    }

    void operator()(const blobref_t& a, const blobref_t& b)
    {
        const std::size_t d1 = a - one.begin();
        const std::size_t d2 = b - two.begin();

        info("overlap: a{} and b{}", d1, d2);
        info("\tblob a #{}: {}", d1, a->as_string());
        info("\tblob b #{}: {}", d2, b->as_string());

        if (!this->in(a, b)) {
            warn("NO CONTAINED CORNERS");
        }

        sel1->insert(d1);
        sel2->insert(d2);
    }

    void dump(JsonEvent& dumper)
    {
        for (const auto ind : *sel1) {
            const auto& br = one[ind];
            dumper(br, 10.0, 1.0, 1, ind);
        }
        for (const auto ind : *sel2) {
            const auto& br = two[ind];
            dumper(br, 20.0, 1.0, 2, ind);
        }
    }
};

static void test_blobs(const blobs_t& blobs)
{
    for (const auto& blob : blobs) {
        const auto& strips = blob.strips();
        debug("blob: {}", blob.as_string());
        for (size_t ind = 0; ind < 5; ++ind) {
            debug("bb[{}]: {} {}", ind, strips[ind].bounds.first, strips[ind].bounds.second);
        }
        CHECK(strips[0].bounds.first == 0);
        CHECK(strips[0].bounds.second == 1);
        CHECK(strips[1].bounds.first == 0);
        CHECK(strips[1].bounds.second == 1);
    }
}

TEST_SUITE("rayclustering") {

TEST_CASE("coordinates") {
    auto raypairs = symmetric_raypairs(width, height, pitch_magnitude);
    Coordinates coords(raypairs);
    CHECK(coords.nlayers() == 5);

    {
        Coordinates empty;
    }
    {
        Coordinates copy;
        copy = coords;
    }
    {
        Coordinates copy(coords);
    }
}

TEST_CASE("blob clustering") {
    auto raypairs = symmetric_raypairs(width, height, pitch_magnitude);
    Coordinates coords(raypairs);
    Tiling tiling(coords, 1e-6);
    (void)tiling;

    std::default_random_engine generator;
    std::vector<Point> pts1 = make_points(generator, 10.0);
    std::vector<Point> pts2 = make_points(generator, 20.0);

    std::vector<measure_t> meas1 = make_measures(coords, pts1);
    std::vector<measure_t> meas2 = make_measures(coords, pts2);

    auto act1 = make_activities(coords, meas1);
    auto act2 = make_activities(coords, meas2);

    auto blobs1 = make_blobs(coords, act1);
    auto blobs2 = make_blobs(coords, act2);

    test_blobs(blobs1);
    test_blobs(blobs2);

    JsonEvent dumper(coords);
    for (const auto& pt : pts1) {
        dumper(pt);
    }
    for (const auto& pt : pts2) {
        dumper(pt);
    }
    Chirp chirp(blobs1, blobs2, coords, dumper);
    associator_t chirpf = chirp;
    associate(blobs1, blobs2, chirpf);

    chirp.dump(dumper);
    // file output (dumper.dump to argv[0]+".json") removed
}

}  // TEST_SUITE("rayclustering")
