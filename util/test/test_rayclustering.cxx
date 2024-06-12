#include "WireCellUtil/RayClustering.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Testing.h"
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
const double gaussian = 3;
const double border = 10;
const double width = 100;
const double height = 100;

#include "raygrid_dump.h"


static std::vector<Point> make_points(std::default_random_engine& generator, double x)
{
    std::vector<Point> points;
    std::uniform_real_distribution<double> position(0, std::max(width, height));
    std::normal_distribution<double> spread(0.0, gaussian);
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

// moved to RayHelpers
// static std::vector<measure_t> make_measures(Coordinates& coords, const std::vector<Point>& points)
// {
//     int nlayers = coords.nlayers();
//     std::vector<measure_t> measures(nlayers);
//     const auto& pitches = coords.pitch_dirs();
//     const auto& centers = coords.centers();
//     const auto& pitch_mags = coords.pitch_mags();

//     for (size_t ipt = 0; ipt < points.size(); ++ipt) {
//         const auto& p = points[ipt];
//         for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
//             const auto& pit = pitches[ilayer];
//             const auto& cen = centers[ilayer];
//             const auto rel = p - cen;
//             const int pit_ind = pit.dot(rel) / pitch_mags[ilayer];
//             if (pit_ind < 0) {
//                 warn("Negative pitch indices not allowed, got {} from ilayer {} ipt {} for point {}", pit_ind, ilayer,
//                      ipt, p);
//                 continue;
//             }
//             if (ilayer <= 1) {
//                 if (pit_ind >= 1 or pit_ind < 0) {
//                     debug("mm: pit_ind={} with ipt={}", pit_ind, ipt);
//                     if (pit_ind == 1) {
//                         debug("\tpit={} cen={} rel={}", pit, cen, rel);
//                     }
//                     continue;
//                 }
//             }
//             measure_t& m = measures[ilayer];
//             if ((int) m.size() <= pit_ind) {
//                 debug("resize for ipt {} ilayer {} from {} to {}", ipt, ilayer, m.size(), pit_ind + 1);
//                 m.resize(pit_ind + 1, 0.0);
//                 debug("done");
//             }

//             debug("adding to pit_ind {} ilayer {} ipt {}", pit_ind, ilayer, ipt);
//             m[pit_ind] += 1.0;
//             debug("valud: {}", m[pit_ind]);
//         }
//     }

//     return measures;
// }

// static activities_t make_activities(Coordinates& coords, std::vector<measure_t>& measures)
// {
//     int nlayers = coords.nlayers();
//     activities_t activities;
//     for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
//         auto& m = measures[ilayer];
//         info("Make activity for layer: {}: {}", ilayer, m.size());
//         std::stringstream ss;
//         for (const auto& v : m) {
//             ss << v << " ";
//         }
//         info("\t{}", ss.str());

//         Activity activity(ilayer, {m.begin(), m.end()});
//         Assert(!activity.empty());
//         activities.push_back(activity);
//     }
//     return activities;
// }

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

        // if at least one corner of b is in side a, return true
        for (const auto& c : b->corners()) {
            int found = 0;
            // go through each layer of blob a
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
            // Assert(this->in(a,b));
        }

        sel1->insert(d1);
        sel2->insert(d2);
    }

    void dump(JsonEvent& dumper)
    {
        int number = 0;
        for (const auto ind : *sel1) {
            const auto& br = one[ind];
            dumper(br, 10.0, 1.0, 1, ind);
            ++number;
        }
        number = 0;
        for (const auto ind : *sel2) {
            const auto& br = two[ind];
            dumper(br, 20.0, 1.0, 2, ind);
            ++number;
        }

        return;
    }
};

static void test_blobs(const blobs_t& blobs)
{
    for (const auto& blob : blobs) {
        const auto& strips = blob.strips();
        std::cerr << "blob: " << blob << "\n";
        for (size_t ind=0; ind<5; ++ind) {
            std::cerr << "\tbb["<<ind<<"]: " << strips[ind].bounds << "\n";
        }

        Assert(strips[0].bounds.first == 0);
        Assert(strips[0].bounds.second == 1);
        Assert(strips[1].bounds.first == 0);
        Assert(strips[1].bounds.second == 1);
    }
}


int main(int argc, char* argv[])
{
    auto raypairs = symmetric_raypairs(width, height, pitch_magnitude);

    Coordinates coords(raypairs);
    Assert(coords.nlayers() == 5);

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

    Tiling tiling(coords, 1e-6);

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

    // fixme: blobs are missing these features:
    // - sort corners
    // - apply X offset
    // - assign some value
    // this is maybe best done converting from Blob to another rep

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

    std::string fout = argv[0];
    fout += ".json";
    dumper.dump(fout);

    return 0;
}
