#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/RaySolving.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Logging.h"

#include "TCanvas.h"
#include "TMarker.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TArrow.h"
#include "TH1F.h"

#include <boost/graph/graphviz.hpp>

#include <math.h>

#include <random>

using namespace WireCell;
using namespace WireCell::Waveform;
using namespace WireCell::RayGrid;
using namespace std;
using spdlog::info;
using spdlog::warn;

const int ndepos = 10;
const int neles = 10;
const double pitch_magnitude = 5;
const double gaussian = 3;
const double border = 10;
const double width = 100;
const double height = 100;

// local helper codes
#include "raygrid.h"
#include "raygrid_draw.h"

// Invent a channel and wire map.  Note, in this test, this mapping is
// used identically for both "channels" and "wires" and it is wholly
// innappropraite for real detectors.  Do not copy-paste.
Grouping::ident_t make_ident(int index, int layer, int face = 0)
{
    return (1 + face) * 10000 + (1 + layer) * 1000 + index;
}

int main(int argc, char* argv[])
{
    Printer print(argv[0]);

    auto raypairs = make_raypairs(width, height, pitch_magnitude);
    const int nlayers = raypairs.size();

    Coordinates coords(raypairs);

    draw_raygrid(print, coords, raypairs);

    const auto& pitches = coords.pitch_dirs();
    const auto& centers = coords.centers();
    const auto& pitch_mags = coords.pitch_mags();

    double nudge=0;
    if (argc > 1) {
        nudge = atof(argv[1]);
    }
    Tiling tiling(coords, nudge);

    std::vector<std::vector<Activity::value_t> > measures(nlayers);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> position(0, std::max(width, height));
    std::normal_distribution<double> spread(0.0, gaussian);
    std::vector<Point> points;
    for (int idepo = 0; idepo < ndepos; ++idepo) {
        Point cp(0, position(generator), position(generator));
        for (int iele = 0; iele < neles; ++iele) {
            Point delta(0, spread(generator), spread(generator));
            points.push_back(cp + delta);
        }
    }

    draw_frame(print.canvas, "Points and Activity");
    for (size_t ipt = 0; ipt < points.size(); ++ipt) {
        const auto& p = points[ipt];
        draw_point(p, 1, 24, ipt + 1);
        for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
            // pimpos normally would help here to find pitch location of arb point
            const auto& pit = pitches[ilayer];
            const auto& cen = centers[ilayer];
            const auto rel = p - cen;
            const int pit_ind = pit.dot(rel) / pitch_mags[ilayer];
            auto& m = measures[ilayer];
            if ((int) m.size() <= pit_ind) {
                m.resize(pit_ind + 1, 0.0);
            }

            m[pit_ind] += 1.0;
        }
    }
    for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
        draw_layer(coords, ilayer, pitch_mags[ilayer], pitches[ilayer], centers[ilayer], measures[ilayer]);
    }
    print();

    blobs_t blobs;

    draw_frame(print.canvas, "Points and Strips");
    for (size_t ipt = 0; ipt < points.size(); ++ipt) {
        const auto& p = points[ipt];
        draw_point(p, 1, 24, ipt + 1);
    }
    activities_t activities;
    for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
        auto& m = measures[ilayer];

        for (size_t ind = 0; ind < m.size(); ++ind) {
            if (m[ind] <= 0.0) {
                continue;
            }
            info("L{} [{}] {}", ilayer, ind, m[ind]);
        }

        Activity activity(ilayer, {m.begin(), m.end()});

        auto strips = activity.make_strips();
        draw_strips(coords, strips);
        activities.push_back(activity);

        auto tot = std::accumulate(m.begin(), m.end(), 0.0);
        info("L{} activity={} in: {} strips", activity.layer(), tot, strips.size());
    }
    print();

    for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
        const auto& activity = activities[ilayer];
        info("Tiling layer {} with {} blobs: {}", ilayer, blobs.size(), activity.as_string());
        if (blobs.empty()) {
            blobs = tiling(activity);
        }
        else {
            blobs = tiling(blobs, activity);
            if (blobs.empty()) {
                warn("lost m'blobs at layer index {} of {}!", ilayer, nlayers);
                return -1;
            }
        }
        drop_invalid(blobs);
        dump(blobs);
        draw_points_blobs(coords, print, points, blobs);
        print();
    }

    draw_points_blobs(coords, print, points, blobs);
    for (const auto& activity : activities) {
        auto strips = activity.make_strips();
        draw_strips(coords, strips, false);
    }
    print();

    // now test solving
    Solving solving;

    // Skip the first two horiz/vert layers bounding sensitivity.  In
    // principle, they can be included but just add fully degenerate
    // terms.
    for (int ilayer = 2; ilayer < nlayers; ++ilayer) {
        Grouping grouping;
        auto& m = measures[ilayer];

        // Load up the "channels" and their associated "wires".  In
        // this test we have a simple one-to-one mapping between both
        // based on a common index.
        for (size_t mind = 0; mind < m.size(); ++mind) {
            const float meas = m[mind];
            if (meas <= 0) {
                continue;
            }
            Grouping::ident_t ident = make_ident(mind, ilayer);
            info("ilayer:{} mind:{} indent:{} meas:{}", ilayer, mind, ident, meas);
            grouping.add('m', ident, {ident}, meas);
        }

        // Load up the "blobs" and their associated "wires" for the current layer
        for (size_t bind = 0; bind < blobs.size(); ++bind) {
            const auto& blob = blobs[bind];
            std::vector<Grouping::ident_t> wids;
            for (const auto& strip : blob.strips()) {
                if (strip.layer != ilayer) {
                    continue;
                }
                for (auto wind = strip.bounds.first; wind < strip.bounds.second; ++wind) {
                    wids.push_back(make_ident(wind, ilayer));
                }
            }
            // note, blob ident must NOT inculde ilayer information
            const float blob_value = 0.0;
            const float blob_weight = 1.0;
            grouping.add('s', bind, wids, blob_value, blob_weight);
        }

        const auto& g = grouping.graph();
        auto labels = [&](std::ostream& out, const auto& vtx) {
            char typ = g[vtx].ntype;
            int lid = g[vtx].ident % 1000;  // l'il ID, erase face and layer
            if (typ == 's') {
                out << "[label=\"s" << lid << "\"]";
            }
            if (typ == 'm') {
                out << "[label=\"m" << lid << "=" << g[vtx].value << "\"]";
            }
            if (typ == 'w') {
                out << "[label=\"w" << lid << "\"]";
            }
        };

        std::string dotfilename = Form("%s-layer%d.dot", argv[0], ilayer);
        std::ofstream dotfile(dotfilename.c_str());
        boost::write_graphviz(dotfile, g, labels);
        cerr << dotfilename << "\n";

        auto clusters = grouping.clusters();
        solving.add(clusters);
    }

    auto solution = solving.solve();

    for (const auto& it : solution) {
        std::cerr << it.first << ": " << it.second << std::endl;
    }

    const auto& g = solving.graph();
    auto labels = [&](std::ostream& out, const auto& vtx) {
        char typ = g[vtx].ntype;
        if (typ == 's') {
            out << "[label=\"s" << g[vtx].ident << "=" << Form("%.1f", g[vtx].value) << "\"]";
        }
        if (typ == 'm') {
            out << "[label=\"m"
                << "=" << g[vtx].value << "\"]";
        }
    };

    std::string dotfilename = Form("%s-solution.dot", argv[0]);
    std::ofstream dotfile(dotfilename.c_str());
    boost::write_graphviz(dotfile, g, labels);
    cerr << dotfilename << "\n";

    // plot something

    draw_points_blobs_solved(coords, print, points, blobs, solution);
    print();

    return 0;
}
