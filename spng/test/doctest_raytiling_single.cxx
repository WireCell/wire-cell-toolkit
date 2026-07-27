#include "WireCellSpng/RayTiling.h"
#include "WireCellSpng/RayGrid.h" // For WireCell::SPNG::RayGrid::Coordinates
#include "WireCellSpng/RayTest.h"

#include "WireCellUtil/Units.h"

// Name collission for "CHECK" between torch and doctest.
#undef CHECK
#include "WireCellUtil/doctest.h"

using namespace WireCell::SPNG::RayGrid;
using namespace WireCell;

TEST_CASE("spng ray grid tiling one point one blob") {

    const double width = 4000; //*units::mm;
    const double height = 4000; // *units::mm;
    const double pitch_magnitude = 5; //*units::mm;
    auto views = symmetric_views(width, height, pitch_magnitude);
    Coordinates coords(views);

    torch::Tensor points = torch::tensor({ {2000.0, 2000.0} });
    auto activities = fill_activity(coords, points);

    // Build the blobs long hand

    auto blobs = trivial_blobs();
    for (size_t view = 2; view < 5; ++view) {
        auto activity = activities[view];
        blobs = apply_activity(coords, blobs, activity);
        std::cerr << "view " << view << " has " << blobs.size(0) << " blobs\n";
    }

}
