// This runs tiling on random points and dumps to a torch pickle file.

#include "WireCellSpng/RayTiling.h"
#include "WireCellSpng/RayTest.h"
#include "WireCellSpng/Stopwatch.h"
#include <fstream>

using namespace WireCell::SPNG::RayGrid;

const double pitch_magnitude = 5;
const double gaussian = 3;
const double width = 4000;
const double height = 4000;

using tensor_map = torch::Dict<std::string, torch::Tensor>;

static void store_tensors(Coordinates& coords, tensor_map& store)
{
    store.insert("coords_pitch_mag", coords.pitch_mag);
    store.insert("coords_pitch_dir", coords.pitch_dir);
    store.insert("coords_center", coords.center);
    store.insert("coords_zero_crossings", coords.zero_crossings);
    store.insert("coords_ray_jump", coords.ray_jump);
    store.insert("coords_a", coords.a);
    store.insert("coords_b", coords.b);
    store.insert("coords_ray_dir", coords.ray_dir);
}
                          

int main(int argc, char* argv[])
{
    // usage: main output
    std::string output = "check_generate_tiling.pt";
    if (argc > 1) {
        output= argv[1];
    }
    int ngroups = 1000;
    if (argc > 2) {
        ngroups = atoi(argv[2]);
    }
    int seed = 42;
    if (argc > 3) {
        seed = atoi(argv[3]);
    }

    torch::AutoGradMode enable_grad(false);

    torch::manual_seed(seed);

    tensor_map to_save;

    auto views = symmetric_views(width, height, pitch_magnitude);
    assert(views.size(0) == 5);
    to_save.insert("views", views);

    Coordinates coords(views);
    store_tensors(coords, to_save);

    std::cerr << "Generating " << ngroups << " point groups\n";
    auto points = random_groups(ngroups, 10, gaussian, {0.0, 0.0}, {width, height});
    to_save.insert("points", points);

    auto activities = fill_activity(coords, points);
    assert(activities.size() == 5);
    {
        // There may be some way to store a list-of-tensor in the tensor map but
        // gemini and I can't get it to work, so punt and enumerate.
        for (size_t layer = 0; layer<5; ++layer) {
            to_save.insert("activity_" + std::to_string(layer), activities[layer]);
        }
    }
    
    auto blobs = trivial_blobs();
    for (size_t view = 2; view < 5; ++view) {
        auto activity = activities[view];

        blobs = apply_activity(coords, blobs, activity);
    }
    to_save.insert("blobs", blobs);


    torch::Tensor crossings = blob_crossings(blobs);
    to_save.insert("crossings", crossings);

    torch::Tensor insides = blob_insides(coords, blobs, crossings);
    to_save.insert("insides", insides);


    std::cerr << "writing " << output << "\n";
    std::ofstream output_file(output, std::ios::binary);
    auto data = torch::pickle_save(to_save);
    output_file.write(data.data(), data.size());

    return 0;
}
