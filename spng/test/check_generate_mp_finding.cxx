// This runs tiling on random points and dumps to a torch pickle file.

#include "WireCellSpng/RayTiling.h"
#include "WireCellSpng/RayTest.h"
#include "WireCellSpng/Stopwatch.h"
#include <fstream>

using namespace WireCell::SPNG::RayGrid;

const double pitch_magnitude = 5;
const double gaussian = 3;


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

torch::Tensor make_indices(torch::Tensor frame) {
    auto reshaped = frame.reshape({-1, frame.size(-1)});
    auto counts = reshaped.sum(-1);
    auto max_size = torch::max(counts).item<long>();
    auto results = torch::full({reshaped.size(0), max_size}, -1);

    auto col_indices = reshaped.nonzero().index({"...", -1});
    int64_t current_offset = 0;
    for (int64_t i = 0; i < reshaped.size(0); ++i) {
        auto count = counts.index({i}).item<int64_t>();
        // std::cout << i << " " << count << std::endl;
        if (count > 0) {
            // Get the indices for the current row
            auto current_col_indices = col_indices.slice(0, current_offset, current_offset + count);

            // Fill the corresponding slice in the result tensor
            results.index_put_({i, torch::indexing::Slice(0, count)}, current_col_indices.to(torch::kInt));
        }
        current_offset += count;
    }

    return results;
}

int main(int argc, char* argv[])
{
    // usage: main output
    std::string output = "check_generate_tiling.pt";
    int ngroups = 10;
    int nsamps = 1;
    int nframes = 1;
    double width = 55;
    double height = 55;
    int seed = 42;
    for (int i = 0; i < argc; ++i) {
        std::cout << i << " " << argv[i] << std::endl;
        if (strcmp(argv[i], "-o") == 0) {
            output = argv[i+1];
        }
        else if (strcmp(argv[i], "--ngroups") == 0) {
            ngroups = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--nsamps") == 0) {
            nsamps = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--nframes") == 0) {
            nframes = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--width") == 0) {
            width = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--height") == 0) {
            height = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--seed") == 0) {
            seed = atoi(argv[++i]);
        }
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
    // auto points = random_groups(ngroups, 10, gaussian, {0.0, 0.0}, {width, height});
    
    std::vector<torch::Tensor> all_frames_2, all_frames_3, all_frames_4;

    for (int iframe = 0; iframe < nframes; ++iframe) {

        std::vector<torch::Tensor> all_points(nsamps);
        std::vector<torch::Tensor> all_activities_2, all_activities_3, all_activities_4;
        for (auto & points : all_points) {
            points = random_groups(ngroups, 10, gaussian, {0., 0.}, {width, height});
            auto activities = fill_activity(coords, points);
            all_activities_2.push_back(activities[2].unsqueeze(0));
            all_activities_3.push_back(activities[3].unsqueeze(0));
            all_activities_4.push_back(activities[4].unsqueeze(0));
            
        }
        auto frame_2 = torch::cat(all_activities_2),
        frame_3 = torch::cat(all_activities_3),
        frame_4 = torch::cat(all_activities_4);
        
        all_frames_2.push_back(frame_2.unsqueeze(0));
        all_frames_3.push_back(frame_3.unsqueeze(0));
        all_frames_4.push_back(frame_4.unsqueeze(0));

        std::cout << frame_2 << std::endl;
        std::cout << frame_2.sum(-1) << std::endl;
        std::cout << torch::max(frame_2.sum(-1)) << std::endl;
    }
    auto mega_frame_2 = torch::cat(all_frames_2);
    auto mega_frame_3 = torch::cat(all_frames_3);
    auto mega_frame_4 = torch::cat(all_frames_4);
    std::cout << mega_frame_2 << std::endl;
    std::cout << mega_frame_2.sum(-1) << std::endl;
    std::cout << torch::max(mega_frame_2.sum(-1)) << std::endl;

    // auto mega_frame_2_reshaped = mega_frame_2.reshape({-1, mega_frame_2.size(-1)});
    // std::cout << "reshaped\n" << mega_frame_2_reshaped << std::endl;

    // // auto sizes = mega_frame_2.sizes().vec();
    // // sizes.pop_back();
    // auto frame_2_counts = mega_frame_2_reshaped.sum(-1);
    // auto max_size = torch::max(frame_2_counts).item<long>();
    // // sizes.push_back(max_size);
    // auto indices_2 = torch::full({mega_frame_2_reshaped.size(0), max_size}, -1);
    
    // std::cout << indices_2 << std::endl;
    // // std::cout << indices_2.reshape({-1, max_size});

    // auto col_indices = mega_frame_2_reshaped.nonzero().index({"...", -1});
    // std::cout << col_indices << std::endl;
    // // Use a for loop to fill the tensor
    // // This loop is necessary to handle the ragged nature of the data
    // int64_t current_offset = 0;
    // for (int64_t i = 0; i < mega_frame_2_reshaped.size(0); ++i) {
    //     auto count = frame_2_counts.index({i}).item<int64_t>();
    //     std::cout << i << " " << count << std::endl;
    //     if (count > 0) {
    //         // Get the indices for the current row
    //         auto current_col_indices = col_indices.slice(0, current_offset, current_offset + count);

    //         // Fill the corresponding slice in the result tensor
    //         indices_2.index_put_({i, torch::indexing::Slice(0, count)}, current_col_indices.to(torch::kInt));
    //     }
    //     current_offset += count;
    // }

    auto indices_2 = make_indices(mega_frame_2);
    std::cout << "Frame 2:\n" << indices_2 << std::endl;
    
    auto indices_3 = make_indices(mega_frame_3);
    std::cout << "Frame 3:\n" << indices_3 << std::endl;
    
    //Make Nl x Nm pairs of indices of the active wires (rays) in each plane
    auto cross = torch::zeros({indices_2.size(0), indices_2.size(-1), indices_3.size(-1), 2}, torch::TensorOptions().dtype(indices_2.dtype()));
    // cross.index_put_({Slice(), Slice(), 1}, m_hi);
    cross.index_put_({"...", 1}, indices_3.unsqueeze(1).expand({-1, indices_2.size(-1), -1}));
    cross = cross.transpose(1,2);
    cross.index_put_({"...", 0}, indices_2.unsqueeze(1).expand({-1, indices_3.size(-1), -1}));
    cross = cross.reshape({-1, 2});

    std::cout << cross << std::endl;

    //Get the ray indices for each plane/view
    auto r1 = cross.index({"...", 0});
    auto r2 = cross.index({"...", 1});

    // auto indices_4 = make_indices(mega_frame_4);
    std::cout << "Frame 4:\n" << mega_frame_4 << std::endl;
    // auto flattened_4 = indices_4.flatten();
    auto view2 = torch::full_like(r1, 2, torch::TensorOptions().dtype(r1.dtype()));
    auto view3 = torch::full_like(r1, 3, torch::TensorOptions().dtype(r1.dtype()));
    auto view4 = torch::full_like(r1, 4, torch::TensorOptions().dtype(r1.dtype()));
    auto view4_short = torch::tensor(4, torch::TensorOptions().dtype(r1.dtype()));

    auto locs = coords.pitch_location(view2, r1, view3, r2, view4);
    //Convert the locations to pitch indices
    auto results = coords.pitch_index(
        locs,
        view4_short
    );
    results = results.reshape({-1, indices_2.size(-1)*indices_3.size(-1)});
    std::cout << "Results\n" << results << std::endl;
    // auto where = torch::where((results >= 0) & (results < std::round(width/pitch_magnitude)));
    // for (auto & w : where) {w = w.unsqueeze(0);}
    // std::cout << "Where\n" << where << std::endl;
    // std::cout << torch::cat(where) << std::endl;
    auto valid_range = (results >= 0) & (results < std::round(width/pitch_magnitude));
    std::cout << valid_range << std::endl;
    std::cout << results.index({valid_range}) << std::endl;

    // auto valid_counts = valid_range.sum(-1);
    // auto valid_max_size = torch::max(valid_counts).item<long>();

    
    auto mega_frame_4_reshaped = mega_frame_4.reshape({-1, mega_frame_4.size(-1)});
    auto mp2 = torch::zeros_like(mega_frame_4_reshaped);
    auto mp3 = torch::zeros_like(mega_frame_4_reshaped);
    // long current_offset = 0;
    for (long i = 0; i < mega_frame_4_reshaped.size(0); ++i) {
        
        auto these_indices = results.index({i, valid_range.index({i})});
        // auto target_row = mega_frame_4_reshaped.index({i, indices});
        mp3.index_put_({i, these_indices}, mega_frame_4_reshaped.index({i, these_indices}));
        mp2.index_put_({i, these_indices}, (mega_frame_4_reshaped.index({i, these_indices}) == 0));
        // std::cout << target_row << std::endl;
        
        // // Get the indices for the current row
        // auto current_col_indices = col_indices.slice(0, current_offset, current_offset + count);
        // // Fill the corresponding slice in the result tensor
        // results.index_put_({i, torch::indexing::Slice(0, count)}, current_col_indices.to(torch::kInt));
        // current_offset += count;
    }
    std::cout << mega_frame_4_reshaped << std::endl;
    std::cout << mp2.reshape({nframes, nsamps, -1}) << std::endl;
    std::cout << mp3.reshape({nframes, nsamps, -1}) << std::endl;

    // std::cout << results.index(torch::where((results >= 0) & (results < std::round(width/pitch_magnitude)))) << std::endl;
    // std::vector<int> results_counts;
    // long index = 0;
    // for (long i = 0; i < results.size(0); ++i) {
        
    // }

    // std::cout << results.index(
    //             {torch::where((results >= 0) & (results < std::round(width/pitch_magnitude)))[0]}
    //         )
    //             << std::endl;

    // to_save.insert("points", points);

    // auto activities = fill_activity(coords, points);
    // assert(activities.size() == 5);
    // {
    //     // There may be some way to store a list-of-tensor in the tensor map but
    //     // gemini and I can't get it to work, so punt and enumerate.
    //     for (size_t layer = 0; layer<5; ++layer) {
    //         to_save.insert("activity_" + std::to_string(layer), activities[layer]);
    //     }
    // }
    
    // auto blobs = trivial_blobs();
    // for (size_t view = 2; view < 5; ++view) {
    //     auto activity = activities[view];

    //     blobs = apply_activity(coords, blobs, activity);
    // }
    // to_save.insert("blobs", blobs);


    // torch::Tensor crossings = blob_crossings(blobs);
    // to_save.insert("crossings", crossings);

    // torch::Tensor insides = blob_insides(coords, blobs, crossings);
    // to_save.insert("insides", insides);


    // std::cerr << "writing " << output << "\n";
    // std::ofstream output_file(output, std::ios::binary);
    // auto data = torch::pickle_save(to_save);
    // output_file.write(data.data(), data.size());

    return 0;
}
