#include "WireCellSpng/RayTiling.h"
#include "WireCellSpng/RayTest.h"
#include "WireCellSpng/Stopwatch.h"
#include <fstream>

using namespace WireCell::Spng::RayGrid;

const double pitch_magnitude = 5;
const double gaussian = 3;
const double width = 4000;
const double height = 4000;


torch::Device parse_device(const std::string& device_name)
{
    torch::Device device(torch::kCPU);
    if (device_name == "cpu") {
        return torch::Device(torch::kCPU);
    }
    if (device_name == "gpu" || device_name == "cuda") {
        if (torch::cuda::is_available()) {
            std::cerr << "CUDA is available for device: "<<device_name<<"\n";
            return torch::Device(torch::kCUDA);
        }
        else {
            std::cerr << "CUDA is NOT available for device: "<<device_name<<"\n";
            throw std::runtime_error("No GPU available");
        }
    }
    std::cerr << "Unknown device: "<<device_name<<"\n";
    throw std::runtime_error("unknown device");
}

int main(int argc, char* argv[])
{
    std::string device_name = "cpu";
    if (argc > 1) {
        device_name = argv[1];
    }
    auto device = parse_device(device_name);


    std::string autograd_name = "";
    if (argc > 2) {
        autograd_name = argv[2];
    }
    bool use_autograd = false;
    if (autograd_name == "yes" || autograd_name == "true") {
        use_autograd = true;
    }
    torch::AutoGradMode enable_grad(use_autograd);

    std::string random_name = "";
    if (argc > 3) {
        random_name = argv[3];
    }
    if (random_name != "" ) {
        int64_t random_seed = atol(random_name.c_str());
        torch::manual_seed(random_seed);
    }

    int nblobs = 1000;
    if (argc > 4) {
        nblobs = std::atoi(argv[4]);
    }

    std::cerr << "Device=" << device << " autograd='" << use_autograd << "', random seed='" << random_name << "'\n";

    Stopwatch sw;

    auto views = symmetric_views(width, height, pitch_magnitude);
    assert(views.size(0) == 5);
    
    std::cerr << "Made symmetric views in " << sw.restart() << " us\n";;
    std::cout << views << std::endl;

    {
        auto bytes = torch::pickle_save(views);
        std::ofstream fout("views.zip", std::ios::out | std::ios::binary);
        fout.write(bytes.data(), bytes.size());
        fout.close();
    }

    Coordinates coords(views);
    std::cerr << "Made coordinates " << sw.restart() << " us\n";


    coords.to(device);
    std::cerr << "Moved coordinates to device " << sw.restart() << " us\n";

    auto points = random_groups(nblobs, 10, gaussian, {0.0, 0.0}, {width, height});
    std::cerr << "Made points in " << sw.restart() << " us, points.shape=" << points.sizes() << "\n";

    //   auto x = torch::ones({3, 3});
    auto bytes = torch::pickle_save(points);
    std::ofstream fout("random_points.zip", std::ios::out | std::ios::binary);
    fout.write(bytes.data(), bytes.size());
    fout.close();
    // torch::save(points, "random_points.pt");

    points = points.to(device);
    std::cerr << "Moved points to device " << sw.restart() << " us\n";

    auto activities = fill_activity(coords, points);
    assert(activities.size() == 5);
    std::cerr << "Made activities in " << sw.restart() << " us\n";
    
    for (size_t iview = 0; iview < activities.size(); ++iview) {
        auto bytes = torch::pickle_save(activities[iview]);
        std::string name = "activities" + std::to_string(iview) + ".zip"; 
        std::ofstream fout(name, std::ios::out | std::ios::binary);
        fout.write(bytes.data(), bytes.size());
        fout.close();
    }

    // activities[2].index_put_({torch::indexing::Slice()}, 0); //torch::Tensor(0)
    // activities[2].index_put_({100}, 1);
    // auto view_2_acc = activities[2].accessor<bool, 1>();

    const size_t ntries = 1001;
    for (size_t tries = 0; tries < ntries; ++tries) {
        auto blobs = trivial_blobs().to(device);
        // if (!tries) {
        //     std::cerr << "trivial has " << blobs.size(0) << " blobs\n";
        // }
        for (size_t view = 2; view < 5; ++view) {
            auto activity = activities[view].to(device);
            // if (!tries && view == 2) {
            //     std::cout << activity.sizes() << std::endl;
            //     // std::cout << activity << std::endl;
            // }
            blobs = apply_activity(coords, blobs, activity);
            // if(!tries)
            //     std::cout << blobs << std::endl;
            if (!tries) {
                std::cerr << "view " << view << " has " << blobs.size(0) << " blobs\n";
                // if (view < 4) {
                //     std::cerr << blobs << std::endl;
                // }
                // torch::save(blobs, "solved_blobs.pt");
                auto bytes = torch::pickle_save(blobs);
                std::ofstream fout("solved_blobs.zip", std::ios::out | std::ios::binary);
                fout.write(bytes.data(), bytes.size());
                fout.close();
            }
            assert (blobs.size(0) > 0);
        }            
        if (!tries) {
            std::cerr << "Made first set of " << blobs.size(0) << " blobs in " << sw.restart() << " us\n";
        }
    }
    double us = sw.restart();

    int ms = us/1000.0;
    std::cerr << "Repeated " << ntries-1 << " in " << ms << " ms, " << ((ntries-1) / (1e-6*us)) <<  " Hz\n";

    return 0;
}
