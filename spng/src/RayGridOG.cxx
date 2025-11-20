#include "WireCellSpng/RayGridOG.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::SPNG {

    torch::Tensor to_spng_views(const WireCell::RayGrid::Coordinates& og,
                                std::vector<int64_t> og_layer_indices)
    {
        const int64_t nlayers = og.nlayers();
        if (nlayers != (int64_t)og_layer_indices.size()) {
            raise<ValueError>("size mismatch between RayGrid and layer indices");
        }

        const auto& centers = og.centers();
        const auto& pitch_dirs = og.pitch_dirs();
        const auto& pitch_mags = og.pitch_mags();
        auto next_rays = centers;
        for (int64_t ilayer = 0; ilayer < nlayers; ++ilayer) {
            next_rays[ilayer] += pitch_dirs[ilayer]*pitch_mags[ilayer];
        }

        auto raygrid_views = torch::zeros({nlayers, 2, 2}, torch::kDouble);

        for (int layer_count=0; layer_count<nlayers; ++layer_count) {
            const int og_layer = og_layer_indices[layer_count];
            raygrid_views.index_put_({layer_count, 0, 0}, centers[og_layer][2]);
            raygrid_views.index_put_({layer_count, 0, 1}, centers[og_layer][1]);
            raygrid_views.index_put_({layer_count, 1, 0}, next_rays[og_layer][2]);
            raygrid_views.index_put_({layer_count, 1, 1}, next_rays[og_layer][1]);
        }
        return raygrid_views;
    }

    torch::Tensor to_spng_views(const WireCell::RayGrid::Coordinates& og)
    {
        std::vector<int64_t> og_layer_indices(og.nlayers());
        std::iota(og_layer_indices.begin(), og_layer_indices.end(), 0);
        return to_spng_views(og, og_layer_indices);
    }
    
}
