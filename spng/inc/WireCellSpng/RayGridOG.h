/// Convert from WCT's original RayGrid to SPNG's Torch-based RayGrid

#pragma once

#include "WireCellUtil/RayGrid.h"
#include "WireCellSpng/RayGrid.h"

namespace WireCell::SPNG {

    /// Return a (N-views, 2 endpoints, 2 coordinates) tensor
    /// Suitable for giving to SPNG::RayGrid::Coordinates().
    torch::Tensor to_spng_views(const WireCell::RayGrid::Coordinates& og);

    /// As above but control the order of the views to follow that in the given
    /// vector of indices into the OG layers.  Eg, the original layer index
    /// given by og_layer_indices[0] will be in the new layer index 0.
    torch::Tensor to_spng_views(const WireCell::RayGrid::Coordinates& og,
                                std::vector<int64_t> og_layer_indices);

    /// Convert 3D {x,y,z} point to SPNG RG 2D {z,y} point.
    torch::Tensor to_spng_2d(const Point& pt);


}

