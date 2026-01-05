#include "WireCellSpng/CellBasis.h"
#include "WireCellSpng/RayGridOG.h"

using namespace torch::indexing;

namespace WireCell::SPNG::CellBasis {

    torch::Tensor wire_endpoints(const IWire::vector& wires_vec)
    {
        int64_t nwires = wires_vec.size();
        torch::Tensor wires_ten = torch::tensor({nwires, 2, 2}, torch::kDouble);
        for (int64_t index=0; index<nwires; ++index) {
            const auto& ray = wires_vec[index]->ray();
            wires_ten.index_put_({index,0,0},  ray.first.z());
            wires_ten.index_put_({index,0,1},  ray.first.y());
            wires_ten.index_put_({index,1,0},  ray.second.z());
            wires_ten.index_put_({index,1,1},  ray.second.y());
        }
            
        return wires_ten;
    }

    torch::Tensor cell_basis(IAnodeFace::pointer face)
    {
        auto rg_views = to_spng_views(face->raygrid());
        RayGrid::Coordinates coords(rg_views);
        auto iplanes = face->planes();

        // (Nwires, 2, 2)
        torch::Tensor uwires_endpoints = wire_endpoints(iplanes[0]->wires());

        // Fixme: u-tail is not necessarily lowest in v-view and etc u-head highest.
        auto utails_in_v = coords.point_indices(uwires_endpoints.index({Slice(), 0, Slice()}),
                                                2+1, RayGrid::Coordinates::Rounding::kCeil);
        auto uheads_in_v = coords.point_indices(uwires_endpoints.index({Slice(), 1, Slice()}),
                                                2+1, RayGrid::Coordinates::Rounding::kFloor);


        return torch::tensor({}); // dummy for now
                             
    }


}
