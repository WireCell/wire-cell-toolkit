#include "WireCellSpng/CellBasis.h"
#include "WireCellSpng/RayGridOG.h"
#include "WireCellSpng/Ragged.h"
#include "WireCellSpng/Util.h"

using namespace torch::indexing;

namespace WireCell::SPNG::CellBasis {

    torch::Tensor wire_endpoints(const IWire::vector& wires_vec)
    {
        int64_t nwires = wires_vec.size();
        torch::Tensor wires_ten = torch::zeros({nwires, 2, 2}, torch::kDouble);
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

        // A tensor of wire endpoints of shape (NwireU, 2, 2) holding U-wire
        // endpoints in the 2D coordinate system of the Ray Grid.  Per-wire rows
        // are kept in original "pitch order".
        torch::Tensor uwires_endpoints = wire_endpoints(iplanes[0]->wires());

        // The RG layers are +2 above their usual numbers because the horiz/vert
        // layers come first.
        const auto u_layer = torch::tensor({2+0}, torch::kLong);
        const auto v_layer = torch::tensor({2+1}, torch::kLong);
        const auto w_layer = torch::tensor({2+2}, torch::kLong);

        // Find the "footprint" of each U-wire in V-view indices.
        //
        // We want the range of physically crossing V-view indices with each
        // U-wire.  On the low side, ceil gives us that.  On the high side, ceil
        // is +1 too far.  But, we will use a closed arange below so this +1 is
        // exactly what we want.
        const auto above = RayGrid::Coordinates::Rounding::kCeil;
        // (NwireU, 1)
        auto utails_in_v = coords.point_indices(uwires_endpoints.index({Slice(), 0, Slice()}),
                                                v_layer.item<int64_t>(), above);
        // (NwireU, 1)
        auto uheads_in_v = coords.point_indices(uwires_endpoints.index({Slice(), 1, Slice()}),
                                                v_layer.item<int64_t>(), above);

        // It is not trivial to know if tail or head U-wire endpoints are on the
        // low/high side as measured in V-view so we must test. But, we need
        // only test one representative.
        torch::Tensor ulo_in_v, uhi_in_v;

        if (utails_in_v[0].item<int64_t>() < uheads_in_v[0].item<int64_t>()) {
            ulo_in_v = utails_in_v;
            uhi_in_v = uheads_in_v;
        }
        else {
            ulo_in_v = uheads_in_v;
            uhi_in_v = utails_in_v;
        }        
        auto uinv_ranges = torch::stack({ulo_in_v, uhi_in_v}, 1);
        const int64_t NwiresV = iplanes[1]->wires().size();

        // Note, these values represent half-open ranges so we allow NwiresV as
        // a value.  But it is possible to get a u-in-v that is even more
        // outside the physical range.
        uinv_ranges.clamp_(0, NwiresV);

        /// The U column of the cell basis tensor
        torch::Tensor cell_u = Ragged::range_index_expansion(uinv_ranges);

        /// The V column of the cell basis tensor
        torch::Tensor cell_v = Ragged::range_value_expansion(uinv_ranges);
        
        // Find the W-view pitch of the U/V crossings
        auto w_pitch = coords.pitch_location(u_layer, cell_u, v_layer, cell_v, w_layer);

        // Find the nearest W-view index just above U/V crossing pitches
        const auto nearest = RayGrid::Coordinates::Rounding::kRound;
        torch::Tensor cell_w = coords.pitch_index(w_pitch, w_layer, nearest);

        // In principle, U/V crossings can have a nearest W-view index that is not physical.
        const int64_t NwiresW = iplanes[2]->wires().size();

        // These are indices, so must be in the physical range.  
        cell_w.clamp_(0, NwiresW-1);
        
        return torch::stack({cell_u, cell_v, cell_w}, 1);
    }


    torch::Tensor wire_channel_index(IWire::vector wires, const IChannel::vector& chans)
    {
        std::unordered_map<int, size_t> ich2ind;
        const size_t nchans = chans.size();
        for (size_t cind=0; cind<nchans; ++cind) {
            ich2ind[chans[cind]->ident()] = cind;
        }

        const size_t nwires = wires.size();
        std::vector<int64_t> out(nwires, -1);

        for (size_t wind=0; wind<nwires; ++wind) {
            auto wire = wires[wind];
            auto it = ich2ind.find(wire->channel());
            if (it != ich2ind.end()) {
                out[wind] = it->second;
            }
        }
        return to_tensor(out);
    }

    torch::Tensor channel_idents(const IChannel::vector& chans)
    {
        const size_t nchans = chans.size();
        std::vector<int> chids(nchans);
        for (size_t ind=0; ind<nchans; ++ind) {
            chids[ind] = chans[ind]->ident();
        }
        return to_tensor(chids);
    }

    torch::Tensor index(torch::Tensor basis, std::vector<torch::Tensor>& indices)
    {
        const int64_t nplanes = indices.size(); 
        std::vector<torch::Tensor> out;
        for (int64_t plane=0; plane < nplanes; ++plane) {
            auto col = basis.select(1, plane);
            out.push_back(indices[plane].index({col}));
        }
        return torch::stack(out, 1);
    }

}
