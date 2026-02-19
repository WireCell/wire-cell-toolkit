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

        // Find the nearest W-view index near U/V crossing pitches
        // const auto nearest = RayGrid::Coordinates::Rounding::kCeil;
        const auto nearest = RayGrid::Coordinates::Rounding::kRound;
        torch::Tensor cell_w = coords.pitch_index(w_pitch, w_layer, nearest);

        // In principle, U/V crossings can have a nearest W-view index that is not physical.
        const int64_t NwiresW = iplanes[2]->wires().size();


        // These are indices, so must be in the physical range.  Indices should
        // be out of range at most by +/- 1 because the U and V crossings that
        // are considered are inside the boundaries defined by their wire
        // endpoints.  If you can tolerate the noise, test that only +/-1
        // clamping happens by uncommenting:
        // if (torch::any(cell_w < 0).item<bool>()) {
        //     std::cerr << "CELL_BASIS most underflowed W wire: "
        //               << torch::min(cell_w).item<int64_t>() <<  "\n";
        // }
        // if (torch::any(cell_w >= NwiresW).item<bool>()) {
        //     std::cerr << "CELL_BASIS most overflowflowed W wire: "
        //               << torch::max(cell_w).item<int64_t>() <<  "\n";
        // }
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
        auto ret = to_tensor(out);
        if (torch::any(ret < 0).item<bool>()) {
            std::cerr << "WIRE_CHANNEL_INDEX given a wire in "<<wires.size()<<" wires with no channel in "<<chans.size()<<" channels\n";
        }
        return ret;
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

    
    std::vector<IChannel::vector> channels_per_view(IAnodePlane::pointer ianode,
                                                    const std::vector<int>& face_idents)
    {
        // Channels ordered along both faces of each view.
        std::vector<IChannel::vector> cpv;
        for (int view = 0; view < 3; ++view) {

            IChannel::vector all_chans;

            for (auto face_ident : face_idents) {
                auto iface = ianode->face(face_ident);
                auto iplane = iface->planes()[view];
                
                const auto& ichans = iplane->channels();
                all_chans.insert(all_chans.end(), ichans.begin(), ichans.end());

                // dump_chans("face chan idents", ichans, face_ident, view);
                // log->debug("view={} face={} nchannels={}",
                //            view, face_ident, ichans.size());
            }
            // All chans: (400+400, 400+400, 480+480) = (800, 800, 960)
            // All wires: (1148+1148, 1148+1148, 480+480) = (2296,2296,960)
            // wc_all is (2296,2296,960), for each wire, its channel INDEX.
            // dump_chans("all chan idents before", all_chans, -1, view);
            cpv.push_back(all_chans);
        }
        return cpv;
    }

    torch::Tensor cell_channel_indices(IAnodePlane::pointer ianode,
                                       const std::vector<int>& face_idents)
    {
        // Channels ordered along both faces of each view.
        std::vector<IChannel::vector> cpv = channels_per_view(ianode, face_idents);

        std::vector<torch::Tensor> cell_channel_indices_by_face;

        // Get wire indices in the cell basis for each face.
        for (auto face_ident : face_idents) {
            auto iface = ianode->face(face_ident);

            // Get wire indices in the cell basis.
            torch::Tensor wire_basis = cell_basis(iface);
            // dump_basis("wire indices", wire_basis, face_ident);

            // Build tensor mapping wire index to channel index.
            std::vector<torch::Tensor> w2c_per_view;

            auto iplanes = iface->planes();
            for (int view = 0; view < 3; ++view) {
                IWire::vector wires = iplanes[view]->wires();
                // dump_wires("w2c wires", wires, face_ident, view);
                const IChannel::vector& chans = cpv[view];
                // dump_chans("all chan idents after", chans, face_ident, view);

                // Get mapping from wire index to channel index
                torch::Tensor w2c = wire_channel_index(wires, chans);
                // dump_basis("w2c indices", w2c, face_ident, view);
                w2c_per_view.push_back(w2c); // indexed by wire 
            }

            // Convert wire indices to channel indices in cell basis.
            torch::Tensor chan_basis = index(wire_basis, w2c_per_view);
            // dump_basis("chan indices", chan_basis, face_ident);
            cell_channel_indices_by_face.push_back(chan_basis);
        }

        return torch::cat(cell_channel_indices_by_face, 0);
    }

    // Strategy of this function is to heavily vectorize construction of MP2 and
    // MP3.  However, full vectorization over the sizes of some detectors can
    // exhaust memory.  To combat that, the vectorization can be chunked over
    // the tick dimension.  This function has a triple nested loop: chunks,
    // batch and output view.  The outer loops have low cardinality and the
    // inner kernel of the full nested loops is "heavily vectorized".  Thus, we
    // expect this function to still be very performant.  We further reduce its
    // cost by not calculating MP2 and final scatter for views that are not
    // requested.  There is a 3-view overhead to make MP3 even if just one
    // output view is requested so it is best to make all output views together.
    std::vector<torch::Tensor> cell_views(const std::vector<torch::Tensor>& uvw_tensors,
                                          const torch::Tensor& cells,
                                          const std::vector<int>& out_views,
                                          const std::vector<std::string>& cell_views,
                                          int64_t chunk_size)
    {
        auto device = cells.device();

        auto U = uvw_tensors[0].to(device);
        auto V = uvw_tensors[1].to(device);
        auto W = uvw_tensors[2].to(device);

        auto nbatch = U.size(0);
        auto ntick  = U.size(2);

        // Fixme: this could be done in configure()
        // This controls the order of feature dimension.
        std::map<std::string, size_t> cell_view_index;
        for (size_t ind=0; ind<cell_views.size(); ++ind) {
            cell_view_index[cell_views[ind]] = ind;
        }

        if (chunk_size <= 0) {
            chunk_size = ntick;
        }

        // 1. Prepare indices for the three columns of cells
        std::vector<torch::Tensor> idx;
        for (int view=0; view<3; ++view) {
            idx.push_back(cells.select(1, view));
        }

        // 2. Initialize output buffers (nbatch, 2, nchanX, ntick)
        std::vector<torch::Tensor> out;
        for (int view : out_views) {
            int64_t nchan = uvw_tensors[view].size(1);
            out.push_back(torch::zeros({nbatch, 2, nchan, ntick}, torch::kBool).to(device));
        }

        // 3. Loop over chunks of ticks.
        for (int64_t t = 0; t < ntick; t += chunk_size) {
            int64_t actual_chunk = std::min(chunk_size, ntick - t);

            // Regardless of the number of out's we need some processing of all three views

            // narrow(dimension, start, length)
            auto U_slice = U.narrow(2, t, actual_chunk);
            auto V_slice = V.narrow(2, t, actual_chunk);
            auto W_slice = W.narrow(2, t, actual_chunk);

            // Indexing: (nbatch, ncell, actual_chunk)
            auto Uc = U_slice.index_select(1, idx[0]);
            auto Vc = V_slice.index_select(1, idx[1]);
            auto Wc = W_slice.index_select(1, idx[2]);

            // Need MP3 for all output views.
            auto MP3 = Uc & Vc & Wc;

            // Only need MP2 for output views
            std::vector<torch::Tensor> MP2;
            for (int view : out_views) {
                switch (view) {
                case 0: MP2.push_back( (~Uc) & Vc & Wc ); break;
                case 1: MP2.push_back( Uc & (~Vc) & Wc ); break;
                case 2: MP2.push_back( Uc & Vc & (~Wc) ); break;
                }
            }

            // 4. Contract and store back to the output buffers
            for (int64_t b = 0; b < nbatch; ++b) {

                // Helper to reduce and place into the correct tick-slice
                auto scatter_chunk = [&](torch::Tensor& out_full, const torch::Tensor& data,
                                         const torch::Tensor& indices, int64_t plane_idx) {
                    // out_full.select(0, b).select(0, plane_idx) is (nchan, ntick)
                    // we narrow the ntick dimension to the current chunk
                    auto target_slice = out_full.select(0, b).select(0, plane_idx).narrow(1, t, actual_chunk);

                    // Workaround: index_reduce_ doesn't support bool, so cast to uint8, reduce, cast back
                    auto target_u8 = target_slice.to(torch::kUInt8);
                    target_u8.index_reduce_(0, indices, data.select(0, b).to(torch::kUInt8), "amax", false);
                    target_slice.copy_(target_u8.to(torch::kBool));
                };

                // Only compute requested output views.
                for (size_t view_index=0; view_index<out_views.size(); ++view_index) {
                    int view = out_views[view_index];
                    scatter_chunk(out[view_index], MP2[view_index], idx[view], cell_view_index["mp2"]);
                    scatter_chunk(out[view_index], MP3, idx[view], cell_view_index["mp3"]);
                }                    
            }
        }

        return out;
    }

}
