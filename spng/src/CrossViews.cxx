#include "WireCellSpng/CrossViews.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Util.h"
#include "WireCellSpng/RayGridOG.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IAnodePlane.h"

WIRECELL_FACTORY(SPNGCrossViews,
                 WireCell::SPNG::CrossViews,
                 WireCell::SPNG::CrossViews::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);

using namespace torch::indexing;

using WireCell::HanaJsonCPP::to_json;
using WireCell::HanaJsonCPP::from_json;

namespace WireCell::SPNG {

    CrossViews::CrossViews()
        : FaninBase<ITorchTensor>("CrossViews", "spng") { }


    void CrossViews::configure(const WireCell::Configuration& config)
    {
        // Propagate and parse config
        WireCell::configure_bases<CrossViews, FaninBase<ITorchTensor>>(this, config);
        from_json(m_cfg, config);

        // Determine which of the three input ports provides the "target"
        // tensor.
        m_targ_index = m_cfg.target_index;
        if (m_targ_index<0 || m_targ_index>2) {
            raise<ValueError>("CrossViews target index out of expected range of [0,2], inclusive");
        }

        // Sanity check the config
        size_t nfaces = m_cfg.face_idents.size();
        if (nfaces > 2) {
            raise<ValueError>("CrossViews can no understand how yhou give it %d faces", nfaces);
        }

        // Initialize the face view data. 
        m_face_view_info.clear();
        m_face_view_info.resize(nfaces);

        auto anode = Factory::find_tn<IAnodePlane>(m_cfg.anode);

        // We don't are about the sem per se but want to turn off autograd
        TorchSemaphore sem(this->context());

        log->debug("anode={} aux1={} aux2={} targ={} with {} faces",
                   anode->ident(), aux1_index(), aux2_index(), targ_index(), nfaces);

        // Transfer WCT/C++ ray grid to a Torch version for each face.
        for (size_t find=0; find<nfaces; ++find) {

            const int face_ident = m_cfg.face_idents[find];
            auto iface = anode->face(face_ident);

            auto iwireface = anode->face(face_ident);
            const auto& coords = iwireface->raygrid();

            std::vector<int64_t> layers = {0, 1,
                                           aux1_index()+2,
                                           aux2_index()+2,
                                           targ_index()+2};
            auto raygrid_views = to_spng_views(coords, layers);

            // Construct ray grid coordinates and move to device.
            auto& face_view = m_face_view_info[find];
            face_view.raygrid.init(raygrid_views);
            face_view.raygrid.to(device());
        }

        // For each input port / plane in U, V, W, find the wire<-->chan maps.
        for (size_t pind=0; pind<3; ++pind) {

            // Find the complete map from channel ID to its tensor row index
            // based on user's face_ident and WAN ordering.  This requires
            // iterating both (channel) faces.
            int nchans = 0;
            std::unordered_map<int, size_t> chid_to_index;
            for (size_t find=0; find<nfaces; ++find) {

                const int face_ident = m_cfg.face_idents[find];
                auto iface = anode->face(face_ident);

                auto iplane = iface->planes()[pind];

                for (const auto& ichan : iplane->channels()) {
                    chid_to_index[ichan->ident()] = nchans++;
                }
            }
                

            // Go through faces again to map the wires in each face to its
            // channels and vice versa
            for (size_t find=0; find<nfaces; ++find) {

                const int face_ident = m_cfg.face_idents[find];
                auto iface = anode->face(face_ident);

                // The maps for this wire face and view
                auto& view_info = m_face_view_info[find].views[pind];

                auto iplane = iface->planes()[pind];
                const auto& iwires = iplane->wires();
                const int nwires = iwires.size();

                view_info.w2c = torch::zeros({nwires},torch::kInt32);
                view_info.c2w = -1 * torch::ones({nchans},torch::kInt32);
                view_info.seg = torch::zeros({nwires},torch::kInt32);

                // Iterate along the wire array
                int n_no_channels = 0;
                for (int wind=0; wind<nwires; ++wind) {
                    auto iwire = iwires[wind];
                    const int chid = iwire->channel();
                    auto cit = chid_to_index.find(chid);
                    if (cit == chid_to_index.end()) {
                        log->warn("no channel for wire index:{}, chID:{} fID:{}, fIND:{}, pIND:{}",
                                  wind, chid, face_ident, find, pind);
                        ++n_no_channels;
                        continue;
                    }
                    const int cind = cit->second;
                    view_info.w2c[wind] = cind;
                    view_info.c2w[cind] = wind;
                    view_info.seg[wind] = iwire->segment();
                }
                if (n_no_channels) {
                    log->warn("no channel for {} out of {} wires in face index={} ident={}, plane index={}",
                              n_no_channels, nwires, find, face_ident, pind);
                }

                view_info.to(device());              
            }
        }
    } // configuration().

    
    WireCell::Configuration CrossViews::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<CrossViews, FaninBase<ITorchTensor> >(this);
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }
    

    bool CrossViews::pre_input(const input_vector& inv, std::vector<torch::Tensor>& inputs, Configuration& md)
    {
        inputs.clear();
        bool batched = false;
        for (size_t ind=0; ind<3; ++ind) {
            update(md, inv[ind]->metadata());

            auto ten = inv[ind]->tensor();
            if (ten.dim() == 2) {
                ten = ten.unsqueeze(0);
            }
            else if (ten.dim() == 3) {
                batched = true;
            }
            else {
                log->critical("unsupported tensor shape on input {}: {}", ind, to_string(ten));
                raise<ValueError>("bad tensor shape");
            }
            if (ten.dtype() != torch::kBool) {
                log->warn("coercing input {} to bool with threshold at 0.0", ind);
                ten = ten > 0;
            }

            // sanity check.  not yet sure what best to do if it fails.
            int nchans_data = ten.size(-2);
            int nchans_view = m_face_view_info[0].views[ind].nchans();
            if (nchans_data != nchans_view) {
                log->warn("channel size mismatch on input {}. input data has {}, view has {}",
                          ind, nchans_data, nchans_view);
                // for now, pray
            }

            // (nbatch, nchannel, ntick)
            log->debug("tensor on input {}: {}", ind, to_string(ten));

            inputs.push_back(ten);
        }
        return batched;
    }


    // Gemini-provided extension to NoTileMPCoincidence's version to encode
    // segment level information.
    /**
     * @brief Encode MPn values as 1-hot bit nibbles per wire-segment.
     *
     * @param input The input tensor (N_rows, N_wires) of type kInt, with values 0-3 (MP value).
     * @param indices The index tensor (N_wires) of type kLong or kInt, mapping wire index to channel index.
     * @param segment The segment tensor (N_wires) of type kInt, with values 0-2.
     * @return torch::Tensor The output tensor (N_rows, N_channels) of type kInt32, where each
     * element is a bitmask resulting from the aggregation.
     */
    static
    torch::Tensor convert_wires_to_channels_extended(
        torch::Tensor input,
        torch::Tensor indices,
        torch::Tensor segment)
    {
        // --- 1. Validate Input Shapes and Types ---
        if (input.dim() != 2 || indices.dim() != 1 || segment.dim() != 1) {
            throw std::runtime_error("Input must be (R, W), indices (W), segment (W).");
        }
        if (input.size(1) != indices.size(0) || indices.size(0) != segment.size(0)) {
            throw std::runtime_error("Wire dimension (W) must be consistent across input, indices, and segment.");
        }
        if (input.dtype() != torch::kInt32 || indices.dtype() != torch::kInt64 || segment.dtype() != torch::kInt32) {
            std::cerr << "Warning: Casting input tensors to required types (Int32/Int64)." << std::endl;
            input = input.to(torch::kInt32);
            indices = indices.to(torch::kInt64); // Indices should be Long for scatter operations
            segment = segment.to(torch::kInt32);
        }
    
        // N_rows: Number of rows/events
        const int64_t N_rows = input.size(0);
        // N_wires: Number of wires
        const int64_t N_wires = input.size(1);
    
        // Determine N_channels: max index + 1. Using max().value() is safer than max().item<int64_t>()
        // because it handles empty tensors gracefully (though indices won't be empty here).
        const int64_t N_channels = (indices.numel() > 0) 
            ? indices.max().values().item<int64_t>() + 1 
            : 0;

        if (N_channels == 0) {
            // Return an empty tensor if no channels are found
            return torch::zeros({N_rows, 0}, input.options());
        }

        // --- 2. Calculate Bit Value to Set ---
        // Goal: Determine the value to add to the aggregation for each (row, wire)
        // The MP value (0-3) determines the bit position (0-3) within the nibble.
        // The segment value (0-2) determines the nibble (0, 1, or 2) where the bit is set.
        // Bit position = MP_value + (segment * 4)
        // The final value to set is 2^(Bit position)

        // Calculate the bit shift amount for each wire: MP_value + segment * 4
        // This is a 1D tensor (N_wires)
        torch::Tensor shift_amount = segment * 4 + input.index({0, "..."}) // MP values from the first row (N_wires)
            .to(torch::kInt32); 

        // The shift_amount is currently 1D (N_wires). We need to broadcast it to (N_rows, N_wires).
        // The MP values are *different* for each row, so we need to re-calculate based on the *entire* input.
        // A simpler way: since segment is constant per wire, we can broadcast it.
    
        // Segment tensor must be broadcast to (N_rows, N_wires)
        torch::Tensor segment_broadcast = segment.unsqueeze(0).expand({N_rows, -1});
    
        // Recalculate the shift amount for the entire (N_rows, N_wires) tensor
        // MP_value is 'input' (N_rows, N_wires)
        shift_amount = segment_broadcast * 4 + input;
    
        // The value to set is 2 ^ (shift_amount)
        // Using torch::pow(2, shift_amount) is the most direct way in LibTorch.
        // The result will be a (N_rows, N_wires) tensor of type Int.
        torch::Tensor values_to_set = torch::pow(2, shift_amount).to(torch::kInt32);


        // --- 3. Prepare for Scatter Aggregation ---
    
        // The 'indices' tensor (N_wires) must be broadcast to (N_rows, N_wires)
        // to match the values_to_set shape for scatter_add_.
        torch::Tensor indices_broadcast = indices.unsqueeze(0).expand({N_rows, N_wires}).to(torch::kInt64);
    
        // Initialize the output tensor (N_rows, N_channels) with zeros (our target)
        // We use kInt32 for the resulting bitmask.
        torch::Tensor output = torch::zeros({N_rows, N_channels}, torch::kInt32);
    
        // --- 4. Perform Aggregation via torch::scatter_add_ ---
        // torch::scatter_add_(dim, index, src) performs:
        // self[index[i][j]][j] += src[i][j]   if dim=0
        // self[i][index[i][j]] += src[i][j]   if dim=1
        // We are aggregating across dim=1 (the wire/channel dimension).
        // output[row][channel_index] += values_to_set[row][wire_index]
    
        // The addition is what we want, because (2^A) + (2^B) acts as a bitwise OR
        // if A != B, which is guaranteed here since each (row, wire) sets only one bit.
        output.scatter_add_(
            1,                     // dimension 1 (the wire/channel dimension)
            indices_broadcast,     // The target channel indices
            values_to_set          // The bit value to add (2^bit_position)
            );

        return output;
    }

    torch::Tensor CrossViews::do_face(const FaceViews& face_info,
                                      std::vector<torch::Tensor>& inputs)
    {
        // input: (nbatch, nchan, ntick)

        std::vector<torch::Tensor> wire_tensors(0);

        const auto nbatch = inputs[0].size(0);
        const auto nticks = inputs[0].size(-1);
        
        // eg, nbatch=1, nticks=6000
        log->debug("do_face nbatch={} nticks={}", nbatch, nticks);

        // Convert from (nbatch,nchan,ntick) to (nbatch*ntick,nwire).  Each wire
        // tensor has same size dim=0.
        for (size_t vind=0; vind<3; ++vind) {

            auto input = inputs[vind];
            // (nbatch, nchan, ntick)
            log->debug("do_face view:{}: input:{}", vind, to_string(input));

            auto ten = input.index({"...", face_info.views[vind].w2c, Slice()});
            // (nbatch, nwire, ntick)
            log->debug("do_face view:{}: index:{}", vind, to_string(ten));

            ten = ten.transpose(-1, -2);
            // (nbatch, ntick, nwire)
            log->debug("do_face view:{}: trans:{}", vind, to_string(ten));

            ten = ten.reshape({-1, ten.size(-1)});
            // (nbatch*ntick, nwire)
            log->debug("do_face view:{}: shape:{}", vind, to_string(ten));
            wire_tensors.push_back(ten);
        }

        const auto& targ_info = face_info.views[targ_index()];

        // (nbatch*ntick, nwire)
        const auto& targ_wires = wire_tensors[targ_index()];
        const auto& aux1_wires = wire_tensors[aux1_index()];
        const auto& aux2_wires = wire_tensors[aux2_index()];

        log->debug("do_face aux1_wires={} sum={}", to_string(aux1_wires), torch::sum(aux1_wires).item<int>());
        log->debug("do_face aux2_wires={} sum={}", to_string(aux2_wires), torch::sum(aux2_wires).item<int>());
        log->debug("do_face targ_wires={} sum={}", to_string(targ_wires), torch::sum(targ_wires).item<int>());

        // Initialize the "cross views" tensor in the wires basis.  Pixel values
        // encode four mutually exclusive values: 1 is "mp1" (name coined here,
        // target is true, both cross-views false), 2 is mp2 (target false, both
        // cross views true), 3 is mp3 (all three views true) and zero otherwise
        // (all three views false).  Shape: (nbatch, nwire, ntick)
        torch::Tensor cross_views_wires = torch::zeros({
                nbatch, targ_info.nwires(), nticks}, torch::kInt32);
        log->debug("do_face cross_views_wires: {}", to_string(cross_views_wires));
        // (nbatch, nwire, ntick)

        const long int targ_nwires = targ_wires.size(1);
        // (nbatch*ntick, nwire)
        cross_views_wires = cross_views_wires.transpose(-1, -2).reshape({-1, targ_nwires});
        log->debug("do_face cross_views_wires: {}", to_string(cross_views_wires));


        // aux1 shape (nbatch*ntick, nwire_1)
        // Loop over nbatch*ntick "rows" of different sizes nwire.
        const long int nrows = targ_wires.size(0);
        for (long int irow = 0; irow < nrows; ++irow) {

            // 1D along each nwire size
            // nwires_aux1
            auto aux1_row = aux1_wires.index({irow});
            auto aux2_row = aux2_wires.index({irow});
            
            // Select indices of true aka "hit" wires.
            auto aux1_true = torch::where(aux1_row)[0];
            auto aux2_true = torch::where(aux2_row)[0];

            // Form crossing pairs of wires
            auto cross = torch::zeros({aux1_true.size(0), aux2_true.size(0), 2},
                                      tensor_options(torch::kInt32));
            cross.index_put_({"...", 1}, aux2_true);
            cross = cross.permute({1, 0, 2});
            cross.index_put_({"...", 0}, aux1_true);
            cross = cross.reshape({-1, 2});
            // Now shape: (ncrossings, 2 wire indices one in each view)
            if (irow == 0) {
                log->debug("row {} of {} cross:{}", irow, nrows, to_string(cross));
            }

            // Get the ray indices for each plane/view
            auto r1 = cross.index({Slice(), 0});
            auto r2 = cross.index({Slice(), 1});

            // Make view tensors in shapes of the ray pairs.  These are fixed
            // indices because we constructed the ray grid coordinates to match.
            auto view1 = torch::full_like(r1, 2, tensor_options(torch::kInt32));
            auto view2 = torch::full_like(r1, 3, tensor_options(torch::kInt32));
            auto view3 = torch::full_like(r1, 4, tensor_options(torch::kInt32));

            // Get the locations of crossing wires within view3
            auto view3_locs = face_info.raygrid.pitch_location(view1, r1, view2, r2, view3);

            // Also just a scalar
            auto view3_short = torch::tensor(4, tensor_options(torch::kInt32));

            // Convert the locations to pitch indices
            // 1D, (possibly many)
            if (irow == 0) {
                log->debug("row {} of {} view3_locs:{}", irow, nrows, to_string(view3_locs));
            }
            auto results = face_info.raygrid.pitch_index( view3_locs, view3_short );

            if (irow == 0) {
                log->debug("row {} of {} raygrid:{}", irow, nrows, to_string(results));
            }

            // Make sure the returned indices are within the active volume
            // 1D, small, varied
            results = std::get<0>(
                at::_unique(
                    results.index(
                        {torch::where((results >= 0) & (results < targ_nwires))[0]}
                        )));

            if (irow == 0) {
                log->debug("row {} of {} unique:{}", irow, nrows, to_string(results));
            }

            auto targ_row = targ_wires.index({irow, results});

            // |-----+--------+-------|
            // | MP# | target | cross |
            // |-----+--------+-------|
            // |   0 | low    | low   |
            // |   1 | HIGH   | low   |
            // |   2 | low    | HIGH  |
            // |   3 | HIGH   | HIGH  |
            // |-----+--------+-------|

            // MP1 and MP3 have target=HIGH.  Start by marking them all MP1
            auto targ_true = torch::nonzero(targ_row);
            if (irow == 0) {
                log->debug("row {} of {} targ_true:{}", irow, nrows, to_string(targ_true));
            }
            cross_views_wires.index_put_({irow, targ_true}, 1);

            // Overwrite some with MP3, MP3 means a high wire in the target
            // plane overlaps with a crossing pair in the other two planes
            cross_views_wires.index_put_({irow, results.index({targ_true})}, 3);
        
            // MP0 and MP2 have target=low
            // overlaps with a crossing pair in the other two planes.
            auto targ_false = torch::nonzero(targ_row == false);
            cross_views_wires.index_put_({irow, results.index({targ_false})}, 2);

        }

        auto cross_views_chans = convert_wires_to_channels_extended(
            cross_views_wires, targ_info.c2w, targ_info.seg);
        return cross_views_chans.reshape({nbatch, nticks, -1}).transpose(-1, -2);

    }

    void CrossViews::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        std::vector<torch::Tensor> inputs;
        Configuration md;
        const bool batched = pre_input(inv, inputs, md);
        
        // Initialize output tensor
        const auto& targ_chans = inputs[targ_index()];
        auto cross_views_chans = to(torch::zeros(targ_chans.sizes(), torch::kInt32));

        // Each wire face 
        const size_t nfaces = m_face_view_info.size();
        for (size_t find=0; find<nfaces; ++find) {
            const auto& face_info = m_face_view_info[find];
            auto cvc = do_face(face_info, inputs);
            cross_views_chans = torch::bitwise_or(cross_views_chans, cvc);
        }

        if (!batched) {
            cross_views_chans = cross_views_chans.squeeze(0);
        }
        out = std::make_shared<SimpleTorchTensor>(cross_views_chans, md);
    };

}
