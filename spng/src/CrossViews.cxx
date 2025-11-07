#include "WireCellSpng/CrossViews.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IAnodePlane.h"

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

        // Sanity check config
        if (m_cfg.face_ident < 0) {
            raise<ValueError>("CrossViews requires a .face_ident configured to select wires");
        }

        if (m_cfg.views.size() != 3 ) {
            raise<ValueError>("CrossViews requires three .views configured");
        }
        for (int vind=0; vind<3; ++vind) {
            auto vcfg = m_cfg.views[vind];
            if (vcfg.plane_ident < 0) {
                vcfg.plane_ident = vind; // default conflation of ident and index
            }
            if (vcfg.face_idents.empty() || vcfg.face_idents.size() > 2) {
                raise<ValueError>("CrossViews requires one or two faces configured for each view");
            }
        }

        // Determine which of the three input ports provides the "target"
        // tensor.
        m_targ_index = m_cfg.target_index;

        // Transfer WCT/C++ ray grid to Torch version.
        auto anode = Factory::find_tn<IAnodePlane>(m_cfg.anode);
        auto iwireface = anode->face(m_cfg.face_ident);
        const auto& coords = iwireface->raygrid();
        const auto& centers = coords.centers();
        const auto& pitch_dirs = coords.pitch_dirs();
        const auto& pitch_mags = coords.pitch_mags();
        auto next_rays = centers;
        for (int ilayer = 0; ilayer < coords.nlayers(); ++ilayer) {
            next_rays[ilayer] += pitch_dirs[ilayer]*pitch_mags[ilayer];
        }
        auto raygrid_views = torch::zeros({5, 2, 2});
        // Hard-wire layer the indices 
        std::vector<int> layers = {0, 1, aux1_index()+2, aux2_index()+2, targ_index()+2};
        for (int layer_count=0; layer_count<5; ++layer_count) {
            const int ilayer = layers[layer_count];
            raygrid_views.index_put_({layer_count, 0, 0}, centers[ilayer][2]);
            raygrid_views.index_put_({layer_count, 0, 1}, centers[ilayer][1]);
            raygrid_views.index_put_({layer_count, 1, 0}, next_rays[ilayer][2]);
            raygrid_views.index_put_({layer_count, 1, 1}, next_rays[ilayer][1]);
        }
        // Construct our ray grid coordinates and move to our device.
        m_raygrid.init(raygrid_views);
        m_raygrid.to(device());


        // Populate the plane info.
        m_view_info.clear();
        m_view_info.resize(3);
        for (int vind=0; vind<3; ++vind) {
            auto pi = m_view_info[vind];

            const auto& vcfg = m_cfg.views[vind];
            const int plane_ident = vcfg.plane_ident;

            // Initialize wires-to-chans
            auto iwireplane = iwireface->plane(plane_ident);
            const auto& iwires = iwireplane->wires();
            const int nwires = iwires.size();
            pi.w2c = torch::zeros({nwires},torch::kInt32);

            // Initialize intermediate channel ID to channel index following
            // users order of channel faces.
            std::unordered_map<int, size_t> chid_to_index;
            int nchans = 0;
            for (int chan_face_ident : vcfg.face_idents) {
                auto chan_face = anode->face(chan_face_ident);
                for (const auto& ichan : chan_face->plane(plane_ident)->channels()) {
                    chid_to_index[ichan->ident()] = nchans++;
                }
            }
            pi.c2w = -1 * torch::ones({nchans},torch::kInt32);

            // Iterate along the wire array
            for (int wind=0; wind<nwires; ++wind) {
                auto iwire = iwires[wind];
                const int chid = iwire->channel();
                auto cit = chid_to_index.find(chid);
                if (cit == chid_to_index.end()) {
                    log->warn("wire:{} has no channel for chid:{} face ID:{} plane ID:{} not given",
                              wind, chid, m_cfg.face_ident, plane_ident);
                    continue;
                }
                const int cind = cit->second;
                pi.w2c[wind] = cind;
                pi.c2w[cind] = wind;
                pi.seg[wind] = iwire->segment();
            }

            // Warn user if extra channels are given.
            for (int cind=0; cind<nchans; ++cind) {
                if (pi.c2w[cind].item<int>() >= 0) {
                    continue;
                }
                log->warn("ch index:{} has no wire in face ID:{} and plane ID:{} not given",
                          cind, m_cfg.face_ident, plane_ident);
            }

            pi.to(device());                
        }

    }
    
    WireCell::Configuration CrossViews::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<CrossViews, FaninBase<ITorchTensor> >(this);
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }
    


    bool CrossViews::pre_input(const input_vector& inv, std::vector<torch::Tensor>& inputs)
    {
        inputs.clear();
        bool batched = false;
        for (int ind=0; ind<3; ++ind) {
            auto ten = inv[ind]->tensor();
            if (ten.dim() == 2) {
                ten = ten.unsqueeze(0);
            }
            else {
                batched = true;
            }
            if (ten.dtype() != torch::kBool) {
                log->warn("coercing input {} to bool with threshold at 0.0", ind);
                ten = ten > 0;
            }

            // sanity check.  not yet sure what best to do if it fails.
            int nchans = ten.size(-2);
            if (nchans != m_view_info[ind].nchans()) {
                log->warn("channel size mismatch on input {}. input data has {}, view has {}",
                          ind, nchans, m_view_info[ind].nchans());
                // for now, pray
            }

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


    void CrossViews::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        std::vector<torch::Tensor> inputs;
        const bool batched = pre_input(inv, inputs);
        
        // Convert from (nbatch,nchan,ntick) to (nbatch*ntick,nwire).  Each wire
        // tensor has same size dim=0.
        Configuration md;
        std::vector<torch::Tensor> wire_tensors;
        for (int vind=0; vind<3; ++vind) {
            /// FIXME: probably not the best way to set the output MD....
            update(md, inv[vind]->metadata());
            

            auto input = inputs[vind];
            auto ten = input.index({Ellipsis, m_view_info[vind].w2c, Ellipsis});
            ten = ten.transpose(-1, -2);
            ten = ten.reshape({-1, ten.size(-1)});
            wire_tensors.push_back(ten);
        }

        // Give semantic aliases.
        const auto& targ_info = m_view_info[targ_index()];
        const auto& targ_wires = wire_tensors[targ_index()];
        const int targ_nwires = targ_wires.size(-2);

        auto nbatch = targ_wires.size(0);
        auto nticks = targ_wires.size(-1);
        
        // Initialize the "cross views" tensor in the wires basis.  Pixel values
        // encode four mutually exclusive values: 1 is "mp1" (name coined here,
        // target is true, both cross-views false), 2 is mp2 (target false, both
        // cross views true), 3 is mp3 (all three views true) and zero otherwise
        // (all three views false).  Shape: (nbatch, nwire, ntick)
        torch::Tensor cross_views_wires = torch::zeros({
                nbatch, targ_info.nwires(), nticks}, torch::kInt);

        

        const auto& aux1_wires = wire_tensors[aux1_index()];
        const auto& aux2_wires = wire_tensors[aux2_index()];
        const int aux1_nwires = aux1_wires.size(-2);

        // Loop over nbatch*ntick "rows" of different sizes nwire.
        for (long int irow = 0; irow < aux1_nwires; ++irow) {

            // 1D along each nwire size
            auto aux1_row = aux1_wires.index({irow});
            auto aux2_row = aux2_wires.index({irow});
            
            // Select indices of true aka "hit" wires.
            auto aux1_true = torch::where(aux1_row)[0];
            auto aux2_true = torch::where(aux2_row)[0];

            // Form crossing pairs of wires
            auto cross = torch::zeros({aux1_true.size(0), aux2_true.size(0), 2},
                                      tensor_options(torch::kInt));
            cross.index_put_({Ellipsis, 1}, aux2_true);
            cross = cross.permute({1, 0, 2});
            cross.index_put_({Ellipsis, 0}, aux1_true);
            cross = cross.reshape({-1, 2});
            // Now shape: (ncrossings, 2 wire indices one in each view)

            // Get the ray indices for each plane/view
            auto r1 = cross.index({Slice(), 0});
            auto r2 = cross.index({Slice(), 1});

            // Make view tensors in shapes of the ray pairs.  These are fixed
            // indices because we constructed the ray grid coordinates to match.
            auto view1 = torch::full_like(r1, 2, tensor_options(torch::kInt));
            auto view2 = torch::full_like(r1, 3, tensor_options(torch::kInt));
            auto view3 = torch::full_like(r1, 4, tensor_options(torch::kInt));

            // Get the locations of crossing wires within view3
            auto view3_locs = m_raygrid.pitch_location(view1, r1, view2, r2, view3);

            // Also just a scalar
            auto view3_short = torch::tensor(4, tensor_options(torch::kInt));

            // Convert the locations to pitch indices
            // 1D, small, varied
            auto results = m_raygrid.pitch_index( view3_locs, view3_short );

            // Make sure the returned indices are within the active volume
            // 1D, small, varied
            results = std::get<0>(
                at::_unique(
                    results.index(
                        {torch::where((results >= 0) & (results < targ_nwires))[0]}
                        )));

            auto targ_row = targ_wires.index({irow, results});

            // Start by setting all target high wires as MP1.
            auto targ_true = torch::nonzero(targ_row);
            cross_views_wires.index_put_({irow, targ_true}, 1);

            // Overwrite some with MP3, MP3 means a high wire in the target
            // plane overlaps with a crossing pair in the other two planes
            cross_views_wires.index_put_({irow, results.index({targ_true})}, 3);
        
            // MP2 means a low wire in the target plane
            // overlaps with a crossing pair in the other two planes.
            auto targ_false = torch::nonzero(targ_row == false);
            cross_views_wires.index_put_({irow, results.index({targ_false})}, 2);

        }

        auto cross_views_chans = convert_wires_to_channels_extended(
            cross_views_wires, targ_info.w2c, targ_info.seg);
        cross_views_chans = cross_views_chans.reshape({nbatch, nticks, -1}).transpose(-1, -2);

        if (!batched) {
            cross_views_chans = cross_views_chans.squeeze(0);
        }
        out = std::make_shared<SimpleTorchTensor>(cross_views_chans, md);
    };

}
