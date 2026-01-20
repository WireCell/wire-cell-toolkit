#include "WireCellSpng/CellViews.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/CellBasis.h"

#include "WireCellAux/WireTools.h"

#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/HanaJsonCPP.h"


#include <algorithm>

WIRECELL_FACTORY(SPNGCellViews,
                 WireCell::SPNG::CellViews,
                 WireCell::SPNG::ITorchTensorSetFilter,
                 WireCell::IConfigurable)

using namespace torch::indexing;
using namespace WireCell::HanaJsonCPP;
using namespace WireCell::Aux;

// For convenience in filter_tensor
using Slice = torch::indexing::Slice;

namespace WireCell::SPNG {
    CellViews::CellViews()
        : TensorSetFilter("CellViews", "spng") {}

    void CellViews::configure(const WireCell::Configuration& config)
    {
        // Propagate and parse config
        this->TensorSetFilter::configure(config);
        from_json(m_config, config);

        auto ianode = WireCell::Factory::find_tn<IAnodePlane>(m_config.anode);


        // A cell basis is defined for each anode face.  We will form a cell
        // basis tensor for each face but since the channel indices are run over
        // both faces for each view we will cat the two channel index cell basis
        // tensors to form a larger one.  After that, we no longer care about
        // faces.

        // To properly build the wire index to channel index and the cell bases
        // we have to do a pair of double loops, each in a different order.
        // The per-view ordered channels
        std::vector<IChannel::vector> channels_per_view;

        // To get ordered channels we need a per view loop over faces.  We must
        // do this first and separately because channels for a given plane
        // (view) wrap around both faces while wire (segments) do not.  We use
        // the same channels for the wire indices of either face.
        for (int view = 0; view < 3; ++view) {

            IChannel::vector all_chans;
            IWire::vector all_wires;

            for (auto face_ident : m_config.face_idents) {
                auto iface = ianode->face(face_ident);
                auto iplane = iface->planes()[view];
                
                const auto& iwires = iplane->wires();
                all_wires.insert(all_wires.end(), iwires.begin(), iwires.end());

                const auto& ichans = iplane->channels();
                all_chans.insert(all_chans.end(), ichans.begin(), ichans.end());
            }

            IChannel::vector wc_all = WireTools::wire_channels(all_wires, all_chans);
            channels_per_view.push_back(wc_all);
        }

        std::vector<torch::Tensor> cell_channel_indices_by_face;

        // Get wire indices in the cell basis for each face we need to loop over
        // views.  We then convert to channel indices on the cell basis
        for (auto face_ident : m_config.face_idents) {
            auto iface = ianode->face(face_ident);

            // Get wire indices in the cell basis.
            torch::Tensor wire_basis = CellBasis::cell_basis(iface);

            // Build tensor indexed by wire indices to return channel indices.
            std::vector<torch::Tensor> chan_indices;
            auto iplanes = iface->planes();

            for (int view = 0; view < 3; ++view) {
                IWire::vector wires = iplanes[view]->wires();
                IChannel::vector chans = channels_per_view[view];

                // Get mapping from wire index to channel index
                torch::Tensor w2c = CellBasis::wire_channel_index(wires, chans);
                chan_indices.push_back(w2c); // indexed by wire 
            }

            // Convert wire indices to channel indices in cell basis.
            torch::Tensor channel_basis = CellBasis::index(wire_basis, chan_indices);
            cell_channel_indices_by_face.push_back(channel_basis);
        }

        m_cell_channel_indices = torch::cat(cell_channel_indices_by_face, 0);

    }

    WireCell::Configuration CellViews::default_configuration() const
    {
        auto cfg = this->TensorSetFilter::default_configuration();
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }



    std::vector<torch::Tensor> CellViews::process_chunked(const std::vector<torch::Tensor>& uvw_tensors)
    {
        auto U = uvw_tensors[0];
        auto V = uvw_tensors[1];
        auto W = uvw_tensors[2];
        auto cells = m_cell_channel_indices;

        auto nbatch = U.size(0);
        auto ntick  = U.size(2);

        int64_t chunk_size = m_config.chunk_size;
        if (chunk_size <= 0) {
            chunk_size = ntick;
        }

        auto nchanU = U.size(1);
        auto nchanV = V.size(1);
        auto nchanW = W.size(1);

        // 1. Prepare indices for the three columns of cells
        auto idxU = cells.select(1, 0);
        auto idxV = cells.select(1, 1);
        auto idxW = cells.select(1, 2);

        // 2. Initialize output buffers (nbatch, 2, nchanX, ntick)
        auto outU = torch::zeros({nbatch, 2, nchanU, ntick}, torch::kBool);
        auto outV = torch::zeros({nbatch, 2, nchanV, ntick}, torch::kBool);
        auto outW = torch::zeros({nbatch, 2, nchanW, ntick}, torch::kBool);

        // 3. Loop over ticks in chunks
        for (int64_t t = 0; t < ntick; t += chunk_size) {
            int64_t actual_chunk = std::min(chunk_size, ntick - t);

            // narrow(dimension, start, length)
            auto U_slice = U.narrow(2, t, actual_chunk);
            auto V_slice = V.narrow(2, t, actual_chunk);
            auto W_slice = W.narrow(2, t, actual_chunk);

            // Indexing: (nbatch, ncell, actual_chunk)
            auto Uc = U_slice.index_select(1, idxU);
            auto Vc = V_slice.index_select(1, idxV);
            auto Wc = W_slice.index_select(1, idxW);

            // Boolean Logic
            auto MP3 = Uc & Vc & Wc;
            auto MP2u = (~Uc) & Vc & Wc;
            auto MP2v = Uc & (~Vc) & Wc;
            auto MP2w = Uc & Vc & (~Wc);

            // 4. Contract and store back to the output buffers
            for (int64_t b = 0; b < nbatch; ++b) {
                // Helper to reduce and place into the correct tick-slice
                auto scatter_chunk = [&](torch::Tensor& out_full, const torch::Tensor& data, 
                                         const torch::Tensor& indices, int64_t plane_idx) {
                    // out_full.select(0, b).select(0, plane_idx) is (nchan, ntick)
                    // we narrow the ntick dimension to the current chunk
                    auto target_slice = out_full.select(0, b).select(0, plane_idx).narrow(1, t, actual_chunk);
                
                    target_slice.index_reduce_(0, indices, data.select(0, b), "amax", false);
                };

                scatter_chunk(outU, MP3,  idxU, 0);
                scatter_chunk(outU, MP2u, idxU, 1);

                scatter_chunk(outV, MP3,  idxV, 0);
                scatter_chunk(outV, MP2v, idxV, 1);

                scatter_chunk(outW, MP3,  idxW, 0);
                scatter_chunk(outW, MP2w, idxW, 1);
            }
        }

        return {outU, outV, outW};
    }

    ITorchTensorSet::pointer CellViews::filter_tensor(const ITorchTensorSet::pointer& in)
    {
        // Extract input tensors
        auto tensors = in->tensors();
        const int ntensors_in = static_cast<int>(tensors->size());

        // Get U, V, W tensors using uvw_index configuration
        std::vector<torch::Tensor> uvw_tensors;
        for (int idx : m_config.uvw_index) {
            if (idx < 0 || ntensors_in) {
                raise<ValueError>("Invalid uvw_index %d for tensor set size %d",
                                  idx, ntensors_in);
            }
            uvw_tensors.push_back((*tensors)[idx]->tensor().to(torch::kBool));
        }

        // Check if batched: shape is either (nchan, ntick) or (nbatch, nchan, ntick)
        const bool batched = uvw_tensors[0].dim() == 3;

        // Ensure batched for uniform processing
        if (!batched) {
            for (auto& ten : uvw_tensors) {
                ten = ten.unsqueeze(0);  // Add batch dimension
            }
        }

        auto got = process_chunked(uvw_tensors);

        ITorchTensor::vector output_tensors;
        for (auto one : got) {
            if (!batched) {
                one = one.squeeze(0);
            }
            output_tensors.push_back(std::make_shared<SimpleTorchTensor>(one));
        }

        // Create output tensor set
        return std::make_shared<SimpleTorchTensorSet>(
            in->ident(),
            in->metadata(),
            output_tensors
        );
    }

}

/* Let me describe the full algorithm. I have a tensor `cells` of shape
 * `(ncells, 3)`. Each of the three columns corresponds to one of three other
 * tensors `U`, `V` and `W` of shape `(nbatch, nchanX, ntick)` where `nchanX` is
 * different for each so we have `nchanU` and `nchanV` and `nchanW`. The columns
 * of the `cells` tensor holds indices into the corresponding "chan" dimension
 * of `U`, `V` or `W` which are Boolean value. The algorithm is in two
 * phases. First, we use `cells` to index the Boolean arrays. With full
 * vectorization this would give three tensors of shape `(nbatch, ncell,
 * ntick)`. Call these `Ucell`, `Vcell` and `Wcell`. We then perform Boolean
 * operations across these three tensors. There are two types of boolean
 * operations. The first one produces a tensor called `MP3_cell` of shape
 * `(nbatch, ntick)` and which is true when the "chan" dimension of `Ucell`,
 * `Vcell` and `Wcell` are all true and false if any one of them is false. Next
 * we perform three more Boolean operations that in turn target each of the
 * `{U,V,W}cell` tensors. A `MP2u_cell` tensor is true only if the `Ucell` is
 * false and both `Vcell` and `Wcell` are true. The `MP2v_cell` and `MP2w_cell`
 * tensors are found by cyclic permuation. We then contract to form a tensor of
 * shape `(nbatch, nmp, nchanX, ntick)` for each `U`, `V` and `W`. The "mp"
 * dimension dim=1 is of size 2 and will hold in its first index values derived
 * from `MP3_cell` and in its second index values derived from the `MP2u_cell`
 * for `U`, etc for `V` and `W`. Deriving `MP3_cell` for `U` involves finding
 * all indices in the original `cells` for column 0 for which the `MP3_cell` is
 * true. We then set these chan dimension indices to true in the
 * output. Likewise we use `MP2u_cell`. Then repeat for `V` and `W`. I am
 * interested to see a fully vectorized form of this algorithm first and without
 * concern for memory usage. But if there is any "trick" to avoid full memory
 * expansion, please apply it. If no trick is available then we should examine
 * ways to mitigate memory usage possibly at the cost of introducing explicit
 * loops. */

