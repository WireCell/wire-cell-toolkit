#include "WireCellSpng/CellViews.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/CellBasis.h"
#include "WireCellSpng/Util.h"

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

        // auto dump_basis = [this](const std::string& name, const torch::Tensor& ten, int face, int view=-1) {
        //     if (view < 0) {
        //         auto imin = torch::amin(ten, {0});
        //         auto imax = torch::amax(ten, {0});

        //         this->log->debug("{}: face={} tensor={} min=[{} {} {}] max=[{} {} {}]",
        //                          name, face, to_string(ten),
        //                          imin[0].item<int>(), imin[1].item<int>(), imin[2].item<int>(),
        //                          imax[0].item<int>(), imax[1].item<int>(), imax[2].item<int>());
        //     }
        //     else {
        //         auto imin = torch::amin(ten);
        //         auto imax = torch::amax(ten);

        //         this->log->debug("{}: face={} view={} tensor={} min={} max={}",
        //                          name, face, view, to_string(ten),
        //                          imin.item<int>(), 
        //                          imax.item<int>());
        //     }
        // };


        // auto dump_chans = [&,this](const std::string& name, const IChannel::vector& chans, int face, int view) {
        //     std::vector<int> chids;
        //     for (const auto& ich : chans) {
        //         chids.push_back(ich->ident());
        //     }
        //     dump_basis(name, to_tensor(chids), face, view);
        //     this->log->debug("first ident={} last ident={}",
        //                      chans.front()->ident(), chans.back()->ident());
        // };

        // auto dump_wires = [&,this](const std::string& name, const IWire::vector& wires, int face, int view) {
        //     this->log->debug("{}: face={} view={} chids: first={} last={}",
        //                      name, face, view,
        //                      wires.front()->channel(),
        //                      wires.back()->channel());
        //     for (const auto& iwire : wires) {
        //         this->log->debug("wire index={} ident={} chid={} seg={} wpid={}",
        //                          iwire->index(), iwire->ident(), iwire->channel(), iwire->segment(), iwire->planeid());
        //     }
        // };

        // Channels ordered along both faces of each view.
        std::vector<IChannel::vector> channels_per_view;
        for (int view = 0; view < 3; ++view) {

            IChannel::vector all_chans;

            for (auto face_ident : m_config.face_idents) {
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
            channels_per_view.push_back(all_chans);
        }

        std::vector<torch::Tensor> cell_channel_indices_by_face;

        // Get wire indices in the cell basis for each face.
        for (auto face_ident : m_config.face_idents) {
            auto iface = ianode->face(face_ident);

            // Get wire indices in the cell basis.
            torch::Tensor wire_basis = CellBasis::cell_basis(iface);
            // dump_basis("wire indices", wire_basis, face_ident);

            // Build tensor mapping wire index to channel index.
            std::vector<torch::Tensor> w2c_per_view;

            auto iplanes = iface->planes();
            for (int view = 0; view < 3; ++view) {
                IWire::vector wires = iplanes[view]->wires();
                // dump_wires("w2c wires", wires, face_ident, view);
                const IChannel::vector& chans = channels_per_view[view];
                // dump_chans("all chan idents after", chans, face_ident, view);

                // Get mapping from wire index to channel index
                torch::Tensor w2c = CellBasis::wire_channel_index(wires, chans);
                // dump_basis("w2c indices", w2c, face_ident, view);
                w2c_per_view.push_back(w2c); // indexed by wire 
            }

            // Convert wire indices to channel indices in cell basis.
            torch::Tensor chan_basis = CellBasis::index(wire_basis, w2c_per_view);
            // dump_basis("chan indices", chan_basis, face_ident);
            cell_channel_indices_by_face.push_back(chan_basis);
        }

        m_cell_channel_indices = to(torch::cat(cell_channel_indices_by_face, 0));
        // log->debug("cell channel indices: {}", to_string(m_cell_channel_indices));
    }

    WireCell::Configuration CellViews::default_configuration() const
    {
        auto cfg = this->TensorSetFilter::default_configuration();
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
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
    std::vector<torch::Tensor> CellViews::process_chunked(const std::vector<torch::Tensor>& uvw_tensors)
    {
        auto U = uvw_tensors[0];
        auto V = uvw_tensors[1];
        auto W = uvw_tensors[2];
        auto cells = m_cell_channel_indices;

        auto nbatch = U.size(0);
        auto ntick  = U.size(2);

        // Fixme: this could be done in configure()
        // This controls the order of feature dimension.
        std::map<std::string, size_t> cell_view_index;
        for (size_t ind=0; ind<m_config.cell_views.size(); ++ind) {
            cell_view_index[m_config.cell_views[ind]] = ind;
        }

        int64_t chunk_size = m_config.chunk_size;
        if (chunk_size <= 0) {
            chunk_size = ntick;
        }

        // auto nchanU = U.size(1);
        // auto nchanV = V.size(1);
        // auto nchanW = W.size(1);

        // 1. Prepare indices for the three columns of cells
        std::vector<torch::Tensor> idx;
        for (int view=0; view<3; ++view) {
            idx.push_back(cells.select(1, view));
        }
        // auto idxU = cells.select(1, 0);
        // auto idxV = cells.select(1, 1);
        // auto idxW = cells.select(1, 2);

        // 2. Initialize output buffers (nbatch, 2, nchanX, ntick)
        std::vector<torch::Tensor> out;
        for (int view : m_config.out_views) {
            int64_t nchan = uvw_tensors[view].size(1);
            out.push_back(to(torch::zeros({nbatch, 2, nchan, ntick}, torch::kBool)));
        }
        // auto outU = torch::zeros({nbatch, 2, nchanU, ntick}, torch::kBool);
        // auto outV = torch::zeros({nbatch, 2, nchanV, ntick}, torch::kBool);
        // auto outW = torch::zeros({nbatch, 2, nchanW, ntick}, torch::kBool);

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
            for (int view : m_config.out_views) {
                switch (view) {
                case 0: MP2.push_back( (~Uc) & Vc & Wc ); break;
                case 1: MP2.push_back( Uc & (~Vc) & Wc ); break;
                case 2: MP2.push_back( Uc & Vc & (~Wc) ); break;
                }
            }
            // auto MP2u = (~Uc) & Vc & Wc;
            // auto MP2v = Uc & (~Vc) & Wc;
            // auto MP2w = Uc & Vc & (~Wc);

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
                for (size_t view_index=0; view_index<m_config.out_views.size(); ++view_index) {
                    int view = m_config.out_views[view_index];
                    scatter_chunk(out[view_index], MP2[view_index], idx[view], cell_view_index["mp2"]);
                    scatter_chunk(out[view_index], MP3, idx[view], cell_view_index["mp3"]);
                }                    
                // scatter_chunk(outU, MP3,  idxU, 0);
                // scatter_chunk(outU, MP2u, idxU, 1);

                // scatter_chunk(outV, MP3,  idxV, 0);
                // scatter_chunk(outV, MP2v, idxV, 1);

                // scatter_chunk(outW, MP3,  idxW, 0);
                // scatter_chunk(outW, MP2w, idxW, 1);
            }
        }

        return out;
    }

    ITorchTensorSet::pointer CellViews::filter_tensor(const ITorchTensorSet::pointer& in)
    {
        // Extract input tensors
        auto tensors = in->tensors();
        const int ntensors_in = static_cast<int>(tensors->size());

        // Get U, V, W tensors using uvw_index configuration
        std::vector<torch::Tensor> uvw_tensors;
        for (int idx : m_config.uvw_index) {
            if (idx < 0 || idx >= ntensors_in) {
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

