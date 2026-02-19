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

        // This horrendous log output is for checking that we get the right
        // channel ordering.  It's probably best to be commented out after initial debugging.
        // if (verbosity() > 2) {
        //     std::vector<IChannel::vector> cpv = CellBasis::channels_per_view(ianode, m_config.face_idents);
        //     for (const auto& ichans : cpv) {
        //         std::stringstream ss;
        //         for (const auto& ichan : ichans) {
        //             ss << " " << ichan->ident();
        //         }
        //         log->debug("{} chids: {}", ichans.size(), ss.str());
        //     }
        // }

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

        m_cell_channel_indices = to(CellBasis::cell_channel_indices(ianode, m_config.face_idents));
    }

    WireCell::Configuration CellViews::default_configuration() const
    {
        auto cfg = this->TensorSetFilter::default_configuration();
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }



    // Translate between our config and the core algorithm call.
    std::vector<torch::Tensor> CellViews::process_chunked(const std::vector<torch::Tensor>& uvw_tensors)
    {
        return CellBasis::cell_views(uvw_tensors,
                                     m_cell_channel_indices,
                                     m_config.out_views,
                                     m_config.cell_views,
                                     m_config.chunk_size);
    };


    ITorchTensorSet::pointer CellViews::filter_tensor(const ITorchTensorSet::pointer& in)
    {
        // Extract input tensors
        auto tensors = in->tensors();
        const int ntensors_in = static_cast<int>(tensors->size());

        // Get U, V, W tensors using uvw_index configuration
        std::vector<torch::Tensor> uvw_tensors;
        std::vector<WireCell::Configuration> input_configs;
        for (int idx : m_config.uvw_index) {
            if (idx < 0 || idx >= ntensors_in) {
                raise<ValueError>("Invalid uvw_index %d for tensor set size %d",
                                  idx, ntensors_in);
            }
            uvw_tensors.push_back((*tensors)[idx]->tensor().to(torch::kBool));
            input_configs.push_back((*tensors)[idx]->metadata());
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
        size_t output_count = 0;
        for (auto one : got) {
            if (!batched) {
                one = one.squeeze(0);
            }
            output_tensors.push_back(std::make_shared<SimpleTorchTensor>(one, input_configs[output_count]));
            ++output_count;
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

