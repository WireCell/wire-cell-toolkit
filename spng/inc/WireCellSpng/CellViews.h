#pragma once

#include "WireCellSpng/TensorSetFilter.h"
#include <torch/torch.h>

namespace WireCell::SPNG {

    struct CellViewsConfig {
        /// The face IDENT numbers.  enumerating one or two per-face blocks of
        /// channels spanned by the input tensors.  These are face IDENT number
        /// as given in a wires-file and not necessarily face index no "which"
        /// numbers.  They are as consumed by IAnodePlane::face(face_ident).
        ///
        /// The resulting channel order MUST match the input per-plane row order.
        std::vector<int> face_idents = {0,1};

        /// The indices for the U, V and W tensors in the input tensor set.
        std::vector<int> uvw_index = {0,1,2};

        /// Required, name of IAnodePlane component corresponding to the channels.
        std::string anode ="";

        /// Name and order of the kinds of "cell views" to produce.  Each cell
        /// view maps to one index on the dim=-3 dimension of the output
        /// tensors.  
        std::vector<std::string> cell_views = {"mp2", "mp3"};

        /// If set, perform chunked processing
        int chunk_size = 0;
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::CellViewsConfig, face_idents, uvw_index, anode, cell_views, chunk_size);

namespace WireCell::SPNG {

    /**
       CellViews is similar to CrossViews+CrossViewsExtract with these differences:

       - It uses a faster algorithm.
       - It combines finding fundamental MP info with extracting mp2 and mp3 images.
       - It consumes and produces a tensor set.

       The input tensor set is expected to have at least 3 tensors selected via
       uvw_index configuration to provide U, V and W plane views.

       Input tensors are shaped either:

       (nchan, ntick)

       or

       (nbatch, nchan, ntick)

       Output tensor set contains exactly 3 tensors in order U, V and W.

       Output tensors are shaped either:

       (ncellview, nchan, ntick)

       or

       (nbatch, ncellview, nchan, ntick)

       The ncellview counts the number of "cell view" images to be produced.
       The cell views are named in the cell_views configuration variable.

     */
    struct CellViews : public TensorSetFilter {
        CellViews();
        virtual ~CellViews() = default;
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        virtual ITorchTensorSet::pointer filter_tensor(const ITorchTensorSet::pointer& in);
        
    private:
        
        std::vector<torch::Tensor> process_chunked(const std::vector<torch::Tensor>& uvw_tensors);

        CellViewsConfig m_config;
        // (3, ncell) holding channel dimension indices, one column for each view.
        torch::Tensor m_cell_channel_indices;

    };
}
