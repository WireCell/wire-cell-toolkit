#ifndef WIRECELL_SPNG_TDMFRAME
#define WIRECELL_SPNG_TDMFRAME

#include "WireCellSpng/TensorIndex.h"
#include "WireCellSpng/ITorchTensor.h"

namespace WireCell::SPNG {

    /// Give some structure to an SPNG TDM "frame" set of tensors.
    struct TdmFrame {

        /// Create an empty frame.  The "frame" tensor will be nullptr and the
        /// part maps will be empty.
        TdmFrame() = default;

        /// Create a frame on a TensorIndex tree node.  The node should have a
        /// value that is a "frame" tensor or be an ancestor node of such a
        /// node.
        ///
        /// Tree nodes are not held so a TdmFrame may outlive the TensorIndex.
        TdmFrame(const TensorIndex::tree_type& node);

        ~TdmFrame() = default;

        // The parent frame tensor held by frame node
        ITorchTensor::pointer frame{nullptr};

        // Maps from label (chamasks) or tag (the rest) to constituent tensors.
        std::map<std::string, ITorchTensor::pointer> traces, summaries, chids, chmasks;

    };

    /// Return the first "frame" type node with a data path that matches.
    TdmFrame frame_at_match(TensorIndex& ti, const std::string& frame_path_match);

    /// Return the index'th "frame" type node.  Default returns first frame.
    TdmFrame frame_at_index(TensorIndex& ti, size_t index=0);
}


#endif
