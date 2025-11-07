#pragma once

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/ITorchTensor.h"
#include "WireCellSpng/RayGrid.h"
#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    /// Per-input view configuration.  One for each of the three input ports. 
    struct CrossViewsViewConfig {
        /// One or two face idents (not indices) that define the span of the
        /// channel dimension of this input tensor.  These are as may be given
        /// to IAnodePlane::face().  Note, for channels, the face is the one
        /// containing the segment-zero wire.  The channels will be enumerated
        /// in this "wire-attachment" aka "channel-index" order across one or
        /// two faces in the order of face_idents.  This order is expected to
        /// match the order of the channel dimension in the input tensor.  Any
        /// channels for which no wires (see "face_ident" below) are associated
        /// will be ignored, with a warning.
        std::vector<int> face_idents = {0,1};

        /// The plane ident number as may be given to IAnodeFace::plane().  This
        /// is not necessarily a plane "index" nor "layer" number though often
        /// ident and index are equated.  If not given, the index number will be
        /// set.  It will select the channels on a given face for this input
        /// tensor.
        int plane_ident=-1;

    };

    struct CrossViewsConfig {
        /// The port number supplying the "target" tensor.
        int target_index=0;           

        /// The face ident as may be given to IAnodePlane::face().  This
        /// determines the wires to consider for each plane.  Any wires on this
        /// face for which there is no channel will be discarded with a warning.
        int face_ident = -1;

        /// The per-input view config.  Three are required.
        std::vector<CrossViewsViewConfig> views;

        /// Required, name of IAnodePlane component.
        std::string anode ="";

    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::CrossViewsViewConfig, face_idents, plane_ident);
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::CrossViewsConfig, target_index, face_ident, views, anode);

namespace WireCell::SPNG {

    /** @brief A node to form a cross-views tensor from three Boolean value tensors.

        This produces the so called "mp2/mp3 multi-plane coincidence".

        The node has 3 input ports which accepts tensors from, in order, U, V
        and W planes.

        The input tensors are 2D channel vs tick images and may optionally be
        3D-batched.

        The channel dimension is expected to be in channel-index aka
        wire-attachment-number order.

     */
    struct CrossViews : public FaninBase<ITorchTensor>, virtual public IConfigurable {

        CrossViews();
        virtual ~CrossViews() = default;

        // IConfigurable - see CrossViewsConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // FaninBase API.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out);

    private:

        CrossViewsConfig m_cfg;

        /// The input port index providing the "target" tensor.  See config
        /// struct.
        int m_targ_index = 0;

        int targ_index() const { return (m_targ_index+0)%3; }
        int aux1_index() const { return (m_targ_index+1)%3; }
        int aux2_index() const { return (m_targ_index+2)%3; }

        bool pre_input(const input_vector& inv, std::vector<torch::Tensor>& inputs);

        /// The ray grid coordinates embody the wire crossing patterns.
        RayGrid::Coordinates m_raygrid;

        /// We need to translate between "row" indices between channel and wire
        /// bases separately for each view.
        struct ViewInfo {
            // wire to channel (row) indices
            torch::Tensor w2c;
            // channel to wire (row) indices
            torch::Tensor c2w;
            // the segment number of the wires
            torch::Tensor seg;            

            int nwires() const { return w2c.size(0); }
            int nchans() const { return c2w.size(0); }

            void to(torch::Device device) {
                w2c = w2c.to(device);
                c2w = c2w.to(device);
                seg = seg.to(device);
            }
        };
        std::vector<ViewInfo> m_view_info;
    };
}
