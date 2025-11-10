#pragma once

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/ITorchTensor.h"
#include "WireCellSpng/RayGrid.h"
#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    /** CrossViews expects to input 3 channel-vs-tick tensors provided in order:

        (U-plane(s), V-plane(s), W-plane(s)).

        The channel dimension of each tensor must be composed of one or two
        per-face blocks.  Within each block, the rows must be in
        wire-attachment-number order.  The face IDENT of each block must be
        given as must its plane IDENT.
    */

    struct CrossViewsConfig {
        /// The port number supplying the "target" tensor.
        int target_index=0;           

        /// The face IDENT numbers.  enumerating one or two per-face blocks of
        /// channels spanned by the input tensors.  These are face IDENT number
        /// as given in a wires-file and not necessarily face index no "which"
        /// numbers.  They are as consumed by IAnodePlane::face(face_ident).
        std::vector<int> face_idents = {0,1};

        /// Required, name of IAnodePlane component corresponding to the channels.
        std::string anode ="";

    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::CrossViewsConfig, target_index, face_idents, anode);

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

        /// The input port index providing the "target" tensor.  See the config
        /// struct for details.  
        int m_targ_index = 0;

        /// The input port may not be zero, cyclical iteration finds the two
        /// "other" ports.
        int targ_index() const { return (m_targ_index+0)%3; }
        int aux1_index() const { return (m_targ_index+1)%3; }
        int aux2_index() const { return (m_targ_index+2)%3; }

        // Handle ingesting the input tensors.  Return true if input is batched.
        bool pre_input(const input_vector& inv, std::vector<torch::Tensor>& inputs, Configuration& md);

        /// A "view" here represents one layer / plane of wires on one face.
        /// Given wrapped wires, channels from both faces must be considered.
        /// Eg, for a DUNE APA, we have two sets of three "views".
        struct ViewInfo {

            // Wire (in a face) to channel (row) indices.
            torch::Tensor w2c;
            // channel to wire (row) indices (in a face).
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

        struct FaceViews {
            // Per-view info for views in the face.
            std::vector<ViewInfo> views = std::vector<ViewInfo>(3);
            
            // The ray grid wire crossings in the face.
            RayGrid::Coordinates raygrid;
        };

        // Do basic CrossViews algorithm on one wire face.
        torch::Tensor do_face(const FaceViews& face_info,
                              std::vector<torch::Tensor>& inputs);

        // A per-face view infos.  Index in this vector is same as index into
        // m_cfg.face_idents.
        using AllViews = std::vector<FaceViews>;
        AllViews m_face_view_info;
    };
}
