#pragma once

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/Torch.h"


namespace WireCell::SPNG {

    struct CrossViewsExtractConfig {
        /** An extraction names an operation to apply to the input cross views tensor.
         *
         * See CrossViews for expectations of the content of the input tensor.
         *
         * Each extraction produces a tensor to an output port.  That is, the
         * length of the extraction list determines the multiplicity of the
         * fanout.
         *
         * The supported extractions:
         *
         * - mp2 :: If a pixel has any MP2 bit, the output pixel is true, else false.
         * - mp3 :: If a pixel has any MP3 bit, the output pixel is true, else false.
         *
         * Default is ["mp2","mp3"].
         */
        std::vector<std::string> extraction = {"mp2","mp3"}; 
    };

}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::CrossViewsExtractConfig, extraction);

namespace WireCell::SPNG {

    /** Extract individual images from a CrossViews tensor.

        See CrossViewsExtractConfig for more documentation.
     */
    struct CrossViewsExtract : FanoutBase<ITorchTensor>, public virtual IConfigurable {
        CrossViewsExtract();
        virtual ~CrossViewsExtract();

        // FanoutBase API.
        virtual void fanout_separate(const input_pointer& in, output_vector& outv);

        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;

    private:

        CrossViewsExtractConfig m_config;

        using Operation = std::function<torch::Tensor(const torch::Tensor&)>;
        std::vector<Operation> m_ops;

    };

}
