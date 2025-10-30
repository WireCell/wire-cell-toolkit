#ifndef WIRECELLSPNG_DNNPOSTPROCESS_H
#define WIRECELLSPNG_DNNPOSTPROCESS_H

#include "WireCellSpng/IDNNROIPostProcess.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        class DNNROIPostProcess : public IDNNROIPostProcess,
                                        public IConfigurable,
                                        public Aux::Logger {
        public:
            DNNROIPostProcess();
            virtual ~DNNROIPostProcess();
            
            void configure(const Configuration& cfg) override;
            std::vector<torch::Tensor> postprocess(
                const std::vector<torch::Tensor>& dnn_output,
                const Configuration& preprocess_metadata) override;
                
        private:
            struct Config {
                double output_scale = 1.0;
                double output_offset = 0.0;
                bool save_debug = false;
            } m_cfg;
        };
        
    }
}

#endif