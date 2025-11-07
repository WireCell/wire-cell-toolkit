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
            
            virtual void configure(const Configuration& cfg) override;
            virtual std::vector<torch::Tensor> postprocess(
                const std::vector<torch::Tensor>& dnn_output,
                const Configuration& preprocess_metadata) override;
                
        private:
            struct Config {
                double output_scale = 1.0;
                double output_offset = 0.0;
                int ntick = 6000;
                int nchunk = 1;
                int tick_per_slice = 4;
            } m_cfg;
        };
        
    }
}

#endif