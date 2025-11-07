#ifndef WIRECELLSPNG_DNNROIPREPROCESS_H
#define WIRECELLSPNG_DNNROIPREPROCESS_H

#include "WireCellSpng/IDNNROIPreProcess.h"
#include "WireCellSpng/ITorchToTensorSet.h"
#include "WireCellSpng/ITorchForward.h"
#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        class DNNROIPreProcess : public IDNNROIPreProcess, 
                                       public IConfigurable,
                                       public Aux::Logger {
        public:
            DNNROIPreProcess();
            ~DNNROIPreProcess() override;

            Configuration default_configuration() const override;
            void configure(const Configuration& cfg) override;
            std::vector<torch::Tensor> preprocess(const ITorchTensorSet::pointer& input) override;
            Configuration get_metadata() const override;

        private:
            struct Config {
                double input_scale{1.0};
                double input_offset{0.0};
                double output_scale{1.0};
                double output_offset{0.0};
                int nchunks{1};
                int nticks{6000};
                int tick_per_slice{4};
            } m_cfg;
            
            Configuration m_metadata;
        };
        
    }
}

#endif