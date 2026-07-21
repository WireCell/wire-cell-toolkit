#ifndef WIRECELLSPNG_DNNPOSTPROCESS_H
#define WIRECELLSPNG_DNNPOSTPROCESS_H

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        class DNNROIPostProcess : public ITorchTensorSetFilter,
                                        public IConfigurable,
                                        public Aux::Logger {
        public:
            DNNROIPostProcess();
            virtual ~DNNROIPostProcess();

            virtual void configure(const Configuration& cfg);
            virtual Configuration default_configuration() const;
            virtual bool operator()(const input_pointer &in, output_pointer& out);
        private:
            struct Config {
                double output_scale{1.0};
                double output_offset{0.0};
                int ntick{6000};
                int nchunk{1};
                int tick_per_slice{4};
            } m_cfg;
        };
        
    }
}

#endif