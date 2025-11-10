#ifndef WIRECELLSPNG_DNNROIPREPROCESS_H
#define WIRECELLSPNG_DNNROIPREPROCESS_H

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        class DNNROIPreProcess : public ITorchTensorSetFilter, 
                                       public IConfigurable,
                                       public Aux::Logger {
        public:
            DNNROIPreProcess();
            virtual ~DNNROIPreProcess();

            virtual Configuration default_configuration() const;
            virtual void configure(const Configuration& cfg);
            virtual bool operator()(const input_pointer &in, output_pointer& out);

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