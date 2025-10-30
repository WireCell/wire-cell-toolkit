#ifndef WIRECELLSPNG_IDNNROIPOSTPROCESS_H
#define WIRECELLSPNG_IDNNROIPOSTPROCESS_H

#include "WireCellUtil/IComponent.h"
#include "WireCellSpng/ITorchTensorSet.h"
#include <torch/torch.h>

namespace WireCell {
    namespace SPNG {

        class IDNNROIPostProcess : public IComponent<IDNNROIPostProcess> {
        public:
            virtual ~IDNNROIPostProcess() {}

            // Process DNN output and return final result
            virtual std::vector<torch::Tensor> postprocess(
                const std::vector<torch::Tensor>& dnn_output,
                const Configuration& preprocess_metadata) = 0;
        };
        
    }
}

#endif