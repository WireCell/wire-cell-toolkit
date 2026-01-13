#ifndef WIRECELLSPNG_IDNNROIPREPROCESS_H
#define WIRECELLSPNG_IDNNROIPREPROCESS_H

#include "WireCellUtil/IComponent.h"
#include "WireCellSpng/ITorchTensorSet.h"
#include <torch/torch.h>

namespace WireCell {
    namespace SPNG {

        class IDNNROIPreProcess : public IComponent<IDNNROIPreProcess> {
        public:
            virtual ~IDNNROIPreProcess() {}

            // Process input tensors and return preprocessed batch ready for DNN
            //virtual torch::Tensor preprocess(const ITorchTensorSet::pointer& input) = 0;
            virtual std::vector<torch::Tensor> preprocess(const ITorchTensorSet::pointer& input) = 0;
            
            // Get preprocessing metadata (e.g., scaling factors, dimensions)
            virtual Configuration get_metadata() const = 0;
        };
        
    } //namespace WireCell
} //namespace WIRECELL_IDNNROIPREPROCESS_H

#endif