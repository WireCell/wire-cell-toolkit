#ifndef WIRECELL_ITORCHTENSORSET
#define WIRECELL_ITORCHTENSORSET

#include "ITorchTensor.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    class ITorchTensorSet : public IData<ITorchTensorSet> {
       public:
        virtual ~ITorchTensorSet() {}

        /// Return some identifier number that is unique to this set.
        virtual int ident() const = 0;

        /// Optional metadata associated with the set of tensors
        virtual Configuration metadata() const { return Configuration(); }

        /// Return the tensors in this set.
        virtual ITorchTensor::shared_vector tensors() const = 0;
    };
}  // namespace WireCell

#endif
