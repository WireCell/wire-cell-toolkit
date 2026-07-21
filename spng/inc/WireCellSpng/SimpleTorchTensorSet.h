#ifndef WIRECELLSPNG_SIMPLETORCHTENSORSET
#define WIRECELLSPNG_SIMPLETORCHTENSORSET

#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell::SPNG {

    class SimpleTorchTensorSet : public WireCell::ITorchTensorSet {
      public:
        explicit SimpleTorchTensorSet(int ident)
            : m_ident(ident)
        {}

        SimpleTorchTensorSet(int ident, Configuration md, ITorchTensor::shared_vector tv)
            : m_ident(ident)
            , m_md(md)
            , m_tv(tv)
        {}

        SimpleTorchTensorSet(int ident, Configuration md, const ITorchTensor::vector& tv)
            : m_ident(ident)
            , m_md(md)
            , m_tv(std::make_shared<ITorchTensor::vector>(tv.begin(), tv.end()))
        {}

        SimpleTorchTensorSet(ITorchTensorSet & input) 
            : m_ident(input.ident()),
              m_md(input.metadata()),
              m_tv(input.tensors())
        {}

        virtual ~SimpleTorchTensorSet() {}

        /// Return some identifier number that is unique to this set.
        virtual int ident() const { return m_ident; }

        /// Optional metadata associated with the tensors
        virtual Configuration metadata() const { return m_md; }

        virtual ITorchTensor::shared_vector tensors() const { return m_tv; }

      private:
        int m_ident;
        Configuration m_md;
        ITorchTensor::shared_vector m_tv;
    };

    class EmptyTorchTensorSet : public WireCell::ITorchTensorSet {
      public:
        
        virtual ~EmptyTorchTensorSet() {}

        virtual int ident() const { return -1; }
        virtual Configuration metadata() const { return Configuration(); }
        virtual ITorchTensor::shared_vector tensors() const { return std::make_shared<ITorchTensor::vector>(); }
    };


}  // namespace WireCell::SPNG

#endif
