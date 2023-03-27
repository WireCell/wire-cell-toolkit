#ifndef WIRECELL_AUX_SIMPLETENSOR
#define WIRECELL_AUX_SIMPLETENSOR

#include "WireCellIface/ITensor.h"

#include <boost/multi_array.hpp>
#include <cstring>

namespace WireCell::Aux {


    class SimpleTensor : public WireCell::ITensor {
      public:

        // Create a simple tensor with a null array and type, just metadata, if any
        explicit SimpleTensor(const Configuration& md = Configuration())
            : m_md(md)
            , m_typeinfo(typeid(void))
        {
        }

        // Create simple tensor, allocating space for data.  If data
        // given it must have at least as many elements as implied by
        // shape and that span will be copied into allocated memory.
        // The SimpleTensor will allocate memory if a null data
        // pointer is given but note the type is required so pass call like:
        //   SimpleTensor(shape, (Type*)nullptr)
        template <typename ElementType>
        SimpleTensor(const shape_t& shape,
                     const ElementType* data,
                     const Configuration& md = Configuration())
            : m_shape(shape)
            , m_md(md)
            , m_typeinfo(typeid(ElementType))
            , m_sizeof(sizeof(ElementType))
        {
            size_t nbytes = element_size();
            for (const auto& s : m_shape) {
                nbytes *= s;
            }
            if (data) {
                const std::byte* bytes = reinterpret_cast<const std::byte*>(data);
                m_store.assign(bytes, bytes+nbytes);
            }
            else {
                m_store.resize(nbytes);
            }
        }
        virtual ~SimpleTensor() {}

        /** Creator may use the underlying data store allocated in
         * contructor in a non-const manner to set the elements.

         Eg, using boost::multi_array_ref:

         SimpleTensor<float> tens({3,4,5});
         auto& d = tens.store();
         boost::multi_array_ref<float, 3> ma(d.data(), {3,4,5});
         md[1][2][3] = 42.0;
        */
        std::vector<std::byte>& store() { return m_store; }
        Configuration& metadata() { return m_md; }

        // ITensor const interface.
        virtual const std::type_info& element_type() const { return m_typeinfo; }
        virtual size_t element_size() const { return m_sizeof; }

        virtual shape_t shape() const { return m_shape; }

        virtual const std::byte* data() const { return m_store.data(); }
        virtual size_t size() const { return m_store.size(); }

        virtual Configuration metadata() const { return m_md; }

      private:
        std::vector<std::byte> m_store;
        std::vector<size_t> m_shape;
        Configuration m_md;
        const std::type_info& m_typeinfo;
        size_t m_sizeof{0};
    };

}  // namespace WireCell::Aux

#endif
