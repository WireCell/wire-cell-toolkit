#ifndef WIRECELL_AUX_SIMPLETENSOR
#define WIRECELL_AUX_SIMPLETENSOR

#include "WireCellIface/ITensor.h"
#include "WireCellUtil/MultiArray.h"
#include <cstring>

namespace WireCell::Aux {


    class SimpleTensor : public WireCell::ITensor {
      public:

        // Create a null tensor.
        SimpleTensor()
            : m_typeinfo(typeid(void))
        {
        }
        // Create a simple tensor with a null array and type, just metadata, if any
        explicit SimpleTensor(const Configuration& md)
            : m_typeinfo(typeid(void))
            , m_mdptr(std::make_unique<Configuration>(md))
        {
        }

        // Create simple tensor, allocating space for data.  If data
        // given it must have at least as many elements as implied by
        // shape and that span will be copied into allocated memory.
        // The SimpleTensor will allocate memory if a null data
        // pointer is given but note the type is required so pass call like:
        //   SimpleTensor(shape, (Type*)nullptr)
        template <typename ElementType>
        SimpleTensor(const shape_t& shape, const ElementType* data)
            : m_typeinfo(typeid(ElementType))
            , m_sizeof(sizeof(ElementType))
            , m_shape(shape)
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
        // Create simple tensor, allocating space for data.  If data
        // given it must have at least as many elements as implied by
        // shape and that span will be copied into allocated memory.
        // The SimpleTensor will allocate memory if a null data
        // pointer is given but note the type is required so pass call like:
        //   SimpleTensor(shape, (Type*)nullptr)
        template <typename ElementType>
        SimpleTensor(const shape_t& shape,
                     const ElementType* data,
                     const Configuration& md)
            : m_typeinfo(typeid(ElementType))
            , m_sizeof(sizeof(ElementType))
            , m_shape(shape)
            , m_mdptr(std::make_unique<Configuration>(md))
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
        Configuration& metadata() {
            if (!m_mdptr) {
                // lazy construction.
                m_mdptr = std::make_unique<Configuration>();
            }
            return *m_mdptr;
        }


        // ITensor const interface.
        virtual const std::type_info& element_type() const { return m_typeinfo; }
        virtual size_t element_size() const { return m_sizeof; }

        virtual shape_t shape() const { return m_shape; }

        virtual const std::byte* data() const { return m_store.data(); }
        virtual size_t size() const { return m_store.size(); }

        virtual Configuration metadata() const { return *m_mdptr; }

      private:
        const std::type_info& m_typeinfo;
        size_t m_sizeof{0};
        std::vector<size_t> m_shape;
        std::vector<std::byte> m_store;
        std::unique_ptr<Configuration> m_mdptr; // avoid constructor if empty
    };

}  // namespace WireCell::Aux

#endif
