#ifndef WIRECELL_POINTCLOUDARRAY
#define WIRECELL_POINTCLOUDARRAY

#include "WireCellUtil/Dtype.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/MultiArray.h"
#include "WireCellUtil/Span.h"

#include <string>
#include <vector>
#include <iterator>

namespace WireCell::PointCloud {


    /** A dense array held in a typeless manner, synergistic with
        ITensor.
    */
    class Array {
        
      public:

        using metadata_t = Configuration;

        /** For scalar, per-point values shape must be (n,) where
            where n is number of points.  N-dimensional, per-point
            data has a shape of (n, d1, ...) where d1 is the size of
            first dimension, etc for any others.
        */
        using shape_t = std::vector<size_t>;

        /** The underlying data may be accessed as a typed, flattened
            span.
         */
        template<typename ElementType>
        using span_t = boost::span<ElementType>;

        /** The underlying data may be accessed as a typed, Boost
            multi array to allow indexing.  This form allows change
            value in the underlying array.
        */
        template<typename ElementType, size_t NumDims>
        using indexed_t = boost::multi_array_ref<ElementType, NumDims>;

        /** The underlying data may be accessed as a typed, Boost
            multi array to allow indexing.  This form is const version
            of above.
        */
        template<typename ElementType, size_t NumDims>
        using const_indexed_t = boost::const_multi_array_ref<ElementType, NumDims>;

        // Default constructor.
        Array();

        // Copy constructor
        Array(const Array& rhs);

        // Move constructor
        Array(Array&& rhs);

        // special case, 1D constructor using vector-like
        template<typename Range>
        explicit Array(const Range& r)
        {
            assign(&*std::begin(r), {r.size()}, false);
        }

        /** Store a point array given in flattened, row-major aka
            C-order.  If share is true then no copy of elements is
            done.  See assure_mutable().
        */
        template<typename ElementType>
        explicit Array(ElementType* elements, shape_t shape, bool share)
        {
            assign(elements, shape, share);
        }
        template<typename ElementType>
        explicit Array(const ElementType* elements, shape_t shape, bool share)
        {
            assign(elements, shape, share);
        }

        template<typename Range>
        explicit Array(const Range& r, shape_t shape, bool share)
        {
            assign(&*std::begin(r), shape, share);
        }

        /// Construct with bytes, shape, element size and whether to share the data or not.
        Array(const std::byte* data, const std::string& dtype, const shape_t& shape);
        Array(std::byte* data, const std::string& dtype, const shape_t& shape, bool share);


        /// Special constructor on initializer list
        template<typename ElementType>
        Array(std::initializer_list<ElementType> il) {
            assign(&*il.begin(), {il.size()}, false);
        }


        // Copy assignment 
        Array& operator=(const Array& rhs);

        // Move assignment 
        Array& operator=(Array&& rhs);

        ~Array();


        /** The various assign() methods will discard any held data and assign new data provided as bytes.
         */

        // Const byte data, will copy.
        void assign(const std::byte* data, const std::string& dtype, const shape_t& shape);

        /// Mutable byte data, shareable.
        void assign(std::byte* data, const std::string& dtype, const shape_t& shape, bool share);

        template<typename ElementType>
        void assign(ElementType* elements, const shape_t& shape, bool share)
        {
            clear();
            m_dtype = WireCell::dtype<ElementType>();
            m_shape = shape;
            size_t nbytes = m_ele_size = sizeof(ElementType);
            for (const auto& s : m_shape) {
                nbytes *= s;
            }
            std::byte* bytes = reinterpret_cast<std::byte*>(elements);
            if (share) {
                m_bytes = span_t<std::byte>(bytes, nbytes);
            }
            else {
                m_store.assign(bytes, bytes+nbytes);
                update_span();
            }
        }

        /// Assign via typed pointer.
        template<typename ElementType>
        void assign(const ElementType* elements, const shape_t& shape, bool share) {
            assign((ElementType*)elements, shape, share);
        }

        /// Assign via iterator.
        template<typename Itr>
        void assign(Itr beg, Itr end, const shape_t& shape, bool share) {
            assign((typename Itr::pointer*)&*beg, std::distance(beg,end), shape, share);
        }

        /// Clear all held data.
        void clear()
        {
            m_store.clear();
            m_dtype = "";
            m_bytes = span_t<std::byte>();
            m_shape.clear();
            m_ele_size = 0;
        }

        /** Assure we are in mutable mode.  If previously sharing user
            data, this results in a copy.
        */
        void assure_mutable()
        {
            if (m_store.empty()) {
                m_store.assign(m_bytes.begin(), m_bytes.end());
                update_span();
            }
        }

        /// The slice() methods return subset of this array starting at given
        /// major axis position and spanning given count major axis elements.

        /// The returned slice may possibly share the underlying array data.
        Array slice(size_t position, size_t count, bool share);
        /// The slice holds a copy.
        Array slice(size_t position, size_t count) const;

        template<typename ElementType, size_t NumDims>
        void check_dims() const {
            if (NumDims != m_shape.size()) {
                raise<ValueError>("ndims mismatch %d != %d",
                                  NumDims, m_shape.size());
            }
        }

        template<typename ElementType>
        void check_size() const {
            if (sizeof(ElementType) != m_ele_size) {
                raise<ValueError>("element size mismatch %d != %d",
                                  sizeof(ElementType), m_ele_size);
            }
        }

        /** Return object that allows type-full indexing of N-D array
            (a Boost const multi array reference).  Note, user must
            assure type is consistent with the originally provided
            data at least in terms of element size.  NumDims must be
            the same as the size of the shape() vector.
            template<typename ElementType, size_t NumDims>
        */
        template<typename ElementType, size_t NumDims>
        indexed_t<ElementType, NumDims> indexed() 
        {
            check_size<ElementType>();
            check_dims<ElementType, NumDims>();

            return indexed_t<ElementType, NumDims>(elements<ElementType>().data(), m_shape);
        }
        template<typename ElementType, size_t NumDims>
        const_indexed_t<ElementType, NumDims> indexed() const
        {
            check_size<ElementType>();
            check_dims<ElementType, NumDims>();
            return const_indexed_t<ElementType, NumDims>(elements<const ElementType>().data(), m_shape);
        }

        /** Return a constant span of flattened array data assuming
            data elements are of given type.
        */
        template<typename ElementType=double>
        span_t<const ElementType> elements() const
        {
            check_size<ElementType>();
            const ElementType* edata = 
                reinterpret_cast<const ElementType*>(m_bytes.data());
            return span_t<const ElementType>(edata, m_bytes.size()/sizeof(ElementType));
        }
        template<typename ElementType=double>
        span_t<ElementType> elements() 
        {
            check_size<ElementType>();
            ElementType* edata = 
                reinterpret_cast<ElementType*>(m_bytes.data());
            return span_t<ElementType>(edata, m_bytes.size()/sizeof(ElementType));
        }

        /** Return element at index as type, no bounds checking is
            done.
         */
        template<typename ElementType=double>
        const ElementType& element(size_t index) const
        {
            const ElementType* edata = reinterpret_cast<const ElementType*>(m_bytes.data());
            return *(edata + index);
        }
        template<typename ElementType=double>
        ElementType& element(size_t index) 
        {
            ElementType* edata = reinterpret_cast<ElementType*>(m_bytes.data());
            return *(edata + index);
        }

        /// Return constant span array as flattened array as bytes.
        span_t<const std::byte> bytes() const
        {
            return m_bytes;
        }

        /// Return an Array like this one but filled with zeros to
        /// given number of elements along major axis.
        Array zeros_like(size_t nmaj);

        /// Append a flat array of bytes.  The number of bytes must be
        /// consistent with the element size and the existing shape.
        void append(const std::byte* data, size_t nbytes);

        /// Append a typed array of objects.  The type must be
        /// consistent in size with current element size and the
        /// number of elements must be consistent with the existing
        /// shape.
        template<typename ElementType>
        void append(const ElementType* data, size_t nele) 
        {
            check_size<ElementType>();

            append(reinterpret_cast<const std::byte*>(data),
                   nele * sizeof(ElementType));
        }

        template<typename Itr>
        void append(Itr beg, Itr end)
        {
            const size_t size = std::distance(beg, end);
            append(&*beg, size);
        }
        
        void append(const Array& arr);

        template<typename Range>
        void append(const Range& r)
        {
            append(std::begin(r), std::end(r));
        }

        /// Return the size of the major axis.
        size_t size_major() const
        {
            if (m_shape.empty()) return 0;
            return m_shape[0];
        }

        /// Return the total number of elements in the typed array.
        /// This is the product of the sizes of each dimension in the
        /// shape.
        size_t num_elements() const;

        shape_t shape() const
        {
            return m_shape;
        }

        std::string dtype() const
        {
            return m_dtype;
        }

        template<typename ElementType>
        bool is_type() const
        {
            if (m_dtype.empty()) {
                return false;
            }
            return WireCell::dtype<ElementType>() == m_dtype;
        }

        size_t element_size() const {
            return m_ele_size;
        }

        metadata_t& metadata() { return m_metadata; }
        const metadata_t& metadata() const { return m_metadata; }
        bool operator==(const Array& other) const;
        bool operator!=(const Array& other) const;

      private:

        shape_t m_shape{};
        size_t m_ele_size{0};
        std::string m_dtype{""};

        // if sharing user data, m_store is empty.
        std::vector<std::byte> m_store{};
        // view of either user's data or our store and through which
        // all access is done.
        span_t<std::byte> m_bytes;

        // Metadata is a passive carry.
        metadata_t m_metadata;

        void update_span() {
            m_bytes = span_t<std::byte>(m_store.data(), m_store.data() + m_store.size());
        }

    };                          // Array

}


#endif
