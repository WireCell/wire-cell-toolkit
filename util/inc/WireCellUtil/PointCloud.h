#ifndef WIRECELL_POINTCLOUD
#define WIRECELL_POINTCLOUD

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Dtype.h"
#include "WireCellUtil/Configuration.h"

// fixme: Only need to keep this ifdef in place unil users upgrade to
// at least boost 1.78.  Can then remove this test from wscript.
#include "WireCellUtil/BuildConfig.h"
#ifdef HAVE_BOOST_CORE_SPAN_HPP
#include "boost/core/span.hpp"
#else
#include "WireCellUtil/boost/core/span.hpp"
#endif  

#include "boost/multi_array.hpp"
#include <set>
#include <map>
#include <string>
#include <vector>
#include <iterator>

namespace WireCell::PointCloud {


    /** The Array and Dataset passively carry metadata */
    using metadata_t = Configuration;

    /** A dense array held in a typeless manner, synergistic with
        ITensor.
    */
    class Array {
        
      public:

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

        // special case, 1D constructor using vector-like
        template<typename Range>
        explicit Array(const Range& r)
        {
            assign(&*std::begin(r), {r.size()}, false);
        }

        // Want defaults for all the rest.
        Array() = default;

        // Copy constructor
        Array(const Array& rhs);

        // Assignment 
        Array& operator=(const Array& rhs);

        //  Move
        Array(Array&&) = default;
        Array& operator=(Array&&) = default;

        ~Array() = default;

        /// Special constructor on initializer list
        template<typename ElementType>
        Array(std::initializer_list<ElementType> il) {
            assign(&*il.begin(), {il.size()}, false);
        }

        /** Discard any held data and assign new data.  See
            constructor for arguments.
        */
        template<typename ElementType>
        void assign(ElementType* elements, shape_t shape, bool share)
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
        void assign(const ElementType* elements, shape_t shape, bool share) {
            assign((ElementType*)elements, shape, share);
        }

        /// Assign via iterator.
        template<typename Itr>
        void assign(Itr beg, Itr end, shape_t shape, bool share) {
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
            if (sizeof(ElementType) != m_ele_size) {
                THROW(ValueError() << errmsg{"element size mismatch in indexed"});
            }
            if (NumDims != m_shape.size()) {
                THROW(ValueError() << errmsg{"ndims mismatch in indexed"});
            }

            return indexed_t<ElementType, NumDims>(elements<ElementType>().data(), m_shape);
        }
        template<typename ElementType, size_t NumDims>
        const_indexed_t<ElementType, NumDims> indexed() const
        {
            if (sizeof(ElementType) != m_ele_size) {
                THROW(ValueError() << errmsg{"element size mismatch in indexed"});
            }
            if (NumDims != m_shape.size()) {
                THROW(ValueError() << errmsg{"ndims mismatch in indexed"});
            }

            return const_indexed_t<ElementType, NumDims>(elements<const ElementType>().data(), m_shape);
        }

        /** Return a constant span of flattened array data assuming
            data elements are of given type.
        */
        template<typename ElementType>
        span_t<const ElementType> elements() const
        {
            if (sizeof(ElementType) != m_ele_size) {
                THROW(ValueError() << errmsg{"element size mismatch in elements"});
            }
            const ElementType* edata = 
                reinterpret_cast<const ElementType*>(m_bytes.data());
            return span_t<const ElementType>(edata, m_bytes.size()/sizeof(ElementType));
        }
        template<typename ElementType>
        span_t<ElementType> elements() 
        {
            if (sizeof(ElementType) != m_ele_size) {
                THROW(ValueError() << errmsg{"element size mismatch in elements"});
            }
            ElementType* edata = 
                reinterpret_cast<ElementType*>(m_bytes.data());
            return span_t<ElementType>(edata, m_bytes.size()/sizeof(ElementType));
        }

        /** Return element at index as type, no bounds checking is
            done.
         */
        template<typename ElementType>
        ElementType element(size_t index) const
        {
            const ElementType* edata = reinterpret_cast<const ElementType*>(m_bytes.data());
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
            if (sizeof(ElementType) != m_ele_size) {
                THROW(ValueError() << errmsg{"element size mismatch in append"});
            }

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

        metadata_t& metadata() { return m_metadata; }
        const metadata_t& metadata() const { return m_metadata; }

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

    };

    /** A collection of indices of points in a point cloud */
    using indices_t = std::vector<size_t>;

    /** An array reference. */
    using ArrayRef = std::reference_wrapper<const Array>;

    /** A list of names, such as used to form a selection. */
    using name_list_t = std::vector<std::string>;

    /** A selection of arrays.  Arrays are held by ref so that any
        updates are available to the user of the selection.
    */
    using selection_t = std::vector<ArrayRef>;

    /** A set of arrays accessed by their name.  All arrays must be
        kept to the same number of elements (points).
    */
    class Dataset {
      public:

        /// A dataset is a glorified map from string to array.
        using store_t = std::map<std::string, Array>;
        using item_t = std::pair<std::string, Array>;

        Dataset() = default;

        // Construct a dataset using an existing store.
        explicit Dataset(const store_t& store) {
            for (const auto& [n,a] : store) {
                add(n,a);
            }
        }

        /// Return the common number of elements along the major axis
        /// of the arrays.  This returns the size of the first
        /// dimension or zero if the dataset is empty.
        size_t size_major() const {
            if (m_store.empty()) {
                return 0;
            }
            return m_store.begin()->second.size_major();
        }
        
        /** Add the array to the dataset.  Return true if its name is
            new.  Throw if number of elments of array does not match
            others in the dataset.  This will copy the array.
        */
        bool add(const std::string& name, const Array& arr);

        /** Add the array to the dataset.  Return true if its name is
            new.  Throw if number of elments of array does not match
            others in the dataset.  This will move the array.
        */
        // bool add(const std::string& name, Array&& arr);

        /// Access the underlying store
        const store_t& store() const { return m_store; }

        /** Return a vector of "references wrappers" around the arrays
            matching the names.  If any name is not provided by this
            then an empty collection is returned.  References may be
            unusual to some so here are two usage examples:

            sel = ds.selection({"one","two"});
            const Array& one = sel[0];
            size_t num_in_two = sel[1].get().num_elements();
        */
        selection_t selection(const name_list_t& names) const;

        /** Return named array or empty if not found */
        const Array& get(const std::string& name) const {
            auto it = m_store.find(name);
            if (it == m_store.end()) {
                static Array empty;
                return empty;
            }
            return it->second;            
        }
        Array& get(const std::string& name) {
            auto it = m_store.find(name);
            if (it == m_store.end()) {
                static Array empty;
                return empty;
            }
            return it->second;            
        }


        // FIXME, add something that returns something that produces a
        // typed tuple of array elements given their common index.

        /// Return the names of the arrays in this.
        using keys_t = std::set<std::string>;
        keys_t keys() const {
            keys_t ret;
            for (const auto& [k,v] : m_store) {
                ret.insert(k);
            }
            return ret;
        }
        bool has(const std::string& key) const {
            for (const auto& [k,v] : m_store) {
                if (k == key) return true;
            }
            return false;
        }

        /// Return a Dataset like this one but filled with zeros to
        /// given number of elements along major axis.
        Dataset zeros_like(size_t nmaj);

        /// Append arrays in tail to this.
        void append(const Dataset& tail) {
            append(tail.store());
        }
        
        /// Append arrays in map to this.
        void append(const std::map<std::string, Array>& tail);

        /** Append zeros to all existing arrays adding nmaj elements
            to the major axis.  This is essentially calling append()
            with the result of zeros_like(nmaj) and thus it triggers
            the update hook.
        */
        void extend(size_t nmaj);

        /** Register a function to be called after a successful append
            is performed on this dataset.  Function is called with the
            right-open range.  Ie, "end" is one past.  NOTE: this is
            different than nanoflan which expects a closed range. 
        */
        using append_callback_f = std::function<void(size_t beg, size_t end)>;
        void register_append(append_callback_f ac)
        {
            m_append_callbacks.push_back(ac);
        }

        metadata_t& metadata() { return m_metadata; }
        const metadata_t& metadata() const { return m_metadata; }

      private:

        store_t m_store;
        std::vector<append_callback_f> m_append_callbacks;
        
        // Metadata is a passive carry.
        metadata_t m_metadata;

    };                          // Dataset
    


}

#endif
