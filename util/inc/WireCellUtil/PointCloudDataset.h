#ifndef WIRECELL_POINTCLOUDDATASET
#define WIRECELL_POINTCLOUDDATASET

#include "WireCellUtil/PointCloudArray.h"

namespace WireCell::PointCloud {

    /** A set of arrays stored and accessed by their name.  All arrays
        maintain equal number of elements (ie, same size_major()).

    */
    class Dataset {
      public:

        using metadata_t = Configuration;

        // copy constructor
        Dataset(const Dataset& other);
        // move constructor
        Dataset(Dataset&& other);

        // copy assignment 
        Dataset& operator=(const Dataset& other);
        // move assignment
        Dataset& operator=(Dataset&& other);

        // Default constructor
        Dataset();

        // Construct a dataset by copying from a set of named arrays.
        explicit Dataset(const std::map<std::string, Array>& arrays);

        ~Dataset();

        /// Return a dataset that consists of a slice of the arrays in this
        /// dataset with count points starting at position.
        Dataset slice(size_t position, size_t count) const;
        Dataset slice(size_t position, size_t count, bool share);

        /// Return the common number of elements along the major axis
        /// of the arrays.  This returns the size of the first
        /// dimension or zero if the dataset is empty.
        size_t size_major() const;        

        /** Add a new named array to the dataset.

            An exception is thrown if the name is already used or if
            the size of the array does not match others in the dataset.
        */
        void add(const std::string& name, const Array& arr);
        void add(const std::string& name, Array&& arr);


        /** A selection is an ordered subset of arrays in the dataset.

            Order reflects that of the list of names.

            If any requested name is not provided by this dataset then
            an empty collection is returned.
        */
        using name_list_t = std::vector<std::string>;

        using array_ptr = std::shared_ptr<Array>;
        using selection_t = std::vector<array_ptr>;
        selection_t selection(const name_list_t& names);

        using const_array_ptr = std::shared_ptr<Array const>;
        using const_selection_t = std::vector<array_ptr>;
        const_selection_t selection(const name_list_t& names) const;
        
        /** Return named array or nullptr if not found.

            Even in non-const case, the array pointer is merely lent.

        */
        array_ptr get(const std::string& name) const;

        /** Return mutable array or nullptr if not found.

            Caller is warned that changing entries directly through
            this array will not trigger the auto update callbacks.
         */
        using mutable_array_ptr = std::shared_ptr<Array>;
        mutable_array_ptr get(const std::string& name);

        // FIXME, add something that returns something that produces a
        // typed tuple of array elements given their common index.

        /// Return the sorted names of the arrays in this.
        using keys_t = std::vector<std::string>;
        keys_t keys() const;
        bool has(const std::string& key) const;

        /// Return a Dataset with same name and type of arrays as this
        /// one but with arrays filled with zeros to the given number
        /// of elements along major axis.
        Dataset zeros_like(size_t nmaj);

        /// Append arrays in tail to this.  This causes array copies.
        void append(const Dataset& tail);

        /// Append as bare, named arrays.  
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

        bool operator==(Dataset const& other) const;
        bool operator!=(Dataset const& other) const;

        bool empty() const { return m_store.empty(); }
        size_t size() const { return m_store.size(); }

      private:


        /// A dataset is a glorified map from string to array.  
        using store_t = std::map<std::string, mutable_array_ptr>;
        store_t m_store;

        std::vector<append_callback_f> m_append_callbacks;
        
        // Metadata is a passive carry.
        metadata_t m_metadata;

    };                          // Dataset
    
    // A "null" filter which incurs a copy.
    inline Dataset copy(const Dataset& ds) { return ds; }


}



#endif
