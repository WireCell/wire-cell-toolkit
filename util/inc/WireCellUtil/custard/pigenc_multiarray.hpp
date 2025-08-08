/**
   Serialize Boost.Multiarray with pigenc
 */
#ifndef PIGENC_MULTIARRAY_HPP
#define PIGENC_MULTIARRAY_HPP

#ifdef __clang__
#  if defined(__has_warning)
#    define HAS_WARNING(warning) __has_warning(warning)
#  else
#    define HAS_WARNING(warning) 1
#  endif
#else
#  define HAS_WARNING(warning) 1
#endif

#if HAS_WARNING("-Wdeprecated-declarations")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "boost/multi_array/storage_order.hpp"

#if HAS_WARNING("-Wdeprecated-declarations")
#pragma GCC diagnostic pop
#endif


#include "pigenc.hpp"

namespace pigenc::multiarray {

    /// Dump an array into a pigenc header object
    template<typename ArrayType>
    void dump(pigenc::Header& pighead, const ArrayType& array)
    {
        using Scalar = typename ArrayType::element;
        std::vector<size_t> shape(array.num_dimensions(),0);
        for (size_t dim=0; dim<array.num_dimensions(); ++dim) {
            shape[dim] = array.shape()[dim];
        }
        const bool fortran_order =  array.storage_order() == boost::fortran_storage_order();
        pighead.set<Scalar>(shape, fortran_order);
    }

    /// Dump an array into a pigenc file object
    template<typename ArrayType>
    void dump(pigenc::File& pig, const ArrayType& array)
    {
        dump(pig.header(), array);
        auto& dat = pig.data();
        dat.clear();
        char* adat = (char*)array.data();
        dat.insert(dat.begin(), adat, adat+pig.header().data_size());
    }

    /// Write multiarray as a pigenc stream (not custard!)
    template<typename ArrayType>
    std::ostream& write(std::ostream& so, const ArrayType& array)
    {
        pigenc::File pig;
        dump<ArrayType>(pig, array);
        pig.write(so);
        return so;
    }

    // Load data from pigenc file to array, return false on if
    // template type is a mismatch.
    template<typename ArrayType>
    bool load(const pigenc::File& pig, ArrayType& array)
    {
        const pigenc::Header& head = pig.header();

        using Scalar = typename ArrayType::element;
        const Scalar* dat = pig.as_type<Scalar>();
        if (!dat) { return false; }

        const auto& shape = head.shape();
        const size_t ndims = shape.size();
        if (ndims != ArrayType::dimensionality) {
            return false;
        }
        if (head.fortran_order()) {
            if (array.storage_order() != boost::fortran_storage_order()) {
                return false;
           }
        } else if (array.storage_order() != boost::c_storage_order()) {
            return false;
        }
        array = ArrayType(dat, shape, array.storage_order());
        return true;
    }


    /// Read multiarray from pigenc stream.  Stream will fail if
    /// ArrayType does not match dtype.  For a type-discovery
    /// mechanism, use a pigenc::File and then load<T>() above.
    template<typename ArrayType>
    std::istream& read(std::istream& si, ArrayType& array)
    {
        pigenc::File pig;
        pig.read(si);
        if (!si) return si;

        if (load(pig, array)) {
            return si;
        }
        si.setstate(std::ios::failbit);
        return si;
    }
    
    

}

#endif
