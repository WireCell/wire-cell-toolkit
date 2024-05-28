/**
   Wire-Cell Toolkit Stream provides methods to send data to a
   std::ostream (which really should be a boost iostream).

   WCT Stream is intended for binary sinks and sources and thus
   applying to std::cin or std::cout is not recomended.

   The Stream API currently supports array-like data by encoding array
   information into Numpy file format using custard/pigenc.
 */

#ifndef WIRECELLUTIL_STREAM
#define WIRECELLUTIL_STREAM

#include "WireCellUtil/custard/pigenc_eigen.hpp"
#include "WireCellUtil/custard/pigenc_stl.hpp"
#include "WireCellUtil/custard/pigenc_multiarray.hpp"
#include "WireCellUtil/MultiArray.h"

#define CUSTARD_BOOST_USE_MINIZ
#include "WireCellUtil/custard/custard_boost.hpp"

#include <boost/iostreams/filtering_stream.hpp>

#include <string>
#include <vector>

namespace WireCell::Stream {

    /// Fill a filtering input stream with filters and a file source
    /// based on the an input file name.  Filters are determined based
    /// on file name extension(s).  Supported extensions include
    /// .tar, .tar.gz, .tar.bz2, .gz, .bz2.  Any unsupported
    /// extensions will result no filtering and a raw binary file
    /// source.
    using custard::input_filters;

    /// Fill a filtering output stream with filters and a file sink
    /// based on the an ouput file name.  See input_filters() for
    /// supported file types.
    using custard::output_filters;

    /// Note, to use these practically, the ostreams need to end in a
    /// tar, zip or other container filter.

    using shape_t = std::vector<size_t>;

    /// Stream array as type-less bytes with all specified.
    std::ostream& write(std::ostream& so, const std::string& fname,
                        const std::byte* bytes, const shape_t& shape,
                        const std::string& dtype_name);

    /// Stream vector to custard stream.  
    template<typename Type>
    std::ostream& write(std::ostream& so, const std::string& fname, const std::vector<Type>& vec) {
        pigenc::File pig;
        pigenc::stl::dump< std::vector<Type> >(pig, vec);
        size_t fsize = pig.header().file_size();
        custard::write(so, fname, fsize);
        if (!so) return so;
        return pig.write(so);
    }
    
    /// Stream vector from custard stream.
    template<typename Type>
    std::istream& read(std::istream& si, std::string& fname, std::vector<Type>& arr) {
        size_t fsize;
        custard::read(si, fname, fsize);
        if (!si) return si;
        pigenc::File pig;
        pig.read(si);
        if (!si) return si;
        bool ok = pigenc::stl::load< std::vector<Type> >(pig, arr);
        if (!ok) si.setstate(std::ios::failbit);
        return si;
    }

    /// Short hand alias as most "ostream" and "istream" will be this
    /// type in practice.
    using boost::iostreams::filtering_istream;
    using boost::iostreams::filtering_ostream;

    /// Stream Eigen array to custard stream.  Must have tar filter!
    template<typename Scalar, int Rows=Eigen::Dynamic, int Cols=Eigen::Dynamic>
    std::ostream& write(std::ostream& so, const std::string& fname, const Eigen::Array<Scalar,Rows,Cols>& arr) {
        pigenc::File pig;
        using ArrType = Eigen::Array<Scalar,Rows,Cols>;
        pigenc::eigen::dump<ArrType>(pig, arr);
        size_t fsize = pig.header().file_size();
        custard::write(so, fname, fsize);
        if (!so) return so;
        return pig.write(so);
    }
    
    /// Stream Eigen array from custard stream.  Must have tar filter!
    template<typename Scalar, int Rows=Eigen::Dynamic, int Cols=Eigen::Dynamic>
    std::istream& read(std::istream& si, std::string& fname, Eigen::Array<Scalar,Rows,Cols>& arr) {
        size_t fsize;
        custard::read(si, fname, fsize);
        if (!si) return si;
        pigenc::File pig;
        pig.read(si);
        if (!si) return si;
        using ArrType = Eigen::Array<Scalar,Rows,Cols>;
        bool ok = pigenc::eigen::load<ArrType>(pig, arr);
        if (!ok) si.setstate(std::ios::failbit);
        return si;
    }

    // Stream Boost.Multiarray to custard.
    template<typename Scalar, long unsigned int NDim>
    std::ostream& write(std::ostream& so, const std::string& fname, const boost::multi_array<Scalar, NDim>& arr) {
        pigenc::File pig;
        using MArrType = boost::multi_array<Scalar, NDim>;
        pigenc::multiarray::dump< MArrType >(pig, arr);
        size_t fsize = pig.header().file_size();
        custard::write(so, fname, fsize);
        if (!so) return so;
        return pig.write(so);
    }
        
    /// Stream multiarray from custard.
    template<typename Scalar, long unsigned int NDim>
    std::istream& read(std::istream& si, std::string& fname, boost::multi_array<Scalar, NDim>& arr) {
        size_t fsize;
        custard::read(si, fname, fsize);
        if (!si) return si;
        pigenc::File pig;
        pig.read(si);
        if (!si) return si;
        using MArrType = boost::multi_array<Scalar, NDim>;
        bool ok = pigenc::multiarray::load< MArrType >(pig, arr);
        if (!ok) si.setstate(std::ios::failbit);
        return si;
    }


}

#endif
