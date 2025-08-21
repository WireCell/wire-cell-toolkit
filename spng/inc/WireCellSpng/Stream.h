/**
   IO stream functions for Torch tensors.

   This APIs mimicks WireCellUtil/Stream.h.

   Example:

     #include <boost/iostreams/filtering_stream.hpp>

     std::string aname;
     torch::tensor aten;

   Writing:

     boost::iostreams::filtering_ostream m_out;
     output_filters(m_out, m_outname);
     if (m_out.size() < 1) { error(); }
     write(m_out, aname, aten);
     m_out.flush();
     m_out.pop();

   Reading:

     input_filters(m_in, m_inname);
     if (m_in.size() < 1) { error(); }
     read(m_in, aname, aten);
     if (!m_in) { error(); }

 */

#ifndef WIRECELLTORCHSTREAM
#define WIRECELLTORCHSTREAM

#include <vector>
#include "WireCellSpng/Torch.h"

// Same namespace as WireCellUtil/Stream.h.  Functions are overloaded.
namespace WireCell::Stream {
    
    // Extract dtype as pigenc dtype string.
    std::string dtype_of(torch::Tensor ten);
    // Convert pigenc dtype string to TensorOptions
    torch::TensorOptions tensoroptions_from_dtype(const std::string& dtype);

    // Tensor::sizes() returns an IntArrayRef which uses long int instead of
    // size.  This does the type conversion.
    using shape_t = std::vector<size_t>;
    shape_t shape_of(torch::Tensor ten);

    std::ostream& write(std::ostream& so, const std::string& fname, torch::Tensor ten);
    std::istream& read(std::istream& si, std::string& fname, torch::Tensor& ten);

}

#endif
