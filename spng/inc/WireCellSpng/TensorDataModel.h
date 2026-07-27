/** This header support functions for producing and operating on tensors that
 * adhere to the data model (TDM).  See TensorIndex.h and docs/datamodel.org for
 * more details.
 */

#ifndef WIRECELL_SPNG_TENSORDATAMODEL
#define WIRECELL_SPNG_TENSORDATAMODEL

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellIface/IFrame.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell::SPNG {

    /// Return metadata attribute or a default value.
    template <typename T>
    T get(ITorchTensor::pointer ten, const std::string& dotpath, const T& def = T())
    {
        return WireCell::get<T>(ten->metadata(), dotpath, def);
    }

    // standard low-level TDM metadata accessors

    /// Return the datapath or ""
    std::string datapath(ITorchTensor::pointer ten) {
        return get<std::string>(ten, "datapath", "");
    }

    /// Return the datatype or ""
    std::string datatype(ITorchTensor::pointer ten) {
        return get<std::string>(ten, "datatype", "");
    }

    /// Return the parent or ""
    std::string parent(ITorchTensor::pointer ten) {
        return get<std::string>(ten, "parent", "");
    }

    /// Return the number of batches or 0.
    size_t batches(ITorchTensor::pointer ten) {
        // json doesn't know "size_t"
        return get<int>(ten, "batches", 0);
    }


    /// Return tensors of matching datapath regex string
    ITorchTensor::vector find_tensors(const ITorchTensorSet::pointer& ts,
                                      const std::string& datapath_regex);


    /// Make a minimally compliant tensor.
    ITorchTensor::pointer make_tensor(const std::string& datatype,
                                      const std::string& datapath,
                                      Configuration metadata = {},
                                      torch::Tensor ten = {},
                                      const std::string& parent = "",
                                      size_t batches = 0);

    

    /// Functions to make frame related tensors.
    ///
    /// These function are expected to take a COMMON base datapath and will
    /// generate a per-tensor datapath following a fixed schema.  Eg:
    ///
    /// @code
    /// <base>/frame
    /// <base>/traces/<tag>
    /// <base>/chids/<tag>
    /// <base>/summaries/<tag>
    /// <base>/chmasks/<label>
    /// @endcode
    ///
    /// If you try hard to work-around this convention, you MUST change the
    /// "parent" metadata of the produced tensors so that it points your
    /// non-standard datapath for the frame.
    
    /// Make a tensor representing the non-constituent parts of a frame.
    ///
    /// See functions below to make constituents and see make_frame_set() for
    /// a do-it-all function converting an IFrame.
    ///
    /// Be sure to observe parentage ordering (parents first) when adding this
    /// and other pointers to an ITorchTensorSet.
    ///
    /// Reminder: TDM says no array part in frame.
    ITorchTensor::pointer make_frame(const std::string& base_datapath,
                                     int ident, double time, double period);

    /// Make a traces tensor.
    ITorchTensor::pointer make_traces(const std::string& base_datapath,
                                      const std::string& tag,
                                      torch::Tensor traces,
                                      int tbin = 0);
    /// Make a chids tensor
    ITorchTensor::pointer make_chids(const std::string& base_datapath,
                                     const std::string& tag,
                                     torch::Tensor chids);

    /// Make a summaries tensor
    ITorchTensor::pointer make_summaries(const std::string& base_datapath,
                                         const std::string& tag,
                                         torch::Tensor summaries);
    
    /// Make a channel masks tensor
    ITorchTensor::pointer make_chmasks(const std::string& base_datapath,
                                       const std::string& label,
                                       torch::Tensor ten);


    /// Make a set of tensors rooted on the base datapath from the given IFrame.
    ITorchTensor::vector make_frame_set(const std::string& base_datapath, IFrame::pointer iframe);

}

#endif
