/**
   Utilities support the tensor data model for Point Cloud Array and Dataset.
 */
#ifndef WIRECELLAUX_TENSORDMDATASET
#define WIRECELLAUX_TENSORDMDATASET

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellIface/ITensorSet.h"

#include <unordered_map>

namespace WireCell::Aux::TensorDM {

    /// PointCloud support

    /// Convert Array to an ITensor following the WCT tensor data
    /// model for "pcarray".  ValueError exception is thrown if array
    /// type is not supported.
    ITensor::pointer as_tensor(const PointCloud::Array& array,
                               const std::string& datapath);

    /// Convert a Dataset to a vector of ITensor following the WCT
    /// tensor data model.  The ITensor representing the "pcdataset"
    /// is element 0.  All other ITensor representing the dataset
    /// arrays have a datapath formed by appending /arrays/<count> to
    /// the given datapath.  
    ITensor::vector as_tensors(const PointCloud::Dataset& dataset,
                               const std::string& datapath);

    /// Convert an ITensor of datatype "pcarray" to an point cloud
    /// array.  A ValueError exception is thrown if the element type
    /// of the tensor is not supported.  If share is true, an unsafe
    /// zero-copy optimization is performed by having the returned
    /// Array share the same array memory as held by the ITensor.  The
    /// caller must assure the ITensor remains alalive for the
    /// lifetime of the Array in order to avoid memory access errors.
    PointCloud::Array as_array(const ITensor::pointer& ten, bool share=false);

    /// Convert a collection of ITensor following the WCT tensor data
    /// model of "pcdataset" to a Dataset.  See as_array() for meaning
    /// of "share".  The tensor found at the given dpath is used else
    /// the first found with metadata attribute "datatype" of
    /// "pcdataset" is taken as the kernel of aggregation.  A
    /// ValueError is thrown if tensors are badly formed.
    PointCloud::Dataset as_dataset(const ITensor::vector& tensors,
                                   const std::string& datapath="",
                                   bool share=false);
    PointCloud::Dataset as_dataset(const TensorIndex& tensors,
                                   const std::string& datapath="",
                                   bool share=false);
    PointCloud::Dataset as_dataset(const TensorIndex& tensors,
                                   ITensor::pointer top,
                                   bool share=false);

    /// Convenience function calling above using the tensors from the
    /// ITensorSet.  Note the dataset level metadata is taken from the
    /// ITensor with datatype "pcdataset" and not from ITensorSet.
    PointCloud::Dataset as_dataset(const ITensorSet::pointer& its,
                                   const std::string& datapath="",
                                   bool share=false);



}

#endif
