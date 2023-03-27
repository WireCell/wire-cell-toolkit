/**
   Utilities support the tensor data model.
 */
#ifndef WIRECELLAUX_TENSORDM
#define WIRECELLAUX_TENSORDM

#include "WireCellIface/ITensorSet.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/PointCloud.h"
#include "WireCellUtil/PointGraph.h"
#include "WireCellUtil/Waveform.h"

namespace WireCell::Aux::TensorDM {

    /// Generic support

    /// Build a tensor set from set of tensors.  The tensor data model
    /// makes no requirements on the ITensorSet itself.  The user may
    /// combine tensors representing multiple objects into one set.
    ITensorSet::pointer as_tensorset(const ITensor::vector& tens,
                                     int ident = 0,
                                     const Configuration& tsetmd = Json::objectValue);


    using located_t = std::map<std::string, ITensor::pointer>;

    /// Index the tensors by their datapaths.
    located_t index_datapaths(const ITensor::vector& tens);

    /// Return first of type
    ITensor::pointer first_of(const ITensor::vector& tens,
                              const std::string& datatype);

    /// Return all tensors with a datapath matching regex pattern.
    ITensor::vector match_at(const ITensor::vector& tens,
                             const std::string& datapath);

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
                                   const std::string datapath="",
                                   bool share=false);

    /// Convenience function calling above using the tensors from the
    /// ITensorSet.  Note the dataset level metadata is taken from the
    /// ITensor with datatype "pcdataset" and not from ITensorSet.
    PointCloud::Dataset as_dataset(const ITensorSet::pointer& its,
                                   const std::string datapath="",
                                   bool share=false);


    /// PointGraph support

    /**
       Convert tensors representing point cloud graph to PointGraph.

       The tensor found at datapath is converted.  If datapath is
       empty the first pcgraph found is used.
     */
    PointGraph as_pointgraph(const ITensor::vector& tens,
                             const std::string& datapath="");

    /** Convert a point graph to tensors.

        The vector first contains the conversion of the "nodes"
        dataset with datapaths <datapath>/nodes/ followed by the
        conversion of the "edges" dataset with datapaths
        <datapath>/edges/.
    */
    ITensor::vector as_tensors(const PointGraph& pcgraph,
                               const std::string datapath);


    /// Frame support


    /** Convert channel masks to JSON rep. */
    Configuration as_config(const Waveform::ChannelMasks& cms);

    /** Convert channel mask map to JSON rep. */
    Configuration as_config(const Waveform::ChannelMaskMap& cmm);

    /** Convert JSON to CMM */
    Waveform::ChannelMaskMap as_cmm(const Configuration& jcmm);


    /// Default waveform sample (non)transform function.
    inline float identity(float x) { return x; }

    /** Convert a single trace to a trace tensor. */
    ITensor::pointer as_trace_tensor(ITrace::pointer trace,
                                     const std::string& datapath, 
                                     bool truncate=false,
                                     std::function<float(float)> transform=identity);

    /** Convert a sequence of traces to a trace tensor. */
    ITensor::pointer as_trace_tensor(const ITrace::vector& traces,
                                     const std::string& datapath,
                                     bool truncate=false,
                                     std::function<float(float)> transform=identity);

    /** Frame tensor mode.
        
    Frames may be mapped to trace tensors in three ways.
    
    - sparse :: Convert each trace to a 1D trace tensor.  This is the
                most lossless
    
    - unified :: Combine all traces into a 2D trace tensor.  This is
                 lossy.  Overlapping trace extent and individual tbin
                 are both lost.  Zero padding is applied.

    - tagged :: Only consider tagged traces, each set of which is
                stored in a 2D trace tensor.  Each are lossy in ways
                common with unified.  In addition, any trace not
                tagged is dropped.
    */
    enum struct FrameTensorMode { sparse, tagged, unified };

    /**
       Convert an IFrame to tensors.

       The tensor representing the frame aggregate will be the first
       element.

       Trace tensors follow in order reflecting total trace ordering.

       Last come trace data tensors.  These are ordered depending on
       mode and trace tag ordering and come in blocks of pcdataset.

       Note, this ordering is only convenience and the definitive
       ordering is set by the attributes traces and tracedata in the
       frame tensor MD.

       The "datapath" specifies the datapath of the frame tensor and
       is used as the root for the datapath of all other tensors.
       Typically, it should be made unique across all calls.  A good
       choice is "frames/<frame-ident>".

       If truncate is true the trace tensors are of type "i2" else
       "f4".  Prior to truncation, the transform is applied to trace
       array values.
    */
    ITensor::vector as_tensors(IFrame::pointer frame,
                               const std::string& datapath,
                               FrameTensorMode mode = FrameTensorMode::sparse,
                               bool truncate=false,
                               std::function<float(float)> transform=identity);

    /**
      Convert an tensors representing frame to IFrame.

      This is the inverse of as_tensors() above but does not require
      the same ordering.

      The datapath specifies the frame tensor which aggregates the
      others.  If empty then the first tensor of type "frame" is
      assumed.

      If transform is given it will be applied to trace samples after
      any inflation to float is performed.  This transform is
      conceptually the inverse of whatever may have been provided to
      as_tensors() above.

    */
    IFrame::pointer as_frame(const ITensor::vector& tens,
                             const std::string& datapath="",
                             std::function<float(float)> transform=identity);


}

#endif
