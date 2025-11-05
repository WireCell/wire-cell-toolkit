/**
   A source of ITensorSet providing a pc-tree with cluster, flash and
   cluster-flash information read from Uboone ROOT TTrees and built upon blobs
   read in by UbooneBlobSources.  This reads TTrees from the same file(s) as
   used to produce its inputs:

   - TC
   - T_flash
   - T_match1

   Despite the name, this is a function node.  It must have a preceding subgraph
   like:

   {4 UbooneBlobSource covering uvw/uv/vw/wu} -> BlobSetMerge -> UbooneClusterSource

   It must read the exact same ROOT files as the preceding UbooneBlobSources.

   Any other subgraph is not expected to produce correct results.

   In addition to producing a grouping/cluster/blob style PC-tree equivalent to
   what WCT clustering outputs, this (optionally) adds PCs "light", "flash" and
   "flashlight" to the grouping node and adds a "flash ID" index into "flash"
   array to the cluster nodes in a scalar PC.

   For live case it will make the CTPC.  It does this by ignoring the activity
   from the slices from the input blobs and instead goes back to the
   corresponding entries in the ROOT tree.  This allows slices that may lack any
   blob to still be represented in the CTPC.

 */

#ifndef WIRECELLROOT_UBOONECLUSTERSETSOURCE
#define WIRECELLROOT_UBOONECLUSTERSETSOURCE

#include "WireCellRoot/UbooneTTrees.h" // for UbooneTFiles

#include "WireCellAux/Logger.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IBlobTensoring.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/PointTree.h"

namespace WireCell::Root {

    class UbooneClusterSource : public Aux::Logger,
                                public IBlobTensoring,
                                public IConfigurable
    {
    public:
        UbooneClusterSource();
        virtual ~UbooneClusterSource();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        
        virtual bool operator()(const WireCell::IBlobSet::pointer& in, output_queue& outq);

    private:

        /** Configuration: anode

            Name the IAnodePlane component describing microboone.

            Required for looking up WCT channels given WCP wire indices (which
            for MB are identical to WCT wire-in-plane numbers).
        */
        IAnodePlane::pointer m_anode;

        /** Configuration: "input" (required)

            A string or array of string giving name(s) of file(s) to read for
            input.  Any missing files will be skipped with a log at warn level.

            This configuration MUST BE IDENTICAL to what is given to the four
            UbooneBlobSource's that precede an instance of this node.
        */
        std::unique_ptr<UbooneTFiles> m_files;

        /** Configuration: kind

            A string "live" or "dead" describing what kind of blobs and slices to produce.

            Both require the Trun and TC TTrees.  In addition, "dead" requires
            the TDC TTree.  If the TTree named T_bad_ch exists it will also be
            loaded for both "live" and "dead".

            This configuration must match that given to the upstream UbooneBlobSources.
        */

        /** Configurations: "light", "flash", "flashlight" (optional)

            The names for the PCs to hold optical data of corresponding type.

            If any are unset or set to the empty string (default) no info is stored.

            Note: "flash" references indices into "light" and "flashlight" into
            both.  Omitting the referred arrays causes these references to
            dangle.
        */
        std::string m_light_name, m_flash_name, m_flashlight_name;

        
        /** Configuration: "sampler" (optional)

            Name an IBlobSampler for producing the "3d" point cloud.

            If not given, blob in PC-tree nodes will have no point clouds.

            If given, this configuration must match that given to the upstream
            UbooneBlobSources.
        */
        IBlobSampler::pointer m_sampler;

        /** Configuration: "datapath" (optional)

            Give a tensor data model path for the resulting ITensorSet
            representation.  If the datapath contains "%d" it will be
            interpolated against the "ident" (execution count of this node).

        */
        std::string m_datapath = "pointtrees/%d/uboone";


        /** Configuration: "keep_slices" (default false)

            If true, the original slices shared by the IBlobs in the live case
            are left untouched.

            Otherwise, the ISlices are replaced with "fresh" ones (re)loaded
            from the ROOT file.  This may result in some ISlices in the output
            ICluster that are not connected to any IBlob.
        */

    private: // methods

        IChannel::pointer get_channel(int chanid) const;

        void extract_live(const IBlobSet::vector& blobsets,
                          ISlice::vector& out_slices, IBlob::vector& out_blobs) const;
        void extract_dead(const IBlobSet::vector& blobsets,
                          ISlice::vector& out_slices, IBlob::vector& out_blobs) const;
        void apply_tagger_flags(WireCell::PointCloud::Tree::Points::node_t* cnode, int cluster_id) const;

    private: // data

        // for logging
        size_t m_calls{0};

        // Collect blob sets
        IBlobSet::vector m_cache;

        // flush graph to output queue
        bool flush(output_queue& outq);

        // Return true if newbs is from a new frame
        bool new_frame(const input_pointer& newbs) const;

        double m_drift_speed{1.101 * units::mm / units::us};
        double m_time_offset{-1600 * units::us + 6 * units::mm/(1.101 * units::mm / units::us)};

        std::vector<double> m_angles{1.0472/*60 degrees*/, -1.0472/*-60 degrees*/, 0.0};

    private: // Constants.

        // The microboone sample period.
        const double m_tick{0.5 * units::us};
        // The value and uncertainty to mark "dead" measures
        const ISlice::value_t m_bodge{0, 1e12};

    };
}

#endif
