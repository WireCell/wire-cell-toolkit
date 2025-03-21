#ifndef WIRECELL_CLUS_MULTIALGBLOBCLUSTERING
#define WIRECELL_CLUS_MULTIALGBLOBCLUSTERING

#include "WireCellAux/Logger.h"
#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellUtil/Bee.h"

namespace WireCell::Clus {

    class MultiAlgBlobClustering
        : public Aux::Logger, public ITensorSetFilter, public IConfigurable, public ITerminal
    {
       public:
        MultiAlgBlobClustering();
        virtual ~MultiAlgBlobClustering() = default;

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual void finalize();

       private:

        /** Config: bee_zip
         *
         * The path to a zip file to fill with Bee JSON
         */
        Bee::Sink m_sink;
        int m_last_ident{-1};
        int m_initial_index{0};  // Default to 0 for backward compatibility
        Bee::Points m_bee_img, m_bee_ld;
        Bee::Patches m_bee_dead;
        // Add new member variables for run/subrun/event
        int m_runNo{0};
        int m_subRunNo{0};
        int m_eventNo{0};
        bool m_use_config_rse{false};  // Flag to determine if we use configured RSE

        void flush(int ident = -1);
        void flush(WireCell::Bee::Points& bpts, int ident);

        bool m_save_deadarea{false};


        // Count how many times we are called
        size_t m_count{0};

        /** Config: "inpath"
         *
         * The datapath for the input point graph data.  This may be a
         * regular expression which will be applied in a first-match
         * basis against the input tensor datapaths.  If the matched
         * tensor is a pcdataset it is interpreted as providing the
         * nodes dataset.  Otherwise the matched tensor must be a
         * pcgraph.
         */
        std::string m_inpath{".*"};

        /** Config: "outpath"
         *
         * The datapath for the resulting pcdataset.  A "%d" will be
         * interpolated with the ident number of the input tensor set.
         */
        std::string m_outpath{""};

        /** Config: "perf"
         *
         * If true, emit time/memory performance measures.  Default is false.
         */
        bool m_perf{false};

        /** Config: "dump_json"
         *
         * If true, dump to files like {live,dead}-summary-<ident>.json a
         * summary of the groupings.  Default is false.  The dumps are rather
         * large despite omitting point cloud point data.  Use of jq or other
         * tool is expected.  These are intended only for debugging / validating.
         */
        bool m_dump_json{false};

        // configurable parameters for dead-live clustering
        int m_dead_live_overlap_offset{2};

        // clustering_examine_x_boundary
        // double m_x_boundary_low_limit{-1*units::cm};
        // double m_x_boundary_high_limit{257*units::cm};

        Configuration m_func_cfgs;

        // the anode to be processed
        std::vector<IAnodePlane::pointer> m_anodes;

        IDetectorVolumes::pointer m_dv;

        // the face to be processed
        int m_face{0};

        // the geometry helper
        IClusGeomHelper::pointer m_geomhelper;
    };
}  // namespace WireCell::Clus

#endif
