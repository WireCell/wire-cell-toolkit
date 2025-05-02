#ifndef WIRECELL_CLUS_MULTIALGBLOBCLUSTERING
#define WIRECELL_CLUS_MULTIALGBLOBCLUSTERING

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellClus/IClusteringMethod.h"
#include "WireCellClus/Facade.h"

#include "WireCellAux/Logger.h"

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ITerminal.h"

#include "WireCellUtil/Bee.h"

#include <vector>

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
        
        // Replace the existing bee points structures with a more flexible approach
        struct BeePointsConfig {
            std::string name;
            std::string detector;
            std::string algorithm;
            std::string pcname;
            std::vector<std::string> coords;
            bool individual;
        };

        // Vector to store configurations for multiple bee points sets
        std::vector<BeePointsConfig> m_bee_points_configs;
        
         // Nested structure to store bee points objects for each configuration, by APA and face
        // First key: bee points set name, second key: "anode_id-face_id" string
        struct ApaBeePoints {
             // Default constructor (add this)
            ApaBeePoints() = default;
            
            // Global points (used when individual == false)
            Bee::Points global;
            
            // Individual points (used when individual == true)
            // Key is "anode_id-face_id" string
            std::map<int, std::map<int , Bee::Points> > by_apa_face; // apa, face
            
        };
    
        std::map<std::string, ApaBeePoints> m_bee_points;

        // New helper function to fill bee points
        void fill_bee_points(const std::string& name, const Facade::Grouping& grouping);
        void fill_bee_points_from_cluster(
            Bee::Points& bpts, const Facade::Cluster& cluster, 
            const std::string& pcname, const std::vector<std::string>& coords);

        void fill_bee_patches_from_grouping(const Facade::Grouping& grouping);
        void fill_bee_patches_from_cluster(const Facade::Cluster& cluster);

        std::map<int, std::map<int, Bee::Patches>> m_bee_dead_patches; 
        // Bee::Patches m_bee_dead; // dead region ...

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
         * The BASE datapath for the input pc tree data.  This may be a
         * regular expression which will be applied in a first-match
         * basis against the input tensor datapaths.  If the matched
         * tensor is a pcdataset it is interpreted as providing the
         * nodes dataset.  Otherwise the matched tensor must be a
         * pcgraph.  See in{live,dead} parameters.
         */
        std::string m_inpath{".*"};

        /** Config: "outpath"
         *
         * The BASE datapath for the resulting pc tree data.  A "%d" will be
         * interpolated with the ident number of the input tensor set.  See
         * out{live,dead} parameters.
         */
        std::string m_outpath{""};

        /** Config: inlive, outlive, indead, outdead.
         *
         * This quad of parameters {in,out}{live,dead} set suffixes to the
         * {in,out}path.  By default, they are "/live" and "/dead" for both
         * in/out pairs.  In exceptional cases, these need setting if the
         * upstream sources or downstream sinks require different patterns.
         */
        // See issue #375.  Note, configure() will expand these to full paths by
        // appending to m_inpath and m_outpath.
        std::string m_inlive{"/live"};
        std::string m_indead{"/dead"};
        std::string m_outlive{"/live"};
        std::string m_outdead{"/dead"};


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

        // Keep track of configured clustering methods with their metadata to
        // assist in debugging/logging.
        struct ClusteringMethod {
            std::string name;
            IClusteringMethod::pointer meth;
            Configuration config;
        };
        std::vector<ClusteringMethod> m_clustering_chain;
        //Configuration m_func_cfgs;

        // the anode to be processed
        std::vector<IAnodePlane::pointer> m_anodes;

        IDetectorVolumes::pointer m_dv;

        // the face to be processed
        int m_face{0};

        // the geometry helper
        // IClusGeomHelper::pointer m_geomhelper;
    };
}  // namespace WireCell::Clus

#endif
