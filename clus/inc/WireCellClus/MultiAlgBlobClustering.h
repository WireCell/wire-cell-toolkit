#ifndef WIRECELL_CLUS_MULTIALGBLOBCLUSTERING
#define WIRECELL_CLUS_MULTIALGBLOBCLUSTERING

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellClus/IEnsembleVisitor.h"
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
            // special name "img" dumps "live" before clustering
            std::string name;   // bee type name
            std::string detector;  // bee geom name
            std::string algorithm; // bee alg name, defaults to type
            std::string pcname;    // PC to take
            std::string grouping;  // grouping to take (default to "live")
            std::string visitor;   // if given, dump just after this visitor runs and any cluster ID enumeration
            std::vector<std::string> coords;
            bool individual;
            int filter{1};// 1 for on, 0 for off, -1 for inverse filter  
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
    
        Facade::Grouping& load_grouping(
            Facade::Ensemble& ensemble,
            const std::string& name,
            const std::string& path,
            const ITensorSet::pointer ints);

        std::map<std::string, ApaBeePoints> m_bee_points;

        // New helper function to fill bee points
        void fill_bee_points(const std::string& name, const Facade::Grouping& grouping);
        void fill_bee_points_from_cluster(
            Bee::Points& bpts, const Facade::Cluster& cluster, 
            const std::string& pcname, const std::vector<std::string>& coords,
            int filter);

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

        /** Config: "groupings"
         *
         * List of groupings to select for processing.
         *
         * Default: ["live","dead"].
         */
        std::vector<std::string> m_groupings = {"live","dead"};

        /** Config: "inpath"
         *
         * The BASE datapath for the input pc tree data.  This may be a regular
         * expression which will be applied in a first-match basis against the
         * input tensor datapaths.  If the matched tensor is a pcdataset it is
         * interpreted as providing the nodes dataset.  Otherwise the matched
         * tensor must be a pcgraph.
         *
         * A "%d" will be interpolated with the ident number of the input tensor
         * set.
         * 
         * See also "insubpaths".
         */
        std::string m_inpath{".*"};

        /** Config: "outpath"
         *
         * The BASE datapath for the resulting pc tree data.
         *
         * A "%d" will be interpolated with the ident number of the input tensor
         * set.
         *
         * See outsubpaths.
         */
        std::string m_outpath{""};

        /** Config: insubpaths, outsubpaths.
         *
         * By default, a grouping of a given NAME is located at an input or
         * output path spelled as: "{inpath,outpath}/NAME".
         *
         * If a grouping NAME is found in either insubpath or outsubpath then
         * this default is overridden.  Both parameters are array of objects,
         * each object has keys "name" and "subpath".  The subpath is a simple
         * string suffix and thus should include a leading "/" if the user
         * wishes to locate the grouping in a "subdirectory".
         */
        // See issue #375 and #416.  
        std::map<std::string, std::string> m_insubpaths, m_outsubpaths;
        std::string inpath(const std::string& name, int ident);
        std::string outpath(const std::string& name, int ident);

        /** Config: "perf"
         *
         * If true, emit time/memory performance measures.  Default is false.
         */
        bool m_perf{false};

        /** Config: "grouping2file_prefix"
         *
         * If not empty, dump the final grouping to a file with this prefix + potential RSE info + .npz.
         */

        std::string m_grouping2file_prefix{};

        /** Config: "dump_json"
         *
         * If true, dump to files like {live,dead}-summary-<ident>.json a
         * summary of the groupings.  Default is false.  The dumps are rather
         * large despite omitting point cloud point data.  Use of jq or other
         * tool is expected.  These are intended only for debugging / validating.
         */
        bool m_dump_json{false};

        /** Config: "cluster_id_order"
         *
         * The various operations can lead to redundantly or non-sequentially
         * numbered cluster idents.  The application of merge_clusters() will
         * cause a resulting cluster to have the ID of its first contributing
         * constituent.  The application of separate() will cause all clusters
         * to have the ID of the original.
         *
         * When this parameter is given, the cluster IDs will be reset after
         * each component operation, on a per-grouping basis, so that they
         * represent a specific sort order.
         *
         * - "tree" :: use the as-is child-order which represents
         *   child-insertion order unless some operation has explicitly
         *   reordered the underlying PC tree.
         *
         * - "size" :: use the "size" of the cluster as determined by
         *   cluster_less() to order the cluster IDs.  This considers the
         *   cluster length, number of children, number of points followed by
         *   per-view min then max bounds and finally cluster center.
         *
         * - "" :: cluster IDs are not modified.
         *
         * By default, no ID rewriting is performed.
         *
         * Notes:
         *
         * - When an ID cluster order is applied, the ID counting starts from 1.
         * - The default (unset) ID is -1.
         */
        std::string m_clusters_id_order;

        // configurable parameters for dead-live clustering
        int m_dead_live_overlap_offset{2};

        // Keep track of configured clustering methods with their metadata to
        // assist in debugging/logging.
        struct EnsembleVisitor {
            std::string name;
            IEnsembleVisitor::pointer meth;
        };
        /** Config: pipeline
         *
         * Array of type/name of instances of IEnsembleVisitor to execute in the pipeline.
         */
        std::vector<EnsembleVisitor> m_pipeline;

        // the anode to be processed
        std::vector<IAnodePlane::pointer> m_anodes;

        IDetectorVolumes::pointer m_dv;

        // the face to be processed
        // int m_face{0};

        // the geometry helper
        // IClusGeomHelper::pointer m_geomhelper;
    };
}  // namespace WireCell::Clus

#endif
