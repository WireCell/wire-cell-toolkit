#ifndef WIRECELL_CLUS_FIDUCIALUTILS
#define WIRECELL_CLUS_FIDUCIALUTILS

#include "WireCellClus/ClusteringFuncsMixins.h"

// headers for "static" data sources
#include "WireCellClus/IPCTransform.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"

// headers for "dynamic" data sources
#include "WireCellClus/Facade_Grouping.h"

// Other data types
#include "WireCellUtil/Point.h"


namespace WireCell::Clus {

    
    /** FiducialUtils is a helper for answering data-dependent questions
     * involving regions of space.
     *
     * This class holds various methods to answer questions that requires a mix
     * of static data (eg wire geometry) and the dynamic data being processed
     * (eg clusters, PCs).
     *
     * It is intended to be constructed with static information and after
     * construction it is to be fed dynamic data at least once.  It may perform
     * initial computation on that dynamic data and store the results in order
     * to optimize the queries performed by its methods.  
     *    
     * See also MakeFiducialUtils which is an IEnsembleVisitor that will perform
     * create a FiducialUtils and add it to a grouping.
     *
     * Note: FiducialUtils is approximately equivalent to WCP's ToyFiducial class.
     */

    class FiducialUtils {
    public:

        /// Bundle the "static" data.  This is copyable. 
        struct StaticData {
            IDetectorVolumes::pointer dv;
            IFiducial::pointer fiducial;
            Clus::IPCTransformSet::pointer pcts;
        };

        /// Bundle the "dynamic" data.
        struct DynamicData {
            const Facade::Grouping& live;
            const Facade::Grouping& dead;
        };

        FiducialUtils();

        // Create it with whatever "static" data sources are required.  Add to
        // this argument list as needed as the methods are populated.
        FiducialUtils(StaticData sd);

        // It is possible to create empty and feed static data after
        // construction.  This will re-initialize any dynamic data as well.
        void feed_static(StaticData sd);

        // This must be called at least once to provide the "dynamic" data.  It
        // may be called multiple times.  Each call will start by clearing and
        // rebuilding derived internal data on which the methods operate.
        void feed_dynamic(const DynamicData& dd);


        // Query methods

        bool inside_fiducial_volume(const Point& p, double offset_x=0,
                                    const std::vector<double>& tolerance_vec = {}) const;

        bool inside_dead_region(const Point& p, const int apa, const int face) const;

        bool check_dead_volume(const Point& p, const Vector& dir,
                               double step = 1.0*units::cm, double offset_x=0) const;

        bool check_signal_processing(const Point& p, const Vector& dir,
                                     double step = 1.0*units::cm, double offset_x=0) const;

        bool check_other_tracks(Facade::Cluster& main_cluster, double offset_x) const;



    private:

        // Holds static user data.  Can be reset which will also clear internal
        // data.
        StaticData m_sd;

        // Bundle the "internal" data derived from static+dynamic
        //
        // This holds whatever data that is derived from "dynamic data" (eg, the
        // groupings) that are needed to satisfy the "query methods".  This
        // should be fully derived data and in particular not include the
        // objects from the "dynamic data" themselves as their life times may be
        // shorter than the lifetime of a FiducialUtils object.
        //
        // This struct must be default-constructable, eg via InternalData{} so
        // make sure all values are initialized here.  the feed_*() methods will
        // start by clearing the InternalData via this default construction.
        struct InternalData {

            // FIXME: this is just an initial place holder to show how the
            // internal data is handled in the implementation.  Delete this
            // "dummy" attribute once some "real" internal data is added.
            size_t dummy{0};

            // save live grouping ...
            Facade::Grouping* live;
        };
        InternalData m_internal;
    };

}
#endif
