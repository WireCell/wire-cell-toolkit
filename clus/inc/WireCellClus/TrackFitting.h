#ifndef WIRECELLCLUS_TRACKFITTING_H
#define WIRECELLCLUS_TRACKFITTING_H

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Logging.h"
#include "WireCellClus/PRGraph.h"

namespace WireCell::Clus {

    /**
     * Dedicated TrackFitting class that can be instantiated and used by 
     * other ensemble visitors without needing to be configured as a component.
     * 
     * This class encapsulates track fitting algorithms that can work on
     * individual clusters or collections of clusters.
     */
    class TrackFitting {
    public:
        
        enum class FittingType {
            Single,
            Multiple
        };

        /**
         * Constructor
         * @param fitting_type The type of fitting to perform (single or multiple tracks)
         */
        explicit TrackFitting(FittingType fitting_type = FittingType::Single);   
        virtual ~TrackFitting() = default;

        /**
         * Set the fitting type
         * @param fitting_type The new fitting type to use
         */
        void set_fitting_type(FittingType fitting_type) { m_fitting_type = fitting_type; }

        /**
         * Get the current fitting type
         * @return The current fitting type
         */
        FittingType get_fitting_type() const { return m_fitting_type; }

        void add_segment(PR::Segment* segment); 

        // after the first round of track fitting, adjust the rough path ...
        WireCell::Point adjust_rough_path(PR::Segment& segment);

        // collect charge
        void collect_charge(double dis_cut, double range_cut);


    private:
        FittingType m_fitting_type;
        
        // cluster and grouping, CTPC is from m_grouping ...
        Facade::Grouping* m_grouping{nullptr};
        std::set<Facade::Cluster*> m_clusters;

        // input segment
        std::set<PR::Segment*> m_segments;

    };

} // namespace WireCell::Clus

#endif // WIRECELLCLUS_TRACKFITTING_H