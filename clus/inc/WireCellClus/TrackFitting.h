#ifndef WIRECELLCLUS_TRACKFITTING_H
#define WIRECELLCLUS_TRACKFITTING_H

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Logging.h"

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
         * Perform track fitting on a single cluster
         * @param cluster The cluster to perform track fitting on
         * @return true if fitting was successful, false otherwise
         */
        bool fit_track(Facade::Cluster& cluster) const;

        /**
         * Perform track fitting on multiple clusters at once
         * @param clusters Vector of clusters to perform track fitting on
         * @return number of successfully fitted tracks
         */
        size_t fit_tracks(std::vector<Facade::Cluster*>& clusters) const;


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

    private:
        FittingType m_fitting_type;
       
    };

} // namespace WireCell::Clus

#endif // WIRECELLCLUS_TRACKFITTING_H