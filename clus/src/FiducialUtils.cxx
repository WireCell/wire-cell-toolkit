#include "WireCellClus/FiducialUtils.h"

namespace WireCell::Clus {


    FiducialUtils::FiducialUtils(FiducialUtils::StaticData sd)
        : m_sd(sd)
    {
    }

    void FiducialUtils::feed_static(FiducialUtils::StaticData sd)
    {
        m_internal = InternalData{}; // clear any previous internal data
        m_sd = sd;
    }

    void FiducialUtils::feed_dynamic(const FiducialUtils::DynamicData& dd) {
        m_internal = InternalData{}; // clear any previous internal data

        // After the above reset, the rest of this method should be filled with
        // code to derive whatever InternalData values from the static (m_sd)
        // and dynamic data (dd) as needed to by the query methods.

        // For now, we simply set our place holder "dummy" to some meaningless value.
        m_internal.dummy = dd.live.nchildren() + dd.dead.nchildren();
        
    }


    // query methods

    
    bool FiducialUtils::inside_dead_region(const Point& p) const
    {
        return m_internal.dummy%2;       // initial bogus query as place holder
    }


    bool FiducialUtils::inside_fiducial_volume(const Point& p, double offset_x,
                                               const std::vector<double>& tolerance_vec) const
    {
        return m_internal.dummy%2 && tolerance_vec.empty();       // initial bogus query as place holder
    }


    bool FiducialUtils::check_dead_volume(const Point& p, const Vector& dir,
                                          double step, double offset_x) const
    {
        return m_internal.dummy%2;       // initial bogus query as place holder
    }


    bool FiducialUtils::check_signal_processing(const Point& p, const Vector& dir,
                                                double step, double offset_x) const
    {
        return m_internal.dummy%2;       // initial bogus query as place holder
    }


    bool FiducialUtils::check_other_tracks(Facade::Cluster& main_cluster, double offset_x) const
    {
        return m_internal.dummy%2;       // initial bogus query as place holder
    }
    
    
}
