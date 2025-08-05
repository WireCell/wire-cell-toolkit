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
        //m_internal.dummy = dd.live.nchildren() + dd.dead.nchildren();
       
        m_internal.live = const_cast<Facade::Grouping*>(&dd.live);
    }


    // query methods

    
    bool FiducialUtils::inside_dead_region(const Point& p_raw, const int apa, const int face, const int minimal_views) const
    {
         // Convert 3D point to time and wire indices
        const auto [tind_u, wind_u] = m_internal.live->convert_3Dpoint_time_ch(p_raw, apa, face,  0);
        const auto [tind_v, wind_v] = m_internal.live->convert_3Dpoint_time_ch(p_raw, apa, face,  1);
        const auto [tind_w, wind_w] = m_internal.live->convert_3Dpoint_time_ch(p_raw, apa, face,  2);

        int dead_view_count = 0;
        
        // Check each plane (U=0, V=1, W=2)
        // Check if this wire at this time is dead
        if (m_internal.live->is_wire_dead(apa, face, 0, wind_u, tind_u)) dead_view_count++;
        if (m_internal.live->is_wire_dead(apa, face, 1, wind_v, tind_v)) dead_view_count++;
        if (m_internal.live->is_wire_dead(apa, face, 2, wind_w, tind_w)) dead_view_count++;

        // Return true if number of dead views >= minimal_views
        return dead_view_count >= minimal_views;
    }


    bool FiducialUtils::inside_fiducial_volume(const Point& p, double offset_x,
                                               const std::vector<double>& tolerance_vec) const
    {
        return m_internal.dummy%2 && tolerance_vec.empty();       // initial bogus query as place holder
    }


    bool FiducialUtils::check_dead_volume(const Point& p, const Vector& dir,
                                          double step) const
    {
        
    }


    bool FiducialUtils::check_signal_processing(const Point& p, const Vector& dir,
                                                double step) const
    {
    }


   
    
}
