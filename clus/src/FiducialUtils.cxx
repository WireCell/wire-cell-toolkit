#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/Facade_Cluster.h"

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


    bool FiducialUtils::inside_fiducial_volume(const Point& p, 
                                               const std::vector<double>& tolerance_vec) const
    {
        return m_internal.dummy%2 && tolerance_vec.empty();       // initial bogus query as place holder
    }


    bool FiducialUtils::check_dead_volume(const Facade::Cluster& main_cluster, const Point& p, const Vector& dir, double step, const double cut_ratio, const int cut_value) const
    {
        if (!inside_fiducial_volume(p)){
            return false;
        }else{
            if (dir.magnitude()==0){
                return true;
            }else{
                Vector normalized_dir = dir;
                normalized_dir *= 1./dir.magnitude();
                Point temp_p = p;
                int num_points = 0;
                int num_points_dead = 0;
                while(inside_fiducial_volume(temp_p)){

                    num_points ++;
                    // for the temp_p, find its apa, face, and raw_temp_p ...
                    auto test_wpid = m_sd.dv->contained_by(temp_p);

                    if (test_wpid.apa() < 0 || test_wpid.face() < 0) {
                       num_points_dead ++;
                    }else{
                        const auto transform = m_sd.pcts->pc_transform(main_cluster.get_scope_transform(main_cluster.get_default_scope()));
                        double cluster_t0 = main_cluster.get_flash().time();
                        auto temp_p_raw = transform->backward(temp_p, cluster_t0, test_wpid.face(), test_wpid.apa());

                        if (inside_dead_region(temp_p_raw, test_wpid.apa(), test_wpid.face())) num_points_dead ++;
                    }
                    if (num_points - num_points_dead >=cut_value) return true;
                    
                    temp_p.x(temp_p.x() + normalized_dir.x() * step);
                    temp_p.y(temp_p.y() + normalized_dir.y() * step);
                    temp_p.z(temp_p.z() + normalized_dir.z() * step);

                }
                
                if (num_points_dead > cut_ratio * num_points){
                    return false;
                }else{
                    return true;
                }
            }
            
        }
    }
    


    bool FiducialUtils::check_signal_processing(const Facade::Cluster& main_cluster, const Point& p, const Vector& dir, double step, const double cut_ratio, const int cut_value) const
    {
        if (dir.magnitude()==0){
            return true;
        }else{
            Vector normalized_dir = dir;
            normalized_dir *= 1./dir.magnitude();
            Point temp_p = p;

            int num_points = 0;
            int num_points_dead = 0;
            
            while(inside_fiducial_volume(temp_p)){
                num_points ++;
                auto test_wpid = m_sd.dv->contained_by(temp_p);
                if (test_wpid.apa() < 0 || test_wpid.face() < 0) {
                    num_points_dead ++;
                }else{
                    // convert temp_p to raw point and find apa etc ...
                    auto transform = m_sd.pcts->pc_transform(main_cluster.get_scope_transform(main_cluster.get_default_scope()));
                    double cluster_t0 = main_cluster.get_flash().time();
                    auto temp_p_raw = transform->backward(temp_p, cluster_t0, test_wpid.face(), test_wpid.apa());

                    auto result_u = m_internal.live->get_closest_points(temp_p_raw,1.2*units::cm,test_wpid.apa(), test_wpid.face(), 0);
                    auto result_v = m_internal.live->get_closest_points(temp_p_raw,1.2*units::cm,test_wpid.apa(), test_wpid.face(), 1);
                    auto result_w = m_internal.live->get_closest_points(temp_p_raw,1.2*units::cm,test_wpid.apa(), test_wpid.face(), 2);
                    if (result_u.size() > 0 || result_v.size() > 0 || result_w.size() > 0 || inside_dead_region(temp_p_raw, test_wpid.apa(), test_wpid.face())) {
                        num_points_dead ++;
                    }
                }

                if (num_points - num_points_dead >=cut_value) return true;
            
                temp_p.x(temp_p.x() + normalized_dir.x() * step);
                temp_p.y(temp_p.y() + normalized_dir.y() * step);
                temp_p.z(temp_p.z() + normalized_dir.z() * step);
            }

            
            if (num_points_dead > cut_ratio * num_points){
                return false;
            }else{
                return true;
            }
        }
        
        return true;
    }

}
