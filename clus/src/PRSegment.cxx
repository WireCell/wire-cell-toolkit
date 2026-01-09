#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRSegmentFunctions.h"

namespace WireCell::Clus::PR {

    void Segment::reset_fit_prop(){
        for (auto& fit : m_fits) {
            fit.reset();
        }
    }

    int Segment::fit_index(int i){
        if (i < 0 || static_cast<size_t>(i) >= m_fits.size()) {
            throw std::out_of_range("Invalid fit index");
        }
        return m_fits[i].index;
    }
    void Segment::fit_index(int i, int idx){
        if (i < 0 || static_cast<size_t>(i) >= m_fits.size()) {
            throw std::out_of_range("Invalid fit index");
        }
        m_fits[i].index = idx;
    }
    bool Segment::fit_flag_skip(int i){
        if (i < 0 || static_cast<size_t>(i) >= m_fits.size()) {
            throw std::out_of_range("Invalid fit index");
        }
        return m_fits[i].flag_fix;
    }
    void Segment::fit_flag_skip(int i, bool flag){
        if (i < 0 || static_cast<size_t>(i) >= m_fits.size()) {
            throw std::out_of_range("Invalid fit index");
        }
        m_fits[i].flag_fix = flag;
    }

    void Segment::set_fit_associate_vec(std::vector<PR::Fit >& tmp_fit_vec, const IDetectorVolumes::pointer& dv,const std::string& cloud_name){        
        // Store fit points in m_fits vector
        m_fits = tmp_fit_vec;

        // for (size_t i = 0; i < tmp_fit_pt_vec.size(); ++i) {
        //     Fit fit;
        //     // Convert WCP::Point to WireCell::Point
        //     fit.point = WireCell::Point(tmp_fit_pt_vec[i].x(), tmp_fit_pt_vec[i].y(), tmp_fit_pt_vec[i].z());
        //     if (i < tmp_fit_index.size()) {
        //         fit.index = tmp_fit_index[i];
        //     }
        //     if (i < tmp_fit_skip.size()) {
        //         if (tmp_fit_skip[i]) {
        //             fit.flag_fix = true;
        //         }
        //     }
        //     m_fits.push_back(fit);
        // }
        
        // Create dynamic point cloud with the fit points
        if (dv) {
            create_segment_fit_point_cloud(shared_from_this(), dv, cloud_name);
        }
        
    }

    bool Segment::reset_global_indices(){
        if (m_pc_global_indices.empty()) {
            return false;
        }
        m_pc_global_indices.clear();
        return true;
    }
    
    bool Segment::reset_global_indices(const std::string& cloud_name){
        auto it = m_pc_global_indices.find(cloud_name);
        if (it == m_pc_global_indices.end()) {
            return false;
        }
        m_pc_global_indices.erase(it);
        return true;
    }

    void Segment::clear_fit(const IDetectorVolumes::pointer& dv, const std::string& cloud_name){
        // Unset the kFit flag
        unset_flags(SegmentFlags::kFit);
        
        // Reset fit vector to match wcpts size and copy point data
        m_fits.clear();
        m_fits.resize(m_wcpts.size());
        for (size_t i = 0; i != m_wcpts.size(); i++) {
            m_fits.at(i).point = m_wcpts.at(i).point;
            // Reset other Fit fields to default values
            m_fits.at(i).dQ = -1;
            m_fits.at(i).dx = 0;
            m_fits.at(i).pu = -1;
            m_fits.at(i).pv = -1;
            m_fits.at(i).pw = -1;
            m_fits.at(i).pt = 0;
            m_fits.at(i).reduced_chi2 = -1;
            m_fits.at(i).index = -1;
            m_fits.at(i).range = -1;
            m_fits.at(i).flag_fix = false;
        }
        
        // Reset direction and particle information
        m_dirsign = 0;
        m_dir_weak = false;
        m_particle_score = 100;
        m_particle_info = nullptr;
        
        // Recreate the dynamic point cloud for fit points
        if (dv) {
            create_segment_fit_point_cloud(shared_from_this(), dv, cloud_name);
        }
        
        reset_global_indices();
    }


  


}
