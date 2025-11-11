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

    void Segment::set_fit_associate_vec(std::vector<WireCell::Point >& tmp_fit_pt_vec, std::vector<int>& tmp_fit_index, std::vector<bool>& tmp_fit_skip, const IDetectorVolumes::pointer& dv,const std::string& cloud_name){        
        // Store fit points in m_fits vector
        m_fits.clear();
        m_fits.reserve(tmp_fit_pt_vec.size());
        
        for (size_t i = 0; i < tmp_fit_pt_vec.size(); ++i) {
            Fit fit;
            // Convert WCP::Point to WireCell::Point
            fit.point = WireCell::Point(tmp_fit_pt_vec[i].x(), tmp_fit_pt_vec[i].y(), tmp_fit_pt_vec[i].z());
            if (i < tmp_fit_index.size()) {
                fit.index = tmp_fit_index[i];
            }
            if (i < tmp_fit_skip.size()) {
                if (tmp_fit_skip[i]) {
                    fit.flag_fix = true;
                }
            }
            m_fits.push_back(fit);
        }
        
        // Create dynamic point cloud with the fit points
        if (dv) {
            create_segment_fit_point_cloud(shared_from_this(), dv, cloud_name);
        }
        
    }

  


}
