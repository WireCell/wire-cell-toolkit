#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/Facade_Cluster.h"
#include <cstdlib>
#include <cstdio>

namespace WireCell::Clus::PR {

    // doc sbnd_xin/docs/pr/40: WCT_PID_WRITE_DEBUG, an env-gated object-level
    // audit of every particle_info() write that sets or clears pdg +-11
    // (same idiom as WCT_SHOWER_TOPO_DEBUG / WCT_DET_DEBUG elsewhere in this
    // file's siblings).  Catches the setter here for all ~40 call sites with
    // zero call-site churn (the __builtin_FILE/LINE default args do the
    // work); the 6 non-const-ref reseat sites and 12 direct set_pdg()
    // mutations are instrumented separately at their own call sites
    // (NeutrinoVertexFinder.cxx, NeutrinoShowerClustering.cxx) since they
    // bypass this setter entirely.
    Segment& Segment::particle_info(std::shared_ptr<Aux::ParticleInfo> pinfo,
                                     const char* _caller_file, int _caller_line) {
        static const bool dbg = std::getenv("WCT_PID_WRITE_DEBUG") != nullptr;
        if (dbg) {
            const int old_pdg = m_particle_info ? m_particle_info->pdg() : 0;
            const int new_pdg = pinfo ? pinfo->pdg() : 0;
            if (std::abs(old_pdg) == 11 || std::abs(new_pdg) == 11) {
                const int cid = cluster() ? cluster()->get_cluster_id() : -1;
                std::fprintf(stderr,
                    "PID_WRITE_DEBUG setter id=%d clus=%d gidx=%zu pdg %d -> %d  at %s:%d\n",
                    m_id, cid, m_graph_index, old_pdg, new_pdg, _caller_file, _caller_line);
            }
        }
        m_particle_info = pinfo;
        return *this;
    }

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

    void Segment::set_fit_associate_vec(std::vector<PR::Fit> tmp_fit_vec, const IDetectorVolumes::pointer& dv, const std::string& cloud_name){
        // Store fit points in m_fits vector (move to avoid a copy)
        m_fits = std::move(tmp_fit_vec);

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
            // Populate apa/face so downstream multi-APA consumers retain face info after re-fit
            if (dv) {
                auto wpid = dv->contained_by(m_wcpts.at(i).point);
                m_fits.at(i).paf = {wpid.apa(), wpid.face()};
            }
        }
        
        // Direction and particle identification are deliberately NOT reset here.
        // prototype (ProtoSegment.cxx:1012-1040): clear_fit() clears only the fit
        // arrays (fit_pt_vec, dQ/dx, pu/pv/pw/pt, reduced_chi2, fit_index_vec,
        // fit_flag_skip) and rebuilds pcloud_fit; particle_type, flag_dir,
        // dir_weak and particle_score survive untouched.  This function has a
        // single caller in each tree -- MyFCN::UpdateInfo (MyFCN.cxx:482,
        // prototype NeutrinoID_improve_vertex.h:943) -- so resetting PID here
        // discarded the identification of every leg of every fitted vertex.
        // That is observable at NeutrinoVertexFinder.cxx:2327, whose
        // `if (!sg->particle_info())` gate (prototype: `get_particle_type()==0`)
        // was therefore always true and re-ran segment_is_shower_topology on
        // legs the prototype leaves alone.  See sbnd_xin/docs/pr/28 sec. 3.3.

        // Recreate the dynamic point cloud for fit points
        if (dv) {
            create_segment_fit_point_cloud(shared_from_this(), dv, cloud_name);
        }
        
        reset_global_indices();
    }


  


}
