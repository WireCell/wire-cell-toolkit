#include "ROI_refinement.h"
#include "PeakFinding.h"
#include <iostream>
#include <set>

using namespace WireCell;
using namespace WireCell::SigProc;

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif
ROI_refinement::ROI_refinement(Waveform::ChannelMaskMap &cmm, int nwire_u, int nwire_v, int nwire_w, float th_factor,
                               float fake_signal_low_th, float fake_signal_high_th, float fake_signal_low_th_ind_factor,
                               float fake_signal_high_th_ind_factor, int pad, int break_roi_loop, float th_peak,
                               float sep_peak, float low_peak_sep_threshold_pre, int max_npeaks, float sigma,
                               float th_percent, bool isWrapped)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
  , th_factor(th_factor)
  , fake_signal_low_th(fake_signal_low_th)
  , fake_signal_high_th(fake_signal_high_th)
  , fake_signal_low_th_ind_factor(fake_signal_low_th_ind_factor)
  , fake_signal_high_th_ind_factor(fake_signal_high_th_ind_factor)
  , pad(pad)
  , break_roi_loop(break_roi_loop)
  , th_peak(th_peak)
  , sep_peak(sep_peak)
  , low_peak_sep_threshold_pre(low_peak_sep_threshold_pre)
  , max_npeaks(max_npeaks)
  , sigma(sigma)
  , th_percent(th_percent)
  , log(Log::logger("sigproc"))
  , isWrapped(isWrapped)
{
    rois_u_tight.resize(nwire_u);
    rois_u_loose.resize(nwire_u);

    rois_v_tight.resize(nwire_v);
    rois_v_loose.resize(nwire_v);

    rois_w_tight.resize(nwire_w);

    for (int i = 0; i != nwire_u; i++) {
        SignalROIList temp_rois;
        rois_u_tight.at(i) = temp_rois;
    }
    for (int i = 0; i != nwire_v; i++) {
        SignalROIList temp_rois;
        rois_v_tight.at(i) = temp_rois;
    }
    for (int i = 0; i != nwire_w; i++) {
        SignalROIList temp_rois;
        rois_w_tight.at(i) = temp_rois;
    }

    for (int i = 0; i != nwire_u; i++) {
        SignalROIList temp_rois;
        rois_u_loose.at(i) = temp_rois;
    }
    for (int i = 0; i != nwire_v; i++) {
        SignalROIList temp_rois;
        rois_v_loose.at(i) = temp_rois;
    }

    for (auto it = cmm["bad"].begin(); it != cmm["bad"].end(); it++) {
        int ch = it->first;
        std::vector<std::pair<int, int>> temps;
        bad_ch_map[ch] = temps;
        for (size_t ind = 0; ind < it->second.size(); ind++) {
            bad_ch_map[ch].push_back(std::make_pair(it->second.at(ind).first, it->second.at(ind).second));
            // std::cout << ch << " " <<  << std::endl;
        }
    }
}

ROI_refinement::~ROI_refinement() { Clear(); }

void ROI_refinement::Clear()
{
    for (int i = 0; i != nwire_u; i++) {
        for (auto it = rois_u_tight.at(i).begin(); it != rois_u_tight.at(i).end(); it++) {
            delete *it;
        }
        rois_u_tight.at(i).clear();
        for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
            delete *it;
        }
        rois_u_loose.at(i).clear();
    }

    for (int i = 0; i != nwire_v; i++) {
        for (auto it = rois_v_tight.at(i).begin(); it != rois_v_tight.at(i).end(); it++) {
            delete *it;
        }
        rois_v_tight.at(i).clear();
        for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
            delete *it;
        }
        rois_v_loose.at(i).clear();
    }

    for (int i = 0; i != nwire_w; i++) {
        for (auto it = rois_w_tight.at(i).begin(); it != rois_w_tight.at(i).end(); it++) {
            delete *it;
        }
        rois_w_tight.at(i).clear();
    }

    front_rois.clear();
    back_rois.clear();
    contained_rois.clear();
}

void ROI_refinement::apply_roi(int plane, Array::array_xxf &r_data)
{
    if (plane == 0) {
        for (int irow = 0; irow != r_data.rows(); irow++) {
            // refresh ...
            Waveform::realseq_t signal(r_data.cols(), 0);
            // loop ROIs and assign data
            for (auto it = rois_u_loose.at(irow).begin(); it != rois_u_loose.at(irow).end(); it++) {
                SignalROI *roi = *it;
                // int start_bin = roi->get_start_bin();
                // int end_bin = roi->get_end_bin();

                int start_bin = roi->get_ext_start_bin();
                int end_bin = roi->get_ext_end_bin();

                float start_content = r_data(irow, start_bin);
                float end_content = r_data(irow, end_bin);
                for (int i = start_bin; i < end_bin + 1; i++) {
                    int content =
                        r_data(irow, i) -
                        ((end_content - start_content) * (i - start_bin) / (end_bin - start_bin) + start_content);
                    signal.at(i) = content;
                    //	  htemp_signal->SetBinContent(i+1,content);
                }
            }
            // load data back ..
            for (int icol = 0; icol != r_data.cols(); icol++) {
                r_data(irow, icol) = signal.at(icol);
            }
        }
    }
    else if (plane == 1) {
        for (int irow = 0; irow != r_data.rows(); irow++) {
            // refresh ...
            Waveform::realseq_t signal(r_data.cols(), 0);
            // loop ROIs and assign data
            for (auto it = rois_v_loose.at(irow).begin(); it != rois_v_loose.at(irow).end(); it++) {
                SignalROI *roi = *it;
                // int start_bin = roi->get_start_bin();
                // int end_bin = roi->get_end_bin();

                int start_bin = roi->get_ext_start_bin();
                int end_bin = roi->get_ext_end_bin();

                float start_content = r_data(irow, start_bin);
                float end_content = r_data(irow, end_bin);
                for (int i = start_bin; i < end_bin + 1; i++) {
                    int content =
                        r_data(irow, i) -
                        ((end_content - start_content) * (i - start_bin) / (end_bin - start_bin) + start_content);
                    signal.at(i) = content;
                    // htemp_signal->SetBinContent(i+1,content);
                }
            }
            // load data back ..
            for (int icol = 0; icol != r_data.cols(); icol++) {
                r_data(irow, icol) = signal.at(icol);
            }
        }
    }
    else {
        for (int irow = 0; irow != r_data.rows(); irow++) {
            // refresh ...
            Waveform::realseq_t signal(r_data.cols(), 0);

            // if (irow==69){
            // 	std::cout << "Xin: " << rois_w_tight.at(irow).size() << " " << std::endl;
            // 	for (auto it = rois_w_tight.at(irow).begin(); it!= rois_w_tight.at(irow).end(); it++){
            // 	  SignalROI *roi =  *it;
            // 	  int start_bin = roi->get_ext_start_bin();
            // 	  int end_bin = roi->get_ext_end_bin();
            // 	  std::cout << "Xin: " << start_bin << " " << end_bin << std::endl;
            // 	}
            // }

            // loop ROIs and assign data
            for (auto it = rois_w_tight.at(irow).begin(); it != rois_w_tight.at(irow).end(); it++) {
                SignalROI *roi = *it;
                // int start_bin = roi->get_start_bin();
                // int end_bin = roi->get_end_bin();

                int start_bin = roi->get_ext_start_bin();
                int end_bin = roi->get_ext_end_bin();

                //	float start_content = r_data(irow,start_bin);
                //      float end_content = r_data(irow,end_bin);
                for (int i = start_bin; i < end_bin + 1; i++) {
                    int content = r_data(
                        irow,
                        i);  // - ((end_content - start_content)*(i-start_bin)/(end_bin-start_bin) + start_content);
                    signal.at(i) = content;
                }
            }
            // load data back ..
            for (int icol = 0; icol != r_data.cols(); icol++) {
                r_data(irow, icol) = signal.at(icol);
            }
        }
    }
}

// helper function for getting ROIs
SignalROIChList &ROI_refinement::get_rois_by_plane(int plane)
{
    if (plane == 0)
        return get_u_rois();
    else if (plane == 1)
        return get_v_rois();
    else
        return get_w_rois();
}

void ROI_refinement::unlink(SignalROI *prev_roi, SignalROI *next_roi)
{
    if (front_rois.find(prev_roi) != front_rois.end()) {
        SignalROISelection &temp_rois = front_rois[prev_roi];
        auto it = find(temp_rois.begin(), temp_rois.end(), next_roi);
        if (it != temp_rois.end()) temp_rois.erase(it);
    }
    if (back_rois.find(next_roi) != back_rois.end()) {
        SignalROISelection &temp_rois = back_rois[next_roi];
        auto it = find(temp_rois.begin(), temp_rois.end(), prev_roi);
        if (it != temp_rois.end()) temp_rois.erase(it);
    }
}

void ROI_refinement::link(SignalROI *prev_roi, SignalROI *next_roi)
{
    if (front_rois.find(prev_roi) != front_rois.end()) {
        SignalROISelection &temp_rois = front_rois[prev_roi];
        auto it = find(temp_rois.begin(), temp_rois.end(), next_roi);
        if (it == temp_rois.end()) temp_rois.push_back(next_roi);
    }
    else {
        SignalROISelection temp_rois;
        temp_rois.push_back(next_roi);
        front_rois[prev_roi] = temp_rois;
    }

    if (back_rois.find(next_roi) != back_rois.end()) {
        SignalROISelection &temp_rois = back_rois[next_roi];
        auto it = find(temp_rois.begin(), temp_rois.end(), prev_roi);
        if (it == temp_rois.end()) temp_rois.push_back(prev_roi);
    }
    else {
        SignalROISelection temp_rois;
        temp_rois.push_back(prev_roi);
        back_rois[next_roi] = temp_rois;
    }
}

void ROI_refinement::load_data(int plane, const Array::array_xxf &r_data, ROI_formation &roi_form)
{
    // fill RMS
    std::vector<float> plane_rms;
    int offset = 0;
    if (plane == 0) {
        offset = 0;
        plane_rms = roi_form.get_uplane_rms();
    }
    else if (plane == 1) {
        offset = nwire_u;
        plane_rms = roi_form.get_vplane_rms();
    }
    else if (plane == 2) {
        offset = nwire_u + nwire_v;
        plane_rms = roi_form.get_wplane_rms();
    }

    // load data ...
    for (int irow = 0; irow != r_data.rows(); irow++) {
        Waveform::realseq_t signal(r_data.cols());
        if (bad_ch_map.find(irow + offset) != bad_ch_map.end()) {
            for (int icol = 0; icol != r_data.cols(); icol++) {
                bool flag = true;
                for (size_t i = 0; i != bad_ch_map[irow + offset].size(); i++) {
                    if (icol >= bad_ch_map[irow + offset].at(i).first &&
                        icol <= bad_ch_map[irow + offset].at(i).second) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    signal.at(icol) = r_data(irow, icol);
                }
                else {
                    signal.at(icol) = 0;
                }
            }
        }
        else {
            for (int icol = 0; icol != r_data.cols(); icol++) {
                signal.at(icol) = r_data(irow, icol);
            }
        }

        int chid = irow + offset;
        // load tight rois
        std::vector<std::pair<int, int>> &uboone_rois = roi_form.get_self_rois(irow + offset);
        for (size_t i = 0; i != uboone_rois.size(); i++) {
            SignalROI *tight_roi =
                new SignalROI(plane, irow + offset, uboone_rois.at(i).first, uboone_rois.at(i).second, signal);
            float threshold = plane_rms.at(irow) * th_factor;
            if (tight_roi->get_above_threshold(threshold).size() == 0) {
                delete tight_roi;
                continue;
            }

            if (plane == 0) {
                rois_u_tight[chid].push_back(tight_roi);
                // [Hongzhao] connectivity map for Induction Planes' tight ROIs was NEVER used after initialization
                // CHANGE(S): Hongzhao commented out the initialization of connectivity map for Induction Planes' tight
                // ROIs
                //       	if (chid>0){
                //       	  //form connectivity map
                //       	  for (auto it = rois_u_tight[chid-1].begin();it!=rois_u_tight[chid-1].end();it++){
                //       	    SignalROI *prev_roi = *it;
                //      	    if (tight_roi->overlap(prev_roi)){
                //      	      if (front_rois.find(prev_roi) == front_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(tight_roi);
                //      		front_rois[prev_roi] = temp_rois;
                //      	      }else{
                //      		front_rois[prev_roi].push_back(tight_roi);
                //      	      }
                //      	      if (back_rois.find(tight_roi) == back_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(prev_roi);
                //      		back_rois[tight_roi] = temp_rois;
                //      	      }else{
                //      		back_rois[tight_roi].push_back(prev_roi);
                //      	      }
                //      	    }
                //       	  }
                //     	}
                // if (chid < nwire_u-1){
                //   // add the code for the next one to be completed
                //   for (auto it = rois_u_tight[chid+1].begin();it!=rois_u_tight[chid+1].end();it++){
                //     SignalROI *next_roi = *it;
                //     if (tight_roi->overlap(next_roi)){
                //       if (back_rois.find(next_roi) == back_rois.end()){
                // 	SignalROISelection temp_rois;
                // 	temp_rois.push_back(tight_roi);
                // 	back_rois[next_roi] = temp_rois;
                //       }else{
                // 	back_rois[next_roi].push_back(tight_roi);
                //       }

                //       if (front_rois.find(tight_roi) == front_rois.end()){
                // 	SignalROISelection temp_rois;
                // 	temp_rois.push_back(next_roi);
                // 	front_rois[tight_roi] = temp_rois;
                //       }else{
                // 	front_rois[tight_roi].push_back(next_roi);
                //       }
                //     }
                //   }
                // }
            }
            else if (plane == 1) {
                rois_v_tight[chid - nwire_u].push_back(tight_roi);
                //      	if (chid>nwire_u){
                //      	  //form connectivity map
                //      	  for (auto it = rois_v_tight[chid -
                //      nwire_u-1].begin();it!=rois_v_tight[chid-nwire_u-1].end();it++){ 	    SignalROI *prev_roi = *it; 	    if
                //      (tight_roi->overlap(prev_roi)){ 	      if (front_rois.find(prev_roi) == front_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(tight_roi);
                //      		front_rois[prev_roi] = temp_rois;
                //      	      }else{
                //      		front_rois[prev_roi].push_back(tight_roi);
                //      	      }
                //      	      if (back_rois.find(tight_roi) == back_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(prev_roi);
                //      		back_rois[tight_roi] = temp_rois;
                //      	      }else{
                //      		back_rois[tight_roi].push_back(prev_roi);
                //      	      }
                //      	    }
                //      	  }
                //      	}

                // if (chid<nwire_u+nwire_v-1){
                //      	  //form connectivity map
                //      	  for (auto it = rois_v_tight[chid -
                //      nwire_u+1].begin();it!=rois_v_tight[chid-nwire_u+1].end();it++){ 	    SignalROI *next_roi = *it; 	    if
                //      (tight_roi->overlap(next_roi)){ 	      if (back_rois.find(next_roi) == back_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(tight_roi);
                //      		back_rois[next_roi] = temp_rois;
                //      	      }else{
                //      		back_rois[next_roi].push_back(tight_roi);
                //      	      }
                //      	      if (front_rois.find(tight_roi) == front_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(next_roi);
                //      		front_rois[tight_roi] = temp_rois;
                //      	      }else{
                //      		front_rois[tight_roi].push_back(next_roi);
                //      	      }
                //      	    }
                //      	  }
                //      	}
            }
            else {
                rois_w_tight[chid - nwire_u - nwire_v].push_back(tight_roi);

                if (chid > nwire_u + nwire_v) {
                    // form connectivity map
                    // [Hongzhao] should avoid fake adjacency in Collection Plane for wrapped configuration
                    // CHANGE(S): Hongzhao added protection to skip fake adjacency
                    if (chid != nwire_u + nwire_v + nwire_w / 2. ||
                        !isWrapped) {  // additional judgement to skip fake adjacency
                        for (auto it = rois_w_tight[chid - nwire_u - nwire_v - 1].begin();
                             it != rois_w_tight[chid - nwire_u - nwire_v - 1].end(); it++) {
                            SignalROI *prev_roi = *it;
                            if (tight_roi->overlap(prev_roi)) {
                                if (front_rois.find(prev_roi) == front_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(tight_roi);
                                    front_rois[prev_roi] = temp_rois;
                                }
                                else {
                                    front_rois[prev_roi].push_back(tight_roi);
                                }
                                if (back_rois.find(tight_roi) == back_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(prev_roi);
                                    back_rois[tight_roi] = temp_rois;
                                }
                                else {
                                    back_rois[tight_roi].push_back(prev_roi);
                                }
                            }
                        }
                    }
                }

                // [Hongzhao] description_tag0
                // In principle, we don't need to explicitly create the connectivity map for channel N & N+1 WHEN we
                // already LOOP OVER all ROIs in channel N to create the connectivity map with channel N-1 (think about
                // the case if N became N+1) BUT, if we only update one ROI, we need to consider both sides (see
                // generate_merge_ROIs()) In practical, channel N+1 here is empty. CHANGE(S): Hongzhao commented out the
                // creation of the connectivity map for channel N & N+1 if (chid<nwire_u+nwire_v+nwire_w-1){
                //      	  //form connectivity map
                //      	  for (auto it =
                //      rois_w_tight[chid-nwire_u-nwire_v+1].begin();it!=rois_w_tight[chid-nwire_u-nwire_v+1].end();it++){
                //      	    SignalROI *next_roi = *it;
                //      	    if (tight_roi->overlap(next_roi)){
                //      	      if (back_rois.find(next_roi) == back_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(tight_roi);
                //      		back_rois[next_roi] = temp_rois;
                //      	      }else{
                //      		back_rois[next_roi].push_back(tight_roi);
                //      	      }
                //      	      if (front_rois.find(tight_roi) == front_rois.end()){
                //      		SignalROISelection temp_rois;
                //      		temp_rois.push_back(next_roi);
                //      		front_rois[tight_roi] = temp_rois;
                //      	      }else{
                //      		front_rois[tight_roi].push_back(next_roi);
                //      	      }
                //      	    }
                //      	  }
                //      	}
            }

        }  // loop over tight rois ...

        if (plane != 2) {
            uboone_rois = roi_form.get_loose_rois(chid);
            for (size_t i = 0; i != uboone_rois.size(); i++) {
                SignalROI *loose_roi =
                    new SignalROI(plane, chid, uboone_rois.at(i).first, uboone_rois.at(i).second, signal);
                float threshold = plane_rms.at(irow) * th_factor;
                if (loose_roi->get_above_threshold(threshold).size() == 0) {
                    delete loose_roi;
                    continue;
                }
                if (plane == 0) {
                    rois_u_loose[chid].push_back(loose_roi);

                    if (chid > 0) {
                        // form connectivity map
                        for (auto it = rois_u_loose[chid - 1].begin(); it != rois_u_loose[chid - 1].end(); it++) {
                            SignalROI *prev_roi = *it;
                            if (loose_roi->overlap(prev_roi)) {
                                if (front_rois.find(prev_roi) == front_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(loose_roi);
                                    front_rois[prev_roi] = temp_rois;
                                }
                                else {
                                    front_rois[prev_roi].push_back(loose_roi);
                                }
                                if (back_rois.find(loose_roi) == back_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(prev_roi);
                                    back_rois[loose_roi] = temp_rois;
                                }
                                else {
                                    back_rois[loose_roi].push_back(prev_roi);
                                }
                            }
                        }
                    }

                    // see [hyu1] description_tag0
                    // CHANGE(S): Hongzhao commented out the creation of the connectivity map for channel N & N+1
                    //  if (chid<nwire_u-1){
                    //    //form connectivity map
                    //    for (auto it=rois_u_loose[chid+1].begin();it!=rois_u_loose[chid+1].end();it++){
                    //      SignalROI *next_roi = *it;
                    //      if (loose_roi->overlap(next_roi)){
                    // if (back_rois.find(next_roi) == back_rois.end()){
                    //   SignalROISelection temp_rois;
                    //   temp_rois.push_back(loose_roi);
                    //   back_rois[next_roi] = temp_rois;
                    // }else{
                    //   back_rois[next_roi].push_back(loose_roi);
                    // }
                    // if (front_rois.find(loose_roi) == front_rois.end()){
                    //   SignalROISelection temp_rois;
                    //   temp_rois.push_back(next_roi);
                    //   front_rois[loose_roi] = temp_rois;
                    // }else{
                    //   front_rois[loose_roi].push_back(next_roi);
                    // }
                    //      }
                    //    }
                    //  }

                    // CHANGE(S): Hongzhao added the code to include wrapped adjacency
                    // CHANGE(S): connect the LAST & FIRST osp_wire_number in an Induction plane
                    if (chid == nwire_u - 1 && isWrapped) {
                        // form connectivity map
                        for (auto it = rois_u_loose[0].begin(); it != rois_u_loose[0].end(); it++) {
                            SignalROI *next_roi = *it;
                            if (loose_roi->overlap(next_roi)) {
                                if (back_rois.find(next_roi) == back_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(loose_roi);
                                    back_rois[next_roi] = temp_rois;
                                }
                                else {
                                    back_rois[next_roi].push_back(loose_roi);
                                }
                                if (front_rois.find(loose_roi) == front_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(next_roi);
                                    front_rois[loose_roi] = temp_rois;
                                }
                                else {
                                    front_rois[loose_roi].push_back(next_roi);
                                }
                            }
                        }
                    }

                    // form contained map ...
                    for (auto it = rois_u_tight[chid].begin(); it != rois_u_tight[chid].end(); it++) {
                        SignalROI *tight_roi = *it;
                        if (tight_roi->overlap(loose_roi)) {
                            if (contained_rois.find(loose_roi) == contained_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(tight_roi);
                                contained_rois[loose_roi] = temp_rois;
                            }
                            else {
                                contained_rois[loose_roi].push_back(tight_roi);
                            }
                            // if (tight_roi->get_start_bin() < loose_roi->get_start_bin())
                            // 	std::cout << tight_roi->get_start_bin() << " " << loose_roi->get_start_bin() << " " <<
                            // tight_roi->get_chid() << " " << tight_roi->get_plane() << std::endl;
                        }
                    }
                }
                else if (plane == 1) {
                    rois_v_loose[chid - nwire_u].push_back(loose_roi);

                    if (chid > nwire_u) {
                        // form connectivity map
                        for (auto it = rois_v_loose[chid - nwire_u - 1].begin();
                             it != rois_v_loose[chid - nwire_u - 1].end(); it++) {
                            SignalROI *prev_roi = *it;
                            if (loose_roi->overlap(prev_roi)) {
                                if (front_rois.find(prev_roi) == front_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(loose_roi);
                                    front_rois[prev_roi] = temp_rois;
                                }
                                else {
                                    front_rois[prev_roi].push_back(loose_roi);
                                }
                                if (back_rois.find(loose_roi) == back_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(prev_roi);
                                    back_rois[loose_roi] = temp_rois;
                                }
                                else {
                                    back_rois[loose_roi].push_back(prev_roi);
                                }
                            }
                        }
                    }

                    // see [hyu1] description_tag0
                    // CHANGE(S): Hongzhao commented out the creation of the connectivity map for channel N & N+1
                    //  if (chid<nwire_u+nwire_v-1){
                    //    //form connectivity map
                    //    for (auto it =
                    //    rois_v_loose[chid-nwire_u+1].begin();it!=rois_v_loose[chid-nwire_u+1].end();it++){
                    //      SignalROI *next_roi = *it;
                    //      if (loose_roi->overlap(next_roi)){
                    // if (back_rois.find(next_roi) == back_rois.end()){
                    //   SignalROISelection temp_rois;
                    //   temp_rois.push_back(loose_roi);
                    //   back_rois[next_roi] = temp_rois;
                    // }else{
                    //   back_rois[next_roi].push_back(loose_roi);
                    // }
                    // if (front_rois.find(loose_roi) == front_rois.end()){
                    //   SignalROISelection temp_rois;
                    //   temp_rois.push_back(next_roi);
                    //   front_rois[loose_roi] = temp_rois;
                    // }else{
                    //   front_rois[loose_roi].push_back(next_roi);
                    // }
                    //      }
                    //    }
                    //  }

                    // CHANGE(S): Hongzhao added the code to include wrapped adjacency
                    // CHANGE(S): connect the LAST & FIRST osp_wire_number in an Induction plane
                    if (chid == nwire_u + nwire_v - 1 && isWrapped) {
                        // form connectivity map
                        for (auto it = rois_v_loose[0].begin(); it != rois_v_loose[0].end(); it++) {
                            SignalROI *next_roi = *it;
                            if (loose_roi->overlap(next_roi)) {
                                if (back_rois.find(next_roi) == back_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(loose_roi);
                                    back_rois[next_roi] = temp_rois;
                                }
                                else {
                                    back_rois[next_roi].push_back(loose_roi);
                                }
                                if (front_rois.find(loose_roi) == front_rois.end()) {
                                    SignalROISelection temp_rois;
                                    temp_rois.push_back(next_roi);
                                    front_rois[loose_roi] = temp_rois;
                                }
                                else {
                                    front_rois[loose_roi].push_back(next_roi);
                                }
                            }
                        }
                    }

                    // form contained map ...
                    for (auto it = rois_v_tight[chid - nwire_u].begin(); it != rois_v_tight[chid - nwire_u].end();
                         it++) {
                        SignalROI *tight_roi = *it;
                        if (tight_roi->overlap(loose_roi)) {
                            if (contained_rois.find(loose_roi) == contained_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(tight_roi);
                                contained_rois[loose_roi] = temp_rois;
                            }
                            else {
                                contained_rois[loose_roi].push_back(tight_roi);
                            }
                        }
                    }
                }
            }
        }

    }  // loop over signal rows
}

void ROI_refinement::CleanUpROIs(int plane)
{
    // clean up ROIs
    std::map<SignalROI *, int> ROIsaved_map;

    if (plane == 0) {
        // int counter = 0;
        for (size_t i = 0; i != rois_u_loose.size(); i++) {
            // counter += rois_u_loose.at(i).size();
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (ROIsaved_map.find(roi) == ROIsaved_map.end()) {
                    if (contained_rois.find(roi) != contained_rois.end()) {
                        // contain good stuff
                        SignalROISelection temp_rois;
                        temp_rois.push_back(roi);
                        ROIsaved_map[roi] = 1;

                        while (temp_rois.size()) {
                            SignalROI *temp_roi = temp_rois.back();
                            temp_rois.pop_back();
                            // save all its neighbour into a temporary holder
                            if (front_rois.find(temp_roi) != front_rois.end()) {
                                for (auto it1 = front_rois[temp_roi].begin(); it1 != front_rois[temp_roi].end();
                                     it1++) {
                                    if (ROIsaved_map.find(*it1) == ROIsaved_map.end()) {
                                        temp_rois.push_back(*it1);
                                        ROIsaved_map[*it1] = 1;
                                    }
                                }
                            }
                            if (back_rois.find(temp_roi) != back_rois.end()) {
                                for (auto it1 = back_rois[temp_roi].begin(); it1 != back_rois[temp_roi].end(); it1++) {
                                    if (ROIsaved_map.find(*it1) == ROIsaved_map.end()) {
                                        temp_rois.push_back(*it1);
                                        ROIsaved_map[*it1] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // remove the bad ones ...
        // int counter2 = 0;
        for (size_t i = 0; i != rois_u_loose.size(); i++) {
            SignalROISelection to_be_removed;
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (ROIsaved_map.find(roi) == ROIsaved_map.end()) {
                    //	counter2 ++;
                    to_be_removed.push_back(roi);
                    // it = rois_u_loose.at(i).erase(it);
                    // check contained map
                    if (contained_rois.find(roi) != contained_rois.end()) {
                        std::cout << "Wrong! " << std::endl;
                    }
                    // check front map
                    if (front_rois.find(roi) != front_rois.end()) {
                        for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                            auto it2 = find(back_rois[*it1].begin(), back_rois[*it1].end(), roi);
                            back_rois[*it1].erase(it2);
                        }
                        front_rois.erase(roi);
                    }
                    // check back map
                    if (back_rois.find(roi) != back_rois.end()) {
                        for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                            auto it2 = find(front_rois[*it1].begin(), front_rois[*it1].end(), roi);
                            front_rois[*it1].erase(it2);
                        }
                        back_rois.erase(roi);
                    }
                }
            }

            for (auto it = to_be_removed.begin(); it != to_be_removed.end(); it++) {
                auto it1 = find(rois_u_loose.at(i).begin(), rois_u_loose.at(i).end(), *it);
                rois_u_loose.at(i).erase(it1);
                delete (*it);
            }
        }
    }
    else if (plane == 1) {
        // int counter = 0;
        for (size_t i = 0; i != rois_v_loose.size(); i++) {
            //  counter += rois_v_loose.at(i).size();
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (ROIsaved_map.find(roi) == ROIsaved_map.end()) {
                    if (contained_rois.find(roi) != contained_rois.end()) {
                        // contain good stuff
                        SignalROISelection temp_rois;
                        temp_rois.push_back(roi);
                        ROIsaved_map[roi] = 1;

                        while (temp_rois.size()) {
                            SignalROI *temp_roi = temp_rois.back();
                            temp_rois.pop_back();
                            // save all its neighbour into a temporary holder
                            if (front_rois.find(temp_roi) != front_rois.end()) {
                                for (auto it1 = front_rois[temp_roi].begin(); it1 != front_rois[temp_roi].end();
                                     it1++) {
                                    if (ROIsaved_map.find(*it1) == ROIsaved_map.end()) {
                                        temp_rois.push_back(*it1);
                                        ROIsaved_map[*it1] = 1;
                                    }
                                }
                            }
                            if (back_rois.find(temp_roi) != back_rois.end()) {
                                for (auto it1 = back_rois[temp_roi].begin(); it1 != back_rois[temp_roi].end(); it1++) {
                                    if (ROIsaved_map.find(*it1) == ROIsaved_map.end()) {
                                        temp_rois.push_back(*it1);
                                        ROIsaved_map[*it1] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // remove the bad ones ...
        // int counter2 = 0;
        for (size_t i = 0; i != rois_v_loose.size(); i++) {
            SignalROISelection to_be_removed;
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (ROIsaved_map.find(roi) == ROIsaved_map.end()) {
                    //	counter2 ++;
                    to_be_removed.push_back(roi);
                    // it = rois_v_loose.at(i).erase(it);
                    // check contained map
                    if (contained_rois.find(roi) != contained_rois.end()) {
                        std::cout << "Wrong! " << std::endl;
                    }
                    // check front map
                    if (front_rois.find(roi) != front_rois.end()) {
                        for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                            auto it2 = find(back_rois[*it1].begin(), back_rois[*it1].end(), roi);
                            back_rois[*it1].erase(it2);
                        }
                        front_rois.erase(roi);
                    }
                    // check back map
                    if (back_rois.find(roi) != back_rois.end()) {
                        for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                            auto it2 = find(front_rois[*it1].begin(), front_rois[*it1].end(), roi);
                            front_rois[*it1].erase(it2);
                        }
                        back_rois.erase(roi);
                    }
                }
            }

            for (auto it = to_be_removed.begin(); it != to_be_removed.end(); it++) {
                auto it1 = find(rois_v_loose.at(i).begin(), rois_v_loose.at(i).end(), *it);
                rois_v_loose.at(i).erase(it1);
                delete (*it);
            }
        }
    }

    // int counter1 = 0;
    // for (int i=0;i!=rois_v_loose.size();i++){
    //   counter1+=rois_v_loose.at(i).size();
    // }

    // std::cout << counter << " " << ROIsaved_map.size() << " " << counter1 << " " << counter2 << std::endl;
}

void ROI_refinement::generate_merge_ROIs(int plane)
{
    // find tight ROIs not contained by the loose ROIs
    if (plane == 0) {
        for (int i = 0; i != nwire_u; i++) {
            std::map<SignalROI *, int> covered_tight_rois;
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (contained_rois.find(roi) != contained_rois.end()) {
                    for (auto it1 = contained_rois[roi].begin(); it1 != contained_rois[roi].end(); it1++) {
                        if (covered_tight_rois.find(*it1) == covered_tight_rois.end()) {
                            covered_tight_rois[*it1] = 1;
                        }
                    }
                }
            }
            SignalROISelection saved_rois;
            for (auto it = rois_u_tight.at(i).begin(); it != rois_u_tight.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (covered_tight_rois.find(roi) == covered_tight_rois.end()) {
                    saved_rois.push_back(roi);
                }
            }
            // if (i == 1212)
            //   std::cout << saved_rois.size() << " " << saved_rois.at(0)->get_start_bin() << " " <<
            //   saved_rois.at(0)->get_end_bin() << std::endl;

            for (auto it = saved_rois.begin(); it != saved_rois.end(); it++) {
                SignalROI *roi = *it;
                // Duplicate them
                SignalROI *loose_roi = new SignalROI(roi);

                rois_u_loose.at(i).push_back(loose_roi);

                // update all the maps
                // contained
                SignalROISelection temp_rois;
                temp_rois.push_back(roi);
                contained_rois[loose_roi] = temp_rois;
                // front map loose ROI
                if (i < nwire_u - 1) {
                    for (auto it1 = rois_u_loose.at(i + 1).begin(); it1 != rois_u_loose.at(i + 1).end(); it1++) {
                        SignalROI *next_roi = *it1;

                        if (loose_roi->overlap(next_roi)) {
                            if (back_rois.find(next_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                back_rois[next_roi] = temp_rois;
                            }
                            else {
                                back_rois[next_roi].push_back(loose_roi);
                            }

                            if (front_rois.find(loose_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(next_roi);
                                front_rois[loose_roi] = temp_rois;
                            }
                            else {
                                front_rois[loose_roi].push_back(next_roi);
                            }
                        }
                    }
                }
                // CHANGE(S): Hongzhao added the code to include wrapped adjacency
                if (i == nwire_u - 1 && isWrapped) {
                    for (auto it1 = rois_u_loose.at(0).begin(); it1 != rois_u_loose.at(0).end(); it1++) {
                        SignalROI *next_roi = *it1;

                        if (loose_roi->overlap(next_roi)) {
                            if (back_rois.find(next_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                back_rois[next_roi] = temp_rois;
                            }
                            else {
                                back_rois[next_roi].push_back(loose_roi);
                            }

                            if (front_rois.find(loose_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(next_roi);
                                front_rois[loose_roi] = temp_rois;
                            }
                            else {
                                front_rois[loose_roi].push_back(next_roi);
                            }
                        }
                    }
                }
                // back map loose ROI
                if (i > 0) {
                    for (auto it1 = rois_u_loose.at(i - 1).begin(); it1 != rois_u_loose.at(i - 1).end(); it1++) {
                        SignalROI *prev_roi = *it1;
                        if (loose_roi->overlap(prev_roi)) {
                            if (front_rois.find(prev_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                front_rois[prev_roi] = temp_rois;
                            }
                            else {
                                front_rois[prev_roi].push_back(loose_roi);
                            }
                            if (back_rois.find(loose_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(prev_roi);
                                back_rois[loose_roi] = temp_rois;
                            }
                            else {
                                back_rois[loose_roi].push_back(prev_roi);
                            }
                        }
                    }
                }
                // CHANGE(S): Hongzhao added the code to include wrapped adjacency
                if (i == 0 && isWrapped) {
                    for (auto it1 = rois_u_loose.at(nwire_u - 1).begin(); it1 != rois_u_loose.at(nwire_u - 1).end();
                         it1++) {
                        SignalROI *prev_roi = *it1;
                        if (loose_roi->overlap(prev_roi)) {
                            if (front_rois.find(prev_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                front_rois[prev_roi] = temp_rois;
                            }
                            else {
                                front_rois[prev_roi].push_back(loose_roi);
                            }
                            if (back_rois.find(loose_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(prev_roi);
                                back_rois[loose_roi] = temp_rois;
                            }
                            else {
                                back_rois[loose_roi].push_back(prev_roi);
                            }
                        }
                    }
                }
            }
        }
    }
    else if (plane == 1) {
        for (int i = 0; i != nwire_v; i++) {
            std::map<SignalROI *, int> covered_tight_rois;
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (contained_rois.find(roi) != contained_rois.end()) {
                    for (auto it1 = contained_rois[roi].begin(); it1 != contained_rois[roi].end(); it1++) {
                        if (covered_tight_rois.find(*it1) == covered_tight_rois.end()) {
                            covered_tight_rois[*it1] = 1;
                        }
                    }
                }
            }
            SignalROISelection saved_rois;
            for (auto it = rois_v_tight.at(i).begin(); it != rois_v_tight.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (covered_tight_rois.find(roi) == covered_tight_rois.end()) {
                    saved_rois.push_back(roi);
                }
            }
            // if (i == 3885-2400)
            //  std::cout << saved_rois.size() << std::endl;
            //   std::cout << saved_rois.size() << " " << saved_rois.at(0)->get_start_bin() << " " <<
            //   saved_rois.at(0)->get_end_bin() << std::endl;

            for (auto it = saved_rois.begin(); it != saved_rois.end(); it++) {
                SignalROI *roi = *it;
                // Duplicate them
                SignalROI *loose_roi = new SignalROI(roi);

                rois_v_loose.at(i).push_back(loose_roi);

                // update all the maps
                // contained
                SignalROISelection temp_rois;
                temp_rois.push_back(roi);
                contained_rois[loose_roi] = temp_rois;
                // front map loose ROI
                if (i < nwire_v - 1) {
                    for (auto it1 = rois_v_loose.at(i + 1).begin(); it1 != rois_v_loose.at(i + 1).end(); it1++) {
                        SignalROI *next_roi = *it1;

                        if (loose_roi->overlap(next_roi)) {
                            if (back_rois.find(next_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                back_rois[next_roi] = temp_rois;
                            }
                            else {
                                back_rois[next_roi].push_back(loose_roi);
                            }

                            if (front_rois.find(loose_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(next_roi);
                                front_rois[loose_roi] = temp_rois;
                            }
                            else {
                                front_rois[loose_roi].push_back(next_roi);
                            }
                        }
                    }
                }
                // CHANGE(S): Hongzhao added the code to include wrapped adjacency
                if (i == nwire_v - 1 && isWrapped) {
                    for (auto it1 = rois_v_loose.at(0).begin(); it1 != rois_v_loose.at(0).end(); it1++) {
                        SignalROI *next_roi = *it1;

                        if (loose_roi->overlap(next_roi)) {
                            if (back_rois.find(next_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                back_rois[next_roi] = temp_rois;
                            }
                            else {
                                back_rois[next_roi].push_back(loose_roi);
                            }

                            if (front_rois.find(loose_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(next_roi);
                                front_rois[loose_roi] = temp_rois;
                            }
                            else {
                                front_rois[loose_roi].push_back(next_roi);
                            }
                        }
                    }
                }
                // back map loose ROI
                if (i > 0) {
                    for (auto it1 = rois_v_loose.at(i - 1).begin(); it1 != rois_v_loose.at(i - 1).end(); it1++) {
                        SignalROI *prev_roi = *it1;
                        if (loose_roi->overlap(prev_roi)) {
                            if (front_rois.find(prev_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                front_rois[prev_roi] = temp_rois;
                            }
                            else {
                                front_rois[prev_roi].push_back(loose_roi);
                            }
                            if (back_rois.find(loose_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(prev_roi);
                                back_rois[loose_roi] = temp_rois;
                            }
                            else {
                                back_rois[loose_roi].push_back(prev_roi);
                            }
                        }
                    }
                }
                // CHANGE(S): Hongzhao added the code to include wrapped adjacency
                if (i == 0 && isWrapped) {
                    for (auto it1 = rois_v_loose.at(nwire_v - 1).begin(); it1 != rois_v_loose.at(nwire_v - 1).end();
                         it1++) {
                        SignalROI *prev_roi = *it1;
                        if (loose_roi->overlap(prev_roi)) {
                            if (front_rois.find(prev_roi) == front_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(loose_roi);
                                front_rois[prev_roi] = temp_rois;
                            }
                            else {
                                front_rois[prev_roi].push_back(loose_roi);
                            }
                            if (back_rois.find(loose_roi) == back_rois.end()) {
                                SignalROISelection temp_rois;
                                temp_rois.push_back(prev_roi);
                                back_rois[loose_roi] = temp_rois;
                            }
                            else {
                                back_rois[loose_roi].push_back(prev_roi);
                            }
                        }
                    }
                }
            }
        }
    }
}

void ROI_refinement::CheckROIs(int plane, ROI_formation &roi_form)
{
    if (plane == 0) {
        std::vector<float> &rms_u = roi_form.get_uplane_rms();

        // for (int i=0;i!=rois_u_loose.size();i++){
        //   for (auto it = rois_u_loose.at(i).begin(); it!= rois_u_loose.at(i).end();it++){
        //     SignalROI *roi = *it;
        //     int chid = roi->get_chid();
        //     if (chid != i)
        // 	std::cout << roi << std::endl;

        //     if (front_rois.find(roi)!=front_rois.end()){
        //  	for (auto it1 = front_rois[roi].begin();it1!=front_rois[roi].end();it1++){
        // 	  SignalROI *roi1 = *it1;
        // 	  int chid1 = roi1->get_chid();
        // 	  if (chid1!=i+1)
        // 	    std::cout << roi1 << " " << chid << " " << chid1 << std::endl;
        // 	}
        //     }
        //   }
        // }

        // loop over u loose
        for (size_t i = 0; i != rois_u_loose.size(); i++) {
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                int chid = roi->get_chid();
                float th;
                th = th_factor * rms_u.at(chid);

                if (front_rois.find(roi) != front_rois.end()) {
                    SignalROISelection temp_rois;
                    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int chid1 = roi1->get_chid();
                        // std::cout << "F " << i << " " << rois_u_loose.size() << " " << roi1 << " " << chid << " " <<
                        // chid1 << std::endl;
                        float th1;
                        th1 = th_factor * rms_u.at(chid1);
                        if (roi->overlap(roi1, th, th1)) {
                        }
                        else {
                            temp_rois.push_back(roi1);
                            // unlink(roi,roi1);
                        }
                    }
                    for (auto it2 = temp_rois.begin(); it2 != temp_rois.end(); it2++) {
                        unlink(roi, *it2);
                    }
                }

                if (back_rois.find(roi) != back_rois.end()) {
                    SignalROISelection temp_rois;
                    for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int chid1 = roi1->get_chid();
                        // std::cout << "B " << roi1 << " " << chid << " " << chid1 << std::endl;
                        float th1;
                        th1 = th_factor * rms_u.at(chid1);
                        if (roi->overlap(roi1, th, th1)) {
                        }
                        else {
                            temp_rois.push_back(roi1);
                            // unlink(roi,roi1);
                        }
                    }
                    for (auto it2 = temp_rois.begin(); it2 != temp_rois.end(); it2++) {
                        unlink(*it2, roi);
                    }
                }
            }
        }
    }
    else if (plane == 1) {
        std::vector<float> &rms_v = roi_form.get_vplane_rms();
        // loop over v loose
        for (size_t i = 0; i != rois_v_loose.size(); i++) {
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                int chid = roi->get_chid() - nwire_u;
                float th;
                th = th_factor * rms_v.at(chid);
                if (front_rois.find(roi) != front_rois.end()) {
                    SignalROISelection temp_rois;
                    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int chid1 = roi1->get_chid() - nwire_u;
                        float th1;
                        th1 = th_factor * rms_v.at(chid1);
                        if (roi->overlap(roi1, th, th1)) {
                        }
                        else {
                            temp_rois.push_back(roi1);
                            // unlink(roi,roi1);
                        }
                    }
                    for (auto it2 = temp_rois.begin(); it2 != temp_rois.end(); it2++) {
                        unlink(roi, *it2);
                    }
                }

                if (back_rois.find(roi) != back_rois.end()) {
                    SignalROISelection temp_rois;
                    for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int chid1 = roi1->get_chid() - nwire_u;
                        float th1;
                        th1 = th_factor * rms_v.at(chid1);
                        if (roi->overlap(roi1, th, th1)) {
                        }
                        else {
                            temp_rois.push_back(roi1);
                            // unlink(roi,roi1);
                        }
                    }
                    for (auto it2 = temp_rois.begin(); it2 != temp_rois.end(); it2++) {
                        unlink(*it2, roi);
                    }
                }
            }
        }
    }
}

void ROI_refinement::CleanUpCollectionROIs()
{
    // deal with tight ROIs,
    // scan with all the tight ROIs to look for peaks above certain threshold, put in a temporary set
    float mean_threshold = fake_signal_low_th;
    float threshold = fake_signal_high_th;  // electrons, about 1/2 of MIP per tick ...
    std::set<SignalROI *> Good_ROIs;
    for (int i = 0; i != nwire_w; i++) {
        for (auto it = rois_w_tight.at(i).begin(); it != rois_w_tight.at(i).end(); it++) {
            SignalROI *roi = *it;
            // std::cout << "Xin: " << roi->get_max_height() << " " << roi->get_average_heights() << std::endl;
            if (roi->get_above_threshold(threshold).size() != 0 || roi->get_average_heights() > mean_threshold)
                Good_ROIs.insert(roi);
        }
    }
    // for a particular ROI if it is not in, or it is not connected with one in the temporary map, then remove it
    std::list<SignalROI *> Bad_ROIs;
    for (int i = 0; i != nwire_w; i++) {
        for (auto it = rois_w_tight.at(i).begin(); it != rois_w_tight.at(i).end(); it++) {
            SignalROI *roi = *it;

            if (Good_ROIs.find(roi) != Good_ROIs.end()) continue;
            if (front_rois.find(roi) != front_rois.end()) {
                SignalROISelection next_rois = front_rois[roi];
                int flag_qx = 0;
                for (size_t i = 0; i != next_rois.size(); i++) {
                    SignalROI *roi1 = next_rois.at(i);
                    if (Good_ROIs.find(roi1) != Good_ROIs.end()) {
                        flag_qx = 1;
                        break;  // continue;
                    }
                }
                if (flag_qx == 1) continue;
            }

            if (back_rois.find(roi) != back_rois.end()) {
                SignalROISelection next_rois = back_rois[roi];
                int flag_qx = 0;
                for (size_t i = 0; i != next_rois.size(); i++) {
                    SignalROI *roi1 = next_rois.at(i);
                    if (Good_ROIs.find(roi1) != Good_ROIs.end()) {
                        flag_qx = 1;
                        break;  // continue;
                    }
                }
                if (flag_qx == 1) continue;
            }

            Bad_ROIs.push_back(roi);
        }
    }

    // remove the ROI and then update the map

    for (auto it = Bad_ROIs.begin(); it != Bad_ROIs.end(); it++) {
        SignalROI *roi = *it;
        int chid = roi->get_chid() - nwire_u - nwire_v;
        // std::cout << chid << std::endl;
        if (front_rois.find(roi) != front_rois.end()) {
            SignalROISelection next_rois = front_rois[roi];
            for (size_t i = 0; i != next_rois.size(); i++) {
                // unlink the current roi
                unlink(roi, next_rois.at(i));
            }
            front_rois.erase(roi);
        }

        if (back_rois.find(roi) != back_rois.end()) {
            SignalROISelection prev_rois = back_rois[roi];
            for (size_t i = 0; i != prev_rois.size(); i++) {
                // unlink the current roi
                unlink(prev_rois.at(i), roi);
            }
            back_rois.erase(roi);
        }
        auto it1 = find(rois_w_tight.at(chid).begin(), rois_w_tight.at(chid).end(), roi);
        if (it1 != rois_w_tight.at(chid).end()) rois_w_tight.at(chid).erase(it1);

        delete roi;
    }
}

void ROI_refinement::CleanUpInductionROIs(int plane)
{
    // deal with loose ROIs
    // focus on the isolated ones first
    float mean_threshold = fake_signal_low_th;
    float threshold = fake_signal_high_th;
    mean_threshold *= fake_signal_low_th_ind_factor;
    threshold *= fake_signal_high_th_ind_factor;

    std::list<SignalROI *> Bad_ROIs;

    // Hongzhao: bad ROIs without good neighbors will be removed by next step
    //  if (plane==0){
    //    for (int i=0;i!=nwire_u;i++){
    //      for (auto it = rois_u_loose.at(i).begin();it!=rois_u_loose.at(i).end();it++){
    // SignalROI* roi = *it;
    // //std::cout << "Xin: " << roi->get_max_height() << " " << roi->get_average_heights() << std::endl;
    // if (front_rois.find(roi)==front_rois.end() && back_rois.find(roi)==back_rois.end()){
    //   if (roi->get_above_threshold(threshold).size()==0 && roi->get_average_heights() < mean_threshold)
    //     Bad_ROIs.push_back(roi);
    // }
    //      }
    //    }
    //    for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    //      SignalROI* roi = *it;
    //      int chid = roi->get_chid();
    //      auto it1 = find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(),roi);
    //      if (it1 != rois_u_loose.at(chid).end())
    // rois_u_loose.at(chid).erase(it1);

    //      if (front_rois.find(roi)!=front_rois.end()){
    // SignalROISelection next_rois = front_rois[roi];
    // for (size_t i=0;i!=next_rois.size();i++){
    //   //unlink the current roi
    //   unlink(roi,next_rois.at(i));
    // }
    // front_rois.erase(roi);
    //      }

    //      if (back_rois.find(roi)!=back_rois.end()){
    // SignalROISelection prev_rois = back_rois[roi];
    // for (size_t i=0;i!=prev_rois.size();i++){
    //   //unlink the current roi
    //   unlink(prev_rois.at(i),roi);
    // }
    // back_rois.erase(roi);
    //      }

    //      delete roi;
    //    }
    //    Bad_ROIs.clear();
    //  }else if (plane==1){
    //    for (int i=0;i!=nwire_v;i++){
    //      for (auto it = rois_v_loose.at(i).begin();it!=rois_v_loose.at(i).end();it++){
    // SignalROI* roi = *it;
    // //	std::cout << "Xin: " << roi->get_max_height() << " " << roi->get_average_heights() << std::endl;
    // if (front_rois.find(roi)==front_rois.end() && back_rois.find(roi)==back_rois.end()){
    //   if (roi->get_above_threshold(threshold).size()==0 && roi->get_average_heights() < mean_threshold)
    //     Bad_ROIs.push_back(roi);
    // }
    //      }
    //    }
    //    for (auto it = Bad_ROIs.begin(); it!=Bad_ROIs.end(); it ++){
    //      SignalROI* roi = *it;
    //      int chid = roi->get_chid()-nwire_u;
    //      auto it1 = find(rois_v_loose.at(chid).begin(), rois_v_loose.at(chid).end(),roi);
    //      if (it1 != rois_v_loose.at(chid).end())
    // rois_v_loose.at(chid).erase(it1);

    //      if (front_rois.find(roi)!=front_rois.end()){
    // SignalROISelection next_rois = front_rois[roi];
    // for (size_t i=0;i!=next_rois.size();i++){
    //   //unlink the current roi
    //   unlink(roi,next_rois.at(i));
    // }
    // front_rois.erase(roi);
    //      }

    //      if (back_rois.find(roi)!=back_rois.end()){
    // SignalROISelection prev_rois = back_rois[roi];
    // for (size_t i=0;i!=prev_rois.size();i++){
    //   //unlink the current roi
    //   unlink(prev_rois.at(i),roi);
    // }
    // back_rois.erase(roi);
    //      }

    //      delete roi;
    //    }
    //  }

    // threshold = fake_signal_low_th;
    std::set<SignalROI *> Good_ROIs;
    if (plane == 0) {
        for (int i = 0; i != nwire_u; i++) {
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                if (roi->get_above_threshold(threshold).size() != 0 || roi->get_average_heights() > mean_threshold)
                    Good_ROIs.insert(roi);
            }
        }
        Bad_ROIs.clear();
        for (int i = 0; i != nwire_u; i++) {
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                SignalROI *roi = *it;

                if (Good_ROIs.find(roi) != Good_ROIs.end()) continue;
                if (front_rois.find(roi) != front_rois.end()) {
                    SignalROISelection next_rois = front_rois[roi];
                    int flag_qx = 0;
                    for (size_t i = 0; i != next_rois.size(); i++) {
                        SignalROI *roi1 = next_rois.at(i);
                        if (Good_ROIs.find(roi1) != Good_ROIs.end()) {
                            flag_qx = 1;
                            break;  // continue;
                        }
                    }
                    if (flag_qx == 1) continue;
                }

                if (back_rois.find(roi) != back_rois.end()) {
                    SignalROISelection next_rois = back_rois[roi];
                    int flag_qx = 0;
                    for (size_t i = 0; i != next_rois.size(); i++) {
                        SignalROI *roi1 = next_rois.at(i);
                        if (Good_ROIs.find(roi1) != Good_ROIs.end()) {
                            flag_qx = 1;
                            break;  // continue;
                        }
                    }
                    if (flag_qx == 1) continue;
                }

                Bad_ROIs.push_back(roi);
            }
        }
        for (auto it = Bad_ROIs.begin(); it != Bad_ROIs.end(); it++) {
            SignalROI *roi = *it;
            int chid = roi->get_chid();
            // std::cout << chid << std::endl;
            if (front_rois.find(roi) != front_rois.end()) {
                SignalROISelection next_rois = front_rois[roi];
                for (size_t i = 0; i != next_rois.size(); i++) {
                    // unlink the current roi
                    unlink(roi, next_rois.at(i));
                }
                front_rois.erase(roi);
            }

            if (back_rois.find(roi) != back_rois.end()) {
                SignalROISelection prev_rois = back_rois[roi];
                for (size_t i = 0; i != prev_rois.size(); i++) {
                    // unlink the current roi
                    unlink(prev_rois.at(i), roi);
                }
                back_rois.erase(roi);
            }
            auto it1 = find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(), roi);
            if (it1 != rois_u_loose.at(chid).end()) rois_u_loose.at(chid).erase(it1);

            delete roi;
        }
    }
    else if (plane == 1) {
        Good_ROIs.clear();
        for (int i = 0; i != nwire_v; i++) {
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                SignalROI *roi = *it;
                // Historically, "mean_threshold" was not included in this condition.
                // It has now been added to ensure consistency between iplane == 0 and iplane == 1.
                if (roi->get_above_threshold(threshold).size() != 0 || roi->get_average_heights() > mean_threshold) Good_ROIs.insert(roi);
            }
        }
        Bad_ROIs.clear();
        for (int i = 0; i != nwire_v; i++) {
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                SignalROI *roi = *it;

                if (Good_ROIs.find(roi) != Good_ROIs.end()) continue;
                if (front_rois.find(roi) != front_rois.end()) {
                    SignalROISelection next_rois = front_rois[roi];
                    int flag_qx = 0;
                    for (size_t i = 0; i != next_rois.size(); i++) {
                        SignalROI *roi1 = next_rois.at(i);
                        if (Good_ROIs.find(roi1) != Good_ROIs.end()) {
                            flag_qx = 1;
                            break;  // continue;
                        }
                    }
                    if (flag_qx == 1) continue;
                }

                if (back_rois.find(roi) != back_rois.end()) {
                    SignalROISelection next_rois = back_rois[roi];
                    int flag_qx = 0;
                    for (size_t i = 0; i != next_rois.size(); i++) {
                        SignalROI *roi1 = next_rois.at(i);
                        if (Good_ROIs.find(roi1) != Good_ROIs.end()) {
                            flag_qx = 1;
                            break;  // continue;
                        }
                    }
                    if (flag_qx == 1) continue;
                }

                Bad_ROIs.push_back(roi);
            }
        }
        for (auto it = Bad_ROIs.begin(); it != Bad_ROIs.end(); it++) {
            SignalROI *roi = *it;
            int chid = roi->get_chid() - nwire_u;
            // std::cout << chid << std::endl;
            if (front_rois.find(roi) != front_rois.end()) {
                SignalROISelection next_rois = front_rois[roi];
                for (size_t i = 0; i != next_rois.size(); i++) {
                    // unlink the current roi
                    unlink(roi, next_rois.at(i));
                }
                front_rois.erase(roi);
            }

            if (back_rois.find(roi) != back_rois.end()) {
                SignalROISelection prev_rois = back_rois[roi];
                for (size_t i = 0; i != prev_rois.size(); i++) {
                    // unlink the current roi
                    unlink(prev_rois.at(i), roi);
                }
                back_rois.erase(roi);
            }
            auto it1 = find(rois_v_loose.at(chid).begin(), rois_v_loose.at(chid).end(), roi);
            if (it1 != rois_v_loose.at(chid).end()) rois_v_loose.at(chid).erase(it1);

            delete roi;
        }
    }
}

void ROI_refinement::ShrinkROI(SignalROI *roi, ROI_formation &roi_form)
{
    // if(proteced_rois.find({roi->get_chid(),roi->get_start_bin()})!=proteced_rois.end()) {
    //   return;
    // }
    // get tight ROI as a inner boundary
    // get the nearby ROIs with threshold as some sort of boundary
    int start_bin = roi->get_start_bin();
    int end_bin = roi->get_end_bin();
    if (start_bin < 0 || end_bin < 0) return;

    int chid = roi->get_chid();
    int plane = roi->get_plane();
    std::vector<float> &contents = roi->get_contents();

    float threshold1 = 0;
    if (plane == 0) {
        threshold1 = roi_form.get_uplane_rms().at(chid) * th_factor;
    }
    else if (plane == 1) {
        threshold1 = roi_form.get_vplane_rms().at(chid - nwire_u) * th_factor;
    }

    int channel_save = 1240;
    int print_flag = 0;

    // std::cout << "check tight ROIs " << std::endl;
    // use to save contents
    Waveform::realseq_t temp_signal(end_bin - start_bin + 1, 0);
    // TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);

    // check tight ROIs
    if (contained_rois.find(roi) != contained_rois.end()) {
        for (auto it = contained_rois[roi].begin(); it != contained_rois[roi].end(); it++) {
            SignalROI *tight_roi = *it;
            int start_bin1 = tight_roi->get_start_bin();
            int end_bin1 = tight_roi->get_end_bin();

            if (chid == channel_save && print_flag)
                std::cout << "Tight "
                             " "
                          << start_bin1 << " " << end_bin1 << std::endl;

            for (int i = start_bin1; i <= end_bin1; i++) {
                if (i - start_bin >= 0 && i - start_bin < int(temp_signal.size())) {
                    // 	  htemp->SetBinContent(i-start_bin+1,1);
                    temp_signal.at(i - start_bin) = 1;
                }
            }
        }
    }

    // std::cout << "check front ROIs " << std::endl;

    // check front ROIs
    if (front_rois.find(roi) != front_rois.end()) {
        for (auto it = front_rois[roi].begin(); it != front_rois[roi].end(); it++) {
            SignalROI *next_roi = *it;
            int start_bin1 = next_roi->get_start_bin();
            int chid1 = next_roi->get_chid();
            int plane1 = next_roi->get_plane();
            float threshold = 0;
            if (plane1 == 0) {
                threshold = roi_form.get_uplane_rms().at(chid1) * th_factor;
            }
            else if (plane1 == 1) {
                threshold = roi_form.get_vplane_rms().at(chid1 - nwire_u) * th_factor;
            }
            std::vector<std::pair<int, int>> contents_above_threshold = next_roi->get_above_threshold(threshold);
            for (size_t i = 0; i != contents_above_threshold.size(); i++) {
                if (chid == channel_save && print_flag)
                    std::cout << "Front " << chid1 << " " << start_bin1 + contents_above_threshold.at(i).first << " "
                              << start_bin1 + contents_above_threshold.at(i).second << std::endl;

                for (int j = contents_above_threshold.at(i).first; j <= contents_above_threshold.at(i).second; j++) {
                    if (j + start_bin1 - start_bin >= 0 && j + start_bin1 - start_bin < int(temp_signal.size())) {
                        if (contents.at(j + start_bin1 - start_bin) > threshold1)
                            temp_signal.at(j + start_bin1 - start_bin) = 1;
                        //	      htemp->SetBinContent(j+start_bin1-start_bin+1,1);
                    }
                }
            }
        }
    }

    // std::cout << "check back ROIs " << std::endl;

    // check back ROIs
    if (back_rois.find(roi) != back_rois.end()) {
        for (auto it = back_rois[roi].begin(); it != back_rois[roi].end(); it++) {
            SignalROI *prev_roi = *it;
            int start_bin1 = prev_roi->get_start_bin();
            int chid1 = prev_roi->get_chid();
            int plane1 = prev_roi->get_plane();
            float threshold = 0;
            if (plane1 == 0) {
                threshold = roi_form.get_uplane_rms().at(chid1) * th_factor;
            }
            else if (plane1 == 1) {
                threshold = roi_form.get_vplane_rms().at(chid1 - nwire_u) * th_factor;
            }
            std::vector<std::pair<int, int>> contents_above_threshold = prev_roi->get_above_threshold(threshold);
            for (size_t i = 0; i != contents_above_threshold.size(); i++) {
                if (chid == channel_save && print_flag)
                    std::cout << "Back " << chid1 << " " << start_bin1 + contents_above_threshold.at(i).first << " "
                              << start_bin1 + contents_above_threshold.at(i).second << std::endl;

                for (int j = contents_above_threshold.at(i).first; j <= contents_above_threshold.at(i).second; j++) {
                    if (j + start_bin1 - start_bin >= 0 && j + start_bin1 - start_bin < int(temp_signal.size())) {
                        if (contents.at(j + start_bin1 - start_bin) > threshold1)
                            temp_signal.at(j + start_bin1 - start_bin) = 1;
                        //	      htemp->SetBinContent(j+start_bin1-start_bin+1,1);
                    }
                }
            }
        }
    }

    // std::cout << "check contents " << std::endl;

    // // consider the 1/2 of the peak as threshold;
    // float max = 0;
    // for (int i=0;i!=contents.size();i++){
    //   if (contents.at(i) > max)
    //     max = contents.at(i);
    // }
    // for (int i=0;i!=contents.size();i++){
    //   // if (contents.at(i) > max/2. && contents.at(i) > threshold1*2 ) htemp->SetBinContent(i+1,1);
    // }

    // get the first bin, and last bin, add pad
    int new_start_bin = start_bin;
    int new_end_bin = end_bin;
    for (size_t i = 0; i != temp_signal.size(); i++) {
        if (temp_signal.at(i) > 0) {
            new_start_bin = i + start_bin;
            break;
        }
    }
    for (int i = int(temp_signal.size()) - 1; i >= 0; i--) {
        if (temp_signal.at(i) > 0) {
            new_end_bin = i + start_bin;
            break;
        }
    }
    new_start_bin -= pad;
    new_end_bin += pad;
    if (new_start_bin < start_bin) new_start_bin = start_bin;
    if (new_end_bin > end_bin) new_end_bin = end_bin;

    if (chid == channel_save && print_flag)
        std::cout << "check contents "
                  << " " << start_bin << " " << end_bin << " " << new_start_bin << " " << new_end_bin << std::endl;

    // create a new ROI
    Waveform::realseq_t signal(end_bin + 1);
    // TH1F *h1 = new TH1F("h1","h1",end_bin+1,0,end_bin+1);
    for (int i = new_start_bin; i <= new_end_bin; i++) {
        signal.at(i) = contents.at(i - start_bin);
        //   h1->SetBinContent(i+1,contents.at(i-start_bin));
    }

    SignalROISelection new_rois;
    if (new_start_bin >= 0 && new_end_bin > new_start_bin) {
        SignalROI *new_roi = new SignalROI(plane, chid, new_start_bin, new_end_bin, signal);
        new_rois.push_back(new_roi);
    }

    // std::cout << "update maps " << std::endl;

    // update the list
    if (chid < nwire_u) {
        auto it = std::find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(), roi);
        rois_u_loose.at(chid).erase(it);
        for (size_t i = 0; i != new_rois.size(); i++) {
            rois_u_loose.at(chid).push_back(new_rois.at(i));
        }
    }
    else if (chid < nwire_u + nwire_v) {
        auto it = std::find(rois_v_loose.at(chid - nwire_u).begin(), rois_v_loose.at(chid - nwire_u).end(), roi);
        rois_v_loose.at(chid - nwire_u).erase(it);
        for (size_t i = 0; i != new_rois.size(); i++) {
            rois_v_loose.at(chid - nwire_u).push_back(new_rois.at(i));
        }
    }

    // update all the maps
    // update front map
    if (front_rois.find(roi) != front_rois.end()) {
        SignalROISelection next_rois = front_rois[roi];
        for (size_t i = 0; i != next_rois.size(); i++) {
            // unlink the current roi
            unlink(roi, next_rois.at(i));
            // loop new rois and link them
            for (size_t j = 0; j != new_rois.size(); j++) {
                if (new_rois.at(j)->overlap(next_rois.at(i))) link(new_rois.at(j), next_rois.at(i));
            }
        }
        front_rois.erase(roi);
    }
    // update back map
    if (back_rois.find(roi) != back_rois.end()) {
        SignalROISelection prev_rois = back_rois[roi];
        for (size_t i = 0; i != prev_rois.size(); i++) {
            // unlink the current roi
            unlink(prev_rois.at(i), roi);
            // loop new rois and link them
            for (size_t j = 0; j != new_rois.size(); j++) {
                if (new_rois.at(j)->overlap(prev_rois.at(i))) link(prev_rois.at(i), new_rois.at(j));
            }
        }
        back_rois.erase(roi);
    }

    // update contained map
    if (contained_rois.find(roi) != contained_rois.end()) {
        SignalROISelection tight_rois = contained_rois[roi];
        for (size_t i = 0; i != tight_rois.size(); i++) {
            for (size_t j = 0; j != new_rois.size(); j++) {
                if (new_rois.at(j)->overlap(tight_rois.at(i))) {
                    if (contained_rois.find(new_rois.at(j)) == contained_rois.end()) {
                        SignalROISelection temp_rois;
                        temp_rois.push_back(tight_rois.at(i));
                        contained_rois[new_rois.at(j)] = temp_rois;
                    }
                    else {
                        contained_rois[new_rois.at(j)].push_back(tight_rois.at(i));
                    }
                }
            }
        }
        contained_rois.erase(roi);
    }

    // delete the old ROI
    delete roi;

    // delete htemp;
    // delete h1;
}

void ROI_refinement::ShrinkROIs(int plane, ROI_formation &roi_form)
{
    // collect all ROIs
    SignalROISelection all_rois;
    if (plane == 0) {
        for (size_t i = 0; i != rois_u_loose.size(); i++) {
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                all_rois.push_back(*it);
            }
        }
    }
    else if (plane == 1) {
        for (size_t i = 0; i != rois_v_loose.size(); i++) {
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                all_rois.push_back(*it);
            }
        }
    }
    for (size_t i = 0; i != all_rois.size(); i++) {
        ShrinkROI(all_rois.at(i), roi_form);
    }
}

void ROI_refinement::BreakROI(SignalROI *roi, float rms)
{
    auto protected_zones = proteced_rois.equal_range({roi->get_chid(), roi->get_start_bin()});
    // if (protected_zones.first!=protected_zones.second) {
    //   return;
    // }

    // main algorithm
    int start_bin = roi->get_start_bin();
    int end_bin = roi->get_end_bin();

    if (start_bin < 0 || end_bin < 0) return;

    // if (roi->get_chid()==1240){
    //   std::cout << "xin: " << start_bin << " " << end_bin << std::endl;
    // }

    Waveform::realseq_t temp_signal(end_bin - start_bin + 1, 0);
    // TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
    std::vector<float> &contents = roi->get_contents();
    for (size_t i = 0; i != temp_signal.size(); i++) {
        temp_signal.at(i) = contents.at(i);
        //    htemp->SetBinContent(i+1,contents.at(i));
    }

    // float th_peak = 3.0;
    // float sep_peak = 6.0;
    float low_peak_sep_threshold = low_peak_sep_threshold_pre;  // electrons
    if (low_peak_sep_threshold < sep_peak * rms) low_peak_sep_threshold = sep_peak * rms;
    std::set<int> saved_boundaries;

    PeakFinding s(max_npeaks, sigma, th_percent);
    const int nfound = s.find_peak(temp_signal);
    // TSpectrum *s = new TSpectrum(200);
    // Int_t nfound = s->Search(htemp,2,"nobackground new",0.1);

    if (nfound == max_npeaks) {
        log->debug("ROI_refinement: local ch index {} (plane {}), found max peaks {} with threshold={}",
                   roi->get_chid(), roi->get_plane(), nfound, th_percent);
    }

    if (nfound > 1) {
        int npeaks = s.GetNPeaks();
        double *peak_pos = s.GetPositionX();
        double *peak_height = s.GetPositionY();
        // const int temp_length = max_npeaks + 5;
        std::vector<int> order_peak_pos;
        int npeaks_threshold = 0;
        for (int j = 0; j != npeaks; j++) {
            order_peak_pos.push_back(*(peak_pos + j) + start_bin);
            if (*(peak_height + j) > th_peak * rms) {
                npeaks_threshold++;
            }
        }

        if (npeaks_threshold > 1) {
            std::sort(order_peak_pos.begin(), order_peak_pos.end());
            float valley_pos[205];
            valley_pos[0] = start_bin;

            // find the first real valley
            float min = 1e9;
            for (int k = 0; k < order_peak_pos[0] - start_bin; k++) {
                if (temp_signal.at(k) < min) {
                    min = temp_signal.at(k);
                    valley_pos[0] = k + start_bin;
                }
            }
            if (valley_pos[0] != start_bin) {
                // temporarily push n a value ...
                order_peak_pos.push_back(0);
                //
                for (int j = npeaks - 1; j >= 0; j--) {
                    order_peak_pos[j + 1] = order_peak_pos[j];
                }
                npeaks++;
                order_peak_pos[0] = start_bin;
                for (int j = start_bin; j != valley_pos[0]; j++) {
                    if (temp_signal.at(j - start_bin) > temp_signal.at(order_peak_pos[0] - start_bin))
                        //	      htemp->GetBinContent(j-start_bin+1) >
                        //htemp->GetBinContent(order_peak_pos[0]-start_bin+1))
                        order_peak_pos[0] = j;
                }
                valley_pos[0] = start_bin;
            }

            //   std::cout << "kaka1 " << npeaks << std::endl;

            for (int j = 0; j != npeaks - 1; j++) {
                float min = 1e9;
                valley_pos[j + 1] = order_peak_pos[j];

                //	std::cout << order_peak_pos[j]-start_bin << " " << order_peak_pos[j+1]-start_bin << std::endl;

                for (int k = order_peak_pos[j] - start_bin; k < order_peak_pos[j + 1] - start_bin; k++) {
                    if (temp_signal.at(k) < min) {
                        min = temp_signal.at(k);
                        valley_pos[j + 1] = k + start_bin;
                    }
                }
            }

            // std::cout << "kaka1 " << npeaks << std::endl;

            // find the end ...
            valley_pos[npeaks] = end_bin;
            min = 1e9;
            for (int k = order_peak_pos[npeaks - 1] - start_bin; k <= end_bin - start_bin; k++) {
                if (temp_signal.at(k) < min) {
                    min = temp_signal.at(k);
                    valley_pos[npeaks] = k + start_bin;
                }
            }

            //  std::cout << "kaka1 " << npeaks << std::endl;

            if (valley_pos[npeaks] != end_bin) {
                npeaks++;
                valley_pos[npeaks] = end_bin;
                // temporarily add an elements ...
                order_peak_pos.push_back(end_bin);
                //	order_peak_pos[npeaks-1] = end_bin;
                for (int j = valley_pos[npeaks - 1]; j != valley_pos[npeaks]; j++) {
                    if (temp_signal.at(j - start_bin) > temp_signal.at(order_peak_pos[npeaks - 1] - start_bin))
                        order_peak_pos[npeaks - 1] = j;
                }
            }

            // std::cout << "kaka1 " << npeaks << std::endl;

            // if (roi->get_chid() == 1195 && roi->get_plane() == WirePlaneType_t(0)){
            // 	for (int j=0;j!=npeaks;j++){
            // 	  std::cout << valley_pos[j] << " " << htemp->GetBinContent(valley_pos[j]-start_bin+1)<< " " <<
            // order_peak_pos[j] << " " << htemp->GetBinContent( order_peak_pos[j]-start_bin+1) << " " <<
            // valley_pos[j+1] << " " << htemp->GetBinContent(valley_pos[j+1] - start_bin+1)<< " " << rms * sep_peak <<
            // std::endl ;
            // 	}
            // }

            // need to organize the peaks and valleys ...
            float valley_pos1[205];
            float peak_pos1[205];
            int npeaks1 = 0;
            // fill in the first valley;
            valley_pos1[0] = valley_pos[0];
            for (int j = 0; j < npeaks; j++) {
                if (npeaks1 > 0) {
                    // std::cout << valley_pos[j]-start_bin << " " << valley_pos1[npeaks1]-start_bin << std::endl;
                    // find the lowest valley except the first peak, except the first one
                    if (temp_signal.at(valley_pos[j] - start_bin) < temp_signal.at(valley_pos1[npeaks1] - start_bin)) {
                        valley_pos1[npeaks1] = valley_pos[j];
                    }
                }

                // if (roi->get_chid() == 1195 && roi->get_plane() == WirePlaneType_t(0)){
                //   std::cout << "c: " << order_peak_pos[j] << " " <<
                //   htemp->GetBinContent(order_peak_pos[j]-start_bin+1) << " " << valley_pos1[npeaks1] << " " <<
                //   htemp->GetBinContent(valley_pos1[npeaks1]-start_bin+1) << std::endl;
                // }

                // find the next peak
                if (temp_signal.at(order_peak_pos[j] - start_bin) - temp_signal.at(valley_pos1[npeaks1] - start_bin) >
                    low_peak_sep_threshold) {
                    peak_pos1[npeaks1] = order_peak_pos[j];
                    npeaks1++;
                    int flag1 = 0;

                    for (int k = j + 1; k != npeaks + 1; k++) {
                        // find the highest peak before ...
                        if (k <= npeaks) {
                            if (temp_signal.at(order_peak_pos[k - 1] - start_bin) >
                                temp_signal.at(peak_pos1[npeaks1 - 1] - start_bin))
                                peak_pos1[npeaks1 - 1] = order_peak_pos[k - 1];
                        }

                        if (temp_signal.at(peak_pos1[npeaks1 - 1] - start_bin) -
                                temp_signal.at(valley_pos[k] - start_bin) >
                            low_peak_sep_threshold) {
                            valley_pos1[npeaks1] = valley_pos[k];
                            j = k - 1;
                            flag1 = 1;
                            break;
                        }
                        // find the next valley
                    }
                    if (flag1 == 0) {
                        valley_pos1[npeaks1] = valley_pos[npeaks];
                        j = npeaks;
                    }

                    // if (roi->get_chid() == 1240 && roi->get_plane() == 0)
                    //   std::cout << "c: " << npeaks << " " << valley_pos1[npeaks1-1] << " " << peak_pos1[npeaks1-1] <<
                    //   " " << valley_pos1[npeaks1] << " " << rms * sep_peak << std::endl;
                }
            }
            // fill the last valley
            valley_pos1[npeaks1] = valley_pos[npeaks];
            // std::cout << roi->get_plane() << " " << roi->get_chid() << " " << npeaks << " " << npeaks1 << " " ;
            // std::cout << start_bin << " " << end_bin << std::endl;

            // sometimes NO individual peaks identified due to low_peak_sep_threshold
            // but still want to keep the entire ROI as good
            if (npeaks1 == 0) {
                npeaks1 = 1;
                valley_pos1[0] = valley_pos[0];
                valley_pos1[1] = valley_pos[npeaks];
            }

            // multi-plane protection
            std::vector<bool> section_protected;
            section_protected.resize(npeaks1, false);
            for (int j = 0; j < npeaks1; ++j) {
                auto section_sta = valley_pos1[j];
                auto section_end = valley_pos1[j + 1];
                auto section_mid = 0.5 * (section_sta + section_end);
                for (auto pzone_iter = protected_zones.first; pzone_iter != protected_zones.second; ++pzone_iter) {
                    auto pzone_sta = pzone_iter->second.first;
                    auto pzone_end = pzone_iter->second.second;
                    // if (section_sta < pzone_sta and  pzone_end < section_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    // if (section_sta < pzone_end and  pzone_end < section_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    // if (pzone_sta < section_sta and section_sta < pzone_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    // if (pzone_sta < section_end and section_end < pzone_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    if (pzone_sta < section_mid and section_mid < pzone_end) {
                        section_protected[j] = true;
                        break;
                    }
                }
            }

            // const int print_chid = 1524;

            // if (roi->get_chid() == print_chid){
            for (int j = 0; j != npeaks1; j++) {
                LogDebug("{" << roi->get_chid() << ", " << roi->get_start_bin() << "} " << valley_pos1[j] << " "
                             << peak_pos1[j] << " " << valley_pos1[j + 1] << " : " << section_protected[j]);
            }
            // }

            for (int j = 0; j != npeaks1; j++) {
                int start_pos = valley_pos1[j];
                int end_pos = valley_pos1[j + 1];

                // if (roi->get_chid()==1195)
                //   std::cout << "b " << start_pos << " " << end_pos << std::endl;

                saved_boundaries.insert(start_pos);
                saved_boundaries.insert(end_pos);
            }

            Waveform::realseq_t temp1_signal = temp_signal;
            temp_signal.clear();
            temp_signal.resize(temp1_signal.size(), 0);
            //      htemp->Reset();
            for (int j = 0; j != npeaks1; j++) {
                // int flag = 0;
                int start_pos = valley_pos1[j];
                double start_content =
                    temp1_signal.at(valley_pos1[j] - start_bin);  // htemp1->GetBinContent(valley_pos1[j]-start_bin+1);
                int end_pos = valley_pos1[j + 1];
                double end_content = temp1_signal.at(
                    valley_pos1[j + 1] - start_bin);  // htemp1->GetBinContent(valley_pos1[j+1]-start_bin+1);

                //	std::cout << "a " << start_pos << " " << end_pos << std::endl;

                if (saved_boundaries.find(start_pos) != saved_boundaries.end() ||
                    saved_boundaries.find(end_pos) != saved_boundaries.end()) {
                    if (section_protected[j]) {
                        for (int k = start_pos; k != end_pos + 1; k++) {
                            temp_signal.at(k - start_bin) = temp1_signal.at(k - start_bin);
                        }
                    }
                    else {
                        for (int k = start_pos; k != end_pos + 1; k++) {
                            double temp_content = temp1_signal.at(k - start_bin) -
                                                  (start_content + (end_content - start_content) * (k - start_pos) /
                                                                       (end_pos - start_pos));
                            temp_signal.at(k - start_bin) = temp_content;
                        }
                    }
                }
            }
            //     delete htemp1;

            // std::cout << "kaka2 " << std::endl;
            // loop through the tight ROIs it contains ...
            // if nothing in the range, add things back ...
            // if (roi->get_chid() == 1151){
            //  	std::cout << "Break:  "  << roi->get_chid() << " " << start_bin << " " << end_bin << " " <<
            //  nfound << std::endl; 	std::cout << npeaks1 << std::endl;
            // 	std::cout << contained_rois[roi].size() << std::endl;

            for (auto it = contained_rois[roi].begin(); it != contained_rois[roi].end(); it++) {
                SignalROI *temp_roi = *it;
                int temp_flag = 0;
                for (int i = temp_roi->get_start_bin(); i <= temp_roi->get_end_bin(); i++) {
                    // std::cout << i-roi->get_start_bin() << " " << i << " " << temp_roi->get_start_bin() << " " <<
                    // temp_roi->get_end_bin() << " " << roi->get_start_bin() << " " << roi->get_end_bin() << " " <<
                    // roi->get_chid() << " " << roi->get_plane() << std::endl;
                    if (i - int(roi->get_start_bin()) >= 0 && i - int(roi->get_start_bin()) < int(temp_signal.size()))
                        if (temp_signal.at(i - roi->get_start_bin()) != 0) {
                            temp_flag = 1;
                            break;
                        }
                }
                //	std::cout << temp_flag << std::endl;
                if (temp_flag == 0) {
                    for (int i = temp_roi->get_start_bin(); i <= temp_roi->get_end_bin(); i++) {
                        // std::cout << i-roi->get_start_bin() << " " << i-temp_roi->get_start_bin() << std::endl;
                        if (i - int(roi->get_start_bin()) >= 0 &&
                            i - int(roi->get_start_bin()) < int(temp_signal.size()))
                            temp_signal.at(i - roi->get_start_bin()) =
                                temp_roi->get_contents().at(i - temp_roi->get_start_bin());
                        //	    htemp->SetBinContent(i-roi->get_start_bin()+1, );
                    }
                }
            }
            // } // if (1151)
        }
    }

    // std::cout << "kaka3 " << std::endl;

    // if (roi->get_chid() == 1151){
    //   std::cout << "Break:  "  << roi->get_chid() << " " << start_bin << " " << end_bin << " " << nfound <<
    //   std::endl; for (int i=0;i<htemp->GetNbinsX();i++){
    //     std::cout << htemp->GetBinContent(i+1) <<std::endl;
    //   }
    // }

    for (int qx = 0; qx != 2; qx++) {
        // Now we should go through the system again and re-adjust the content
        std::vector<std::pair<int, int>> bins;
        for (int i = 0; i < int(temp_signal.size()); i++) {
            if (temp_signal.at(i) < th_factor * rms) {
                int start = i;
                int end = i;
                for (int j = i + 1; j < int(temp_signal.size()); j++) {
                    if (temp_signal.at(j) < th_factor * rms) {
                        end = j;
                    }
                    else {
                        break;
                    }
                }
                bins.push_back(std::make_pair(start, end));
                // if (roi->get_chid() == 1240)
                //   std::cout << qx << " " <<  start+start_bin << " " << end + start_bin << " " << 3*rms  << std::endl;
                i = end;
            }
        }

        // std::cout << "kaka4 " << std::endl;

        std::vector<int> saved_b;
        for (int i = 0; i != int(bins.size()); i++) {
            int start = bins.at(i).first;
            int end = bins.at(i).second;
            // find minimum or zero
            float min = 1e9;
            int bin_min = start;
            for (int j = start; j <= end; j++) {
                if (fabs(temp_signal.at(j)) < 1e-3) {
                    bin_min = j;
                    break;
                }
                else {
                    if (temp_signal.at(j) < min) {
                        min = temp_signal.at(j);
                        bin_min = j;
                    }
                }
            }

            // if (roi->get_chid() == 1308)
            // 	std::cout << bin_min+start_bin << std::endl;

            saved_b.push_back(bin_min);
        }

        {
            // multi-plane protection
            std::vector<bool> section_protected;
            section_protected.resize(saved_b.size(), false);
            for (unsigned int j = 0; j < saved_b.size() - 1; ++j) {
                auto section_sta = saved_b[j] + start_bin;
                auto section_end = saved_b[j + 1] + start_bin;
                auto section_mid = 0.5 * (section_sta + section_end);
                for (auto pzone_iter = protected_zones.first; pzone_iter != protected_zones.second; ++pzone_iter) {
                    auto pzone_sta = pzone_iter->second.first;
                    auto pzone_end = pzone_iter->second.second;
                    // if (section_sta < pzone_sta and pzone_end < section_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    // if (section_sta < pzone_end and pzone_end < section_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    // if (pzone_sta < section_sta and section_sta < pzone_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    // if (pzone_sta < section_end and section_end < pzone_end) {
                    //   section_protected[j] = true;
                    //   break;
                    // }
                    if (pzone_sta < section_mid and section_mid < pzone_end) {
                        section_protected[j] = true;
                        break;
                    }
                }
            }

            Waveform::realseq_t temp1_signal = temp_signal;
            temp_signal.clear();
            temp_signal.resize(temp1_signal.size(), 0);

            for (int j = 0; j < int(saved_b.size()) - 1; j++) {
                int start_pos = saved_b.at(j);
                float start_content = temp1_signal.at(start_pos);
                int end_pos = saved_b.at(j + 1);
                float end_content = temp1_signal.at(end_pos);

                if (section_protected[j]) {
                    for (int k = start_pos; k != end_pos + 1; k++) {
                        temp_signal.at(k) = temp1_signal.at(k);
                    }
                }
                else {
                    for (int k = start_pos; k != end_pos + 1; k++) {
                        double temp_content =
                            temp1_signal.at(k) -
                            (start_content + (end_content - start_content) * (k - start_pos) / (end_pos - start_pos));
                        temp_signal.at(k) = temp_content;
                    }
                }
            }
        }
    }

    // std::cout << "kaka6 " << std::endl;

    // get back to the original content
    for (int i = 0; i != int(temp_signal.size()); i++) {
        contents.at(i) = temp_signal.at(i);  // htemp->GetBinContent(i+1);
    }

    //  delete s;
    //  delete htemp;
}

void ROI_refinement::BreakROI1(SignalROI *roi)
{
    int start_bin = roi->get_start_bin();
    int end_bin = roi->get_end_bin();
    if (start_bin < 0 || end_bin < 0) return;

    Waveform::realseq_t temp_signal(end_bin - start_bin + 1, 0);
    // TH1F *htemp = new TH1F("htemp","htemp",end_bin-start_bin+1,start_bin,end_bin+1);
    std::vector<float> &contents = roi->get_contents();
    for (size_t i = 0; i != temp_signal.size(); i++) {
        temp_signal.at(i) = contents.at(i);
        //   htemp->SetBinContent(i+1,contents.at(i));
    }

    // now create many new ROIs
    std::vector<int> bins;
    for (size_t i = 0; i != temp_signal.size(); i++) {
        if (fabs(temp_signal.at(i)) < 1e-3) bins.push_back(i + start_bin);
    }
    int chid = roi->get_chid();
    int plane = roi->get_plane();
    SignalROISelection new_rois;

    // if (chid == 1274)
    //   std::cout << "BreakROI1: " << chid << " " << roi->get_start_bin() << " " << roi->get_end_bin() << " " <<
    //   bins.size()  << " " << htemp->GetBinContent(1) << " " << htemp->GetBinContent(end_bin-start_bin+1) <<
    //   std::endl;

    Waveform::realseq_t signal(end_bin + 1, 0);
    //  TH1F *h1 = new TH1F("h1","h1",end_bin+1,0,end_bin+1);
    for (int i = 0; i != int(bins.size()) - 1; i++) {
        int start_bin1 = bins.at(i);
        int end_bin1 = bins.at(i + 1);
        // if (chid == 1274)
        //   std::cout << start_bin1 << " " << end_bin1 << std::endl;
        signal.clear();
        signal.resize(end_bin + 1, 0);
        //    h1->Reset();
        for (int j = start_bin1; j <= end_bin1; j++) {
            signal.at(j) = temp_signal.at(j - start_bin);
            //      h1->SetBinContent(j+1,htemp->GetBinContent(j-start_bin+1));
        }
        if (start_bin1 >= 0 && end_bin1 > start_bin1) {
            SignalROI *sub_roi = new SignalROI(plane, chid, start_bin1, end_bin1, signal);
            new_rois.push_back(sub_roi);
        }
    }

    // update the list
    if (chid < nwire_u) {
        auto it = std::find(rois_u_loose.at(chid).begin(), rois_u_loose.at(chid).end(), roi);
        rois_u_loose.at(chid).erase(it);
        for (size_t i = 0; i != new_rois.size(); i++) {
            rois_u_loose.at(chid).push_back(new_rois.at(i));
        }
    }
    else if (chid < nwire_u + nwire_v) {
        auto it = std::find(rois_v_loose.at(chid - nwire_u).begin(), rois_v_loose.at(chid - nwire_u).end(), roi);
        rois_v_loose.at(chid - nwire_u).erase(it);
        for (size_t i = 0; i != new_rois.size(); i++) {
            rois_v_loose.at(chid - nwire_u).push_back(new_rois.at(i));
        }
    }

    // update all the maps
    // update front map
    if (front_rois.find(roi) != front_rois.end()) {
        SignalROISelection next_rois = front_rois[roi];
        for (size_t i = 0; i != next_rois.size(); i++) {
            // unlink the current roi
            unlink(roi, next_rois.at(i));
            // loop new rois and link them
            for (size_t j = 0; j != new_rois.size(); j++) {
                if (new_rois.at(j)->overlap(next_rois.at(i))) link(new_rois.at(j), next_rois.at(i));
            }
        }
        front_rois.erase(roi);
    }
    // update back map
    if (back_rois.find(roi) != back_rois.end()) {
        SignalROISelection prev_rois = back_rois[roi];
        for (size_t i = 0; i != prev_rois.size(); i++) {
            // unlink the current roi
            unlink(prev_rois.at(i), roi);
            // loop new rois and link them
            for (size_t j = 0; j != new_rois.size(); j++) {
                if (new_rois.at(j)->overlap(prev_rois.at(i))) link(prev_rois.at(i), new_rois.at(j));
            }
        }
        back_rois.erase(roi);
    }

    // update contained map
    if (contained_rois.find(roi) != contained_rois.end()) {
        SignalROISelection tight_rois = contained_rois[roi];
        for (size_t i = 0; i != tight_rois.size(); i++) {
            for (size_t j = 0; j != new_rois.size(); j++) {
                if (new_rois.at(j)->overlap(tight_rois.at(i))) {
                    if (contained_rois.find(new_rois.at(j)) == contained_rois.end()) {
                        SignalROISelection temp_rois;
                        temp_rois.push_back(tight_rois.at(i));
                        contained_rois[new_rois.at(j)] = temp_rois;
                    }
                    else {
                        contained_rois[new_rois.at(j)].push_back(tight_rois.at(i));
                    }
                }
            }
        }
        contained_rois.erase(roi);
    }

    // delete the old ROI
    delete roi;
    //  delete h1;
    //  delete htemp;
}

void ROI_refinement::BreakROIs(int plane, ROI_formation &roi_form)
{
    SignalROISelection all_rois;

    if (plane == 0) {
        std::vector<float> &rms_u = roi_form.get_uplane_rms();
        for (size_t i = 0; i != rois_u_loose.size(); i++) {
            for (auto it = rois_u_loose.at(i).begin(); it != rois_u_loose.at(i).end(); it++) {
                BreakROI(*it, rms_u.at(i));
                all_rois.push_back(*it);
            }
        }
    }
    else if (plane == 1) {
        std::vector<float> &rms_v = roi_form.get_vplane_rms();
        for (size_t i = 0; i != rois_v_loose.size(); i++) {
            for (auto it = rois_v_loose.at(i).begin(); it != rois_v_loose.at(i).end(); it++) {
                BreakROI(*it, rms_v.at(i));
                all_rois.push_back(*it);
            }
        }
    }

    for (size_t i = 0; i != all_rois.size(); i++) {
        // if (all_rois.at(i)->get_chid()==1151){
        //   std::cout << all_rois.at(i)->get_chid() << " " << all_rois.at(i)->get_start_bin() << " " <<
        //   all_rois.at(i)->get_end_bin() << std::endl; for (int j=0;j!=all_rois.at(i)->get_contents().size();j++){
        // 	std::cout << j << " " << all_rois.at(i)->get_contents().at(j) << std::endl;
        //   }
        // }

        BreakROI1(all_rois.at(i));
    }
}

void ROI_refinement::refine_data(int plane, ROI_formation &roi_form)
{
    // if (plane==2) std::cout << "Xin: " << rois_w_tight.at(69).size() << " " << std::endl;

    // std::cout << "Clean up loose ROIs" << std::endl;
    CleanUpROIs(plane);
    // std::cout << "Generate more loose ROIs from isolated good tight ROIs" << std::endl;
    generate_merge_ROIs(plane);

    // if (plane==2)  std::cout << "Xin: " << rois_w_tight.at(69).size() << " " << std::endl;

    for (int qx = 0; qx != break_roi_loop; qx++) {
        // std::cout << "Break loose ROIs" << std::endl;
        BreakROIs(plane, roi_form);
        // std::cout << "Clean up ROIs 2nd time" << std::endl;
        CheckROIs(plane, roi_form);
        CleanUpROIs(plane);
    }

    // if (plane==2)  std::cout << "Xin: " << rois_w_tight.at(69).size() << " " << std::endl;

    //  std::cout << "Shrink ROIs" << std::endl;
    ShrinkROIs(plane, roi_form);
    // std::cout << "Clean up ROIs 3rd time" << std::endl;
    CheckROIs(plane, roi_form);
    CleanUpROIs(plane);

    //  if (plane==2)  std::cout << "Xin: " << rois_w_tight.at(69).size() << " " << std::endl;

    // Further reduce fake hits
    // std::cout << "Remove fake hits " << std::endl;
    if (plane == 2) {
        CleanUpCollectionROIs();
    }
    else {
        CleanUpInductionROIs(plane);
    }

    ExtendROIs(plane);
    // TestROIs();
    // if (plane==2)  std::cout << "Xin: " << rois_w_tight.at(69).size() << " " << std::endl;
}

// Basically rewrite the refine_data()
// But now do the ROI operation step by step
void ROI_refinement::refine_data_debug_mode(int plane, ROI_formation &roi_form, const std::string &cmd)
{
    if (cmd == "CleanUpROIs") {
        // std::cout << "Clean up loose ROIs" << std::endl;
        CleanUpROIs(plane);
        // std::cout << "Generate more loose ROIs from isolated good tight ROIs" << std::endl;
        generate_merge_ROIs(plane);
    }

    else if (cmd == "BreakROIs") {
        // std::cout << "Break loose ROIs" << std::endl;
        BreakROIs(plane, roi_form);
        // std::cout << "Clean up ROIs 2nd time" << std::endl;
        CheckROIs(plane, roi_form);
        CleanUpROIs(plane);
    }

    else if (cmd == "ShrinkROIs") {
        //  std::cout << "Shrink ROIs" << std::endl;
        ShrinkROIs(plane, roi_form);
        // std::cout << "Clean up ROIs 3rd time" << std::endl;
        CheckROIs(plane, roi_form);
        CleanUpROIs(plane);
    }

    else if (cmd == "ExtendROIs") {
        // Further reduce fake hits
        // std::cout << "Remove fake hits " << std::endl;
        if (plane == 2) {
            CleanUpCollectionROIs();
        }
        else {
            CleanUpInductionROIs(plane);
        }

        ExtendROIs(plane);
        // TestROIs();
    }
}

namespace {
    WireCell::Waveform::real_t feature_val(const std::vector<WireCell::Waveform::real_t> &vec, int start, int end)
    {
        LogDebug("feature_val: " << vec.size() << " " << start << " " << end);
        WireCell::Waveform::real_t ret = -1e9;
        if (vec.size() == 0) return ret;
        if (start < 0) start = 0;
        if (end < 0) end = 0;
        if (start >= (int) vec.size()) start = vec.size() - 1;
        if (end >= (int) vec.size()) end = vec.size() - 1;
        // size_t mid = (start + end) / 2;
        // if (vec[start] > ret) ret = vec[start];
        // if (vec[mid] > ret) ret = vec[mid];
        // if (vec[end] > ret) ret = vec[end];
        // LogDebug("feature_val: " << vec[start] << " " << vec[mid] << " " << vec[end] << " " << ret);
        ret = *(std::max_element(vec.begin(), vec.end()));
        LogDebug("feature_val: " << ret);
        return ret;
    }
}  // namespace

// mp3: 3 plane protection based on cleanup_roi_traces
void ROI_refinement::MP3ROI(const int plane, const IAnodePlane::pointer anode, const IAnodeFace::pointer face,
                            const std::map<int, int>& map_ch, ROI_formation& roi_form,
                            const double mp_th1, const double mp_th2,
                            const int tick_resolution, const int wire_resolution,
                            const int nbounds_layers, const std::vector<int> iplane2layer)
{
    //log->info("ROI_refinement::MP3ROI:");
    LogDebug("mp_th1: " << mp_th1 << ", mp_th2: " << mp_th2);
    std::set<int> print_chids = {1441, 875};

    if (plane == 2) return;
    int ref_planes[2] = {0, 2};
    if (plane == 0) ref_planes[0] = 1;

    std::map<int, std::vector<WireCell::RayGrid::coordinate_t>> map_tick_coord[3];
    std::map<int, std::map<int, SignalROI *>> map_tick_pitch_roi[3];

    for (int iplane = 0; iplane < 3; ++iplane) {
        auto rois_chan = get_rois_by_plane(iplane);
        for (auto rois : rois_chan) {
            for (auto roi : rois) {
                auto chid = map_ch.at(roi->get_chid());
                WireCell::RayGrid::coordinate_t coord;
                auto ch = anode->channel(chid);
                auto wires = ch->wires();
                for (auto wire : wires) {
                    if (face->which() != wire->planeid().face()) continue;
                    auto pit_id = wire->index();
                    coord.grid = pit_id;
                    coord.layer = iplane2layer[iplane] + nbounds_layers;
                    // coord.layer = iplane + nbounds_layers;
                    // integer division choose floor, including end bin seems better
                    for (int tick = roi->get_start_bin() / tick_resolution; tick <= roi->get_end_bin() / tick_resolution;
                         ++tick) {
                        int content_id = tick * tick_resolution - roi->get_start_bin();
                        const auto feat_val = feature_val(roi->get_contents(), content_id, content_id+tick_resolution);

                        //  if (print_chids.find(chid)!=print_chids.end())
                        LogDebug(tick * tick_resolution << ", " << chid << " : {" << roi->get_chid() << ", "
                                                        << roi->get_start_bin() << ", " << roi->get_end_bin()
                                                        << " } : " << feat_val);

                        if (feat_val < mp_th2) continue;

                        if (map_tick_pitch_roi[iplane].find(tick) == map_tick_pitch_roi[iplane].end()) {
                            std::map<int, SignalROI *> mtmp;
                            mtmp[pit_id] = roi;
                            map_tick_pitch_roi[iplane][tick] = mtmp;
                        }
                        else {
                            map_tick_pitch_roi[iplane][tick][pit_id] = roi;
                        }

                         LogDebug(iplane << ", " << tick*tick_resolution << ", " << pit_id);

                        if (feat_val <
                            //  mp_th1 * roi_form.get_rms_by_plane(iplane).at(roi->get_chid())
                            mp_th1)
                            continue;

                         LogDebug(iplane << ", " << tick*tick_resolution << ", " << coord.grid);

                        if (map_tick_coord[iplane].find(tick) == map_tick_coord[iplane].end()) {
                            std::vector<WireCell::RayGrid::coordinate_t> vtmp;
                            vtmp.push_back(coord);
                            map_tick_coord[iplane][tick] = vtmp;
                        }
                        else {
                            map_tick_coord[iplane][tick].push_back(coord);
                        }
                    }
                }
            }
        }

        LogDebug(" map_tick_coord:     " << map_tick_coord[iplane].size()
                << " map_tick_pitch_roi: " << map_tick_pitch_roi[iplane].size());
    }

    WireCell::RayGrid::layer_index_t layer = iplane2layer[plane] + nbounds_layers;
    // WireCell::RayGrid::layer_index_t layer = plane + nbounds_layers;
    for (auto tc1 : map_tick_coord[ref_planes[0]]) {
        for (auto tc2 : map_tick_coord[ref_planes[1]]) {
            if (tc2.first != tc1.first) continue;
            auto tick = tc1.first;
            for (auto c1 : tc1.second) {
                for (auto c2 : tc2.second) {
                    auto pitch = face->raygrid().pitch_location(c1, c2, layer);
                    auto index = face->raygrid().pitch_index(pitch, layer);
                    LogDebug(tick*tick_resolution << ", " << c1.grid << ", " << c2.grid);
                    if (map_tick_pitch_roi[plane].find(tick) == map_tick_pitch_roi[plane].end()) continue;
                    for (auto pitch_roi : map_tick_pitch_roi[plane][tick]) {
                        auto pitch = pitch_roi.first;
                        if (abs(pitch - index) < wire_resolution) {
                            int sta = tick * tick_resolution;
                            int end = sta + tick_resolution;
                            std::pair<int, int> key = {pitch_roi.second->get_chid(), pitch_roi.second->get_start_bin()};

                            if (proteced_rois.find(key) == proteced_rois.end()) {
                                proteced_rois.insert({key, {sta, end}});
                                // std::cout << "{" << key.first << ", " << key.second << "} : {"
                                //           << sta << ", " << end << "}" << std::endl;
                            }
                            else {
                                auto last_iter = --proteced_rois.upper_bound(key);
                                auto last_end = last_iter->second.second;
#ifdef __DEBUG__
                                if (print_chids.find(map_ch.at(key.first)) != print_chids.end()) {
                                    LogDebug("{" << key.first << ", " << key.second << "} : {"
                                                 << last_iter->second.first << ", " << last_iter->second.second
                                                 << "} -> {" << sta << ", " << end << "}");
                                }
#endif
                                if (sta - last_end <= tick_resolution) {
                                    last_iter->second.second = end;
                                }
                                else {
                                    proteced_rois.insert({key, {sta, end}});

                                    LogDebug("split: {" << key.first << ", " << key.second << "} ");
                                }
                            }
#ifdef __DEBUG__
                            if (print_chids.find(map_ch.at(key.first)) != print_chids.end()) {
                                LogDebug(tick * tick_resolution
                                         << " { " << map_ch.at(map_tick_pitch_roi[plane][tick][pitch]->get_chid())
                                         << ", " << map_ch.at(pitch_roi.second->get_chid()) << ", "
                                         << map_ch.at(map_tick_pitch_roi[ref_planes[0]][tick][c1.grid]->get_chid())
                                         << ", "
                                         << map_ch.at(map_tick_pitch_roi[ref_planes[1]][tick][c2.grid]->get_chid())
                                         << " } ");
                            }
#endif
                        }
                    }
                }
            }
        }
    }

#ifdef __DEBUG__
    {
        LogDebug(" total rois:          " << get_rois_by_plane(plane).size() << " at plane " << plane
                                          << " protected_rois:      " << proteced_rois.size());
        std::multimap<int, std::pair<std::pair<int, int>, std::pair<int, int>>> tmp;
        for (auto roi : proteced_rois) {
            tmp.insert({map_ch.at(roi.first.first), roi});
        }
        for (auto ch : tmp) {
            LogDebug(ch.first << " : {" << ch.second.first.first << ", " << ch.second.first.second << "} : {"
                              << ch.second.second.first << ", " << ch.second.second.second << "}, ");
        }
    }
#endif
}

// mp2: 2 plane protection based on cleaup ROI
void ROI_refinement::MP2ROI(const int target_plane, const IAnodePlane::pointer anode, const IAnodeFace::pointer face,
                            const std::map<int, int>& map_roichid_anodechid, ROI_formation& roi_form,
                            const double mp_th1, const double mp_th2,
                            const int tick_resolution, const int wire_resolution,
                            const int nbounds_layers, const std::vector<int> iplane2layer)
{
    //log->info("ROI_refinement::MP2ROI:");
    LogDebug("mp_th1: " << mp_th1 << ", mp_th2: " << mp_th2);
    std::set<int> print_chids = {1441, 875};

    // reverse map
    std::map<int, int> map_anodechid_roichid;
    for (auto r2a : map_roichid_anodechid) {
        map_anodechid_roichid[r2a.second] = r2a.first;
    }

    std::map<int, int> map_wireid_roichid[3];
    // Reusing the conversion vector (default: {0,1,2}) from OmnibusSigProc for simplicity.
    // For a more flexible design, consider making this configurable.
    std::vector<int> ilayer2plane = iplane2layer;
    for (auto chident : anode->channels()) {  // Anode chiid
        auto ilayer = anode->resolve(chident).index();
        auto iplane = ilayer2plane[ilayer];
        auto roichid = map_anodechid_roichid[chident];
        auto ch = anode->channel(chident);
        auto wires = ch->wires();
        for (auto wire : wires) {
            if (face->which() != wire->planeid().face()) continue;
            auto wireid = wire->index();
            map_wireid_roichid[iplane][wireid] = roichid;
        }
    }

    if (target_plane == 2) return;
    int ref_planes[2] = {0, 2};
    if (target_plane == 0) ref_planes[0] = 1;

    std::map<int, std::vector<WireCell::RayGrid::coordinate_t>> map_tick_coord[3];
    std::map<int, std::map<int, SignalROI *>> map_tick_pitch_roi[3];

    for (int iplane = 0; iplane < 3; ++iplane) {
        auto rois_chan = get_rois_by_plane(iplane);
        for (auto rois : rois_chan) {
            for (auto roi : rois) {
                auto chid = map_roichid_anodechid.at(roi->get_chid());
                WireCell::RayGrid::coordinate_t coord;
                auto ch = anode->channel(chid);
                auto wires = ch->wires();
                for (auto wire : wires) {
                    if (face->which() != wire->planeid().face()) continue;
                    auto pit_id = wire->index();
                    coord.grid = pit_id;
                    coord.layer = iplane2layer[iplane] + nbounds_layers;
                    // coord.layer = iplane + nbounds_layers;
                    for (int tick = roi->get_start_bin() / tick_resolution; tick <= roi->get_end_bin() / tick_resolution;
                         ++tick) {
                        int content_id = tick * tick_resolution - roi->get_start_bin();
                        const auto feat_val = feature_val(roi->get_contents(), content_id, content_id+tick_resolution);

                        if (feat_val < mp_th2) continue;

                        if (map_tick_pitch_roi[iplane].find(tick) == map_tick_pitch_roi[iplane].end()) {
                            std::map<int, SignalROI *> mtmp;
                            mtmp[pit_id] = roi;
                            map_tick_pitch_roi[iplane][tick] = mtmp;
                        }
                        else {
                            map_tick_pitch_roi[iplane][tick][pit_id] = roi;
                        }

                        //  LogDebug(iplane << ", " << tick*tick_resolution << ", " << pit_id);

                        if (feat_val <
                            //  mp_th1 * roi_form.get_rms_by_plane(iplane).at(roi->get_chid())
                            mp_th1)
                            continue;

                        //  LogDebug(iplane << ", " << tick*tick_resolution << ", " << coord.grid);

                        if (map_tick_coord[iplane].find(tick) == map_tick_coord[iplane].end()) {
                            std::vector<WireCell::RayGrid::coordinate_t> vtmp;
                            vtmp.push_back(coord);
                            map_tick_coord[iplane][tick] = vtmp;
                        }
                        else {
                            map_tick_coord[iplane][tick].push_back(coord);
                        }
                    }
                }
            }
        }

        LogDebug(" map_tick_coord:     " << map_tick_coord[iplane].size()
                << " map_tick_pitch_roi: " << map_tick_pitch_roi[iplane].size());
    }

    WireCell::RayGrid::layer_index_t layer = iplane2layer[target_plane] + nbounds_layers;
    // WireCell::RayGrid::layer_index_t layer = target_plane + nbounds_layers;
    for (auto tc1 : map_tick_coord[ref_planes[0]]) {
        for (auto tc2 : map_tick_coord[ref_planes[1]]) {
            if (tc2.first != tc1.first) continue;
            auto tick = tc1.first;
            for (auto c1 : tc1.second) {
                for (auto c2 : tc2.second) {
                    auto pitch = face->raygrid().pitch_location(c1, c2, layer);
                    auto index = face->raygrid().pitch_index(pitch, layer);
                    // LogDebug(tick*tick_resolution << ", " << c1.grid << ", " << c2.grid);

                    for (auto pitch_target = index - wire_resolution + 1; pitch_target < index + wire_resolution;
                         ++pitch_target) {
                        if (map_wireid_roichid[target_plane].find(pitch_target) ==
                            map_wireid_roichid[target_plane].end()) {
                            LogDebug("no roichid found for pitch_target: " << pitch_target << " layer: " << layer
                                                                           << " pitch: " << pitch
                                                                           << " index: " << index);
                            continue;
                        }
                        std::pair<int, int> key = {map_wireid_roichid[target_plane][pitch_target], 0};
                        int sta = tick * tick_resolution;
                        int end = sta + tick_resolution;
                        mp_rois.insert({key, {sta, end}});
                    }
                }
            }
        }
    }

#ifdef __DEBUG__
    {
        LogDebug(" total rois:          " << get_rois_by_plane(target_plane).size() << " at target_plane "
                                          << target_plane << " protected_rois:      " << mp_rois.size());
        std::multimap<int, std::pair<std::pair<int, int>, std::pair<int, int>>> tmp;
        for (auto roi : mp_rois) {
            tmp.insert({map_roichid_anodechid.at(roi.first.first), roi});
        }
        for (auto ch : tmp) {
            LogDebug(ch.first << " : {" << ch.second.first.first << ", " << ch.second.first.second << "} : {"
                              << ch.second.second.first << ", " << ch.second.second.second << "}, ");
        }
    }
#endif
}

void ROI_refinement::TestROIs()
{
    for (int chid = 0; chid != nwire_u; chid++) {
        for (auto it = rois_u_loose.at(chid).begin(); it != rois_u_loose.at(chid).end(); it++) {
            SignalROI *roi = *it;
            // loop through front
            for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                SignalROI *roi1 = *it1;
                if (find(rois_u_loose.at(chid + 1).begin(), rois_u_loose.at(chid + 1).end(), roi1) ==
                    rois_u_loose.at(chid + 1).end())
                    std::cout << chid << " u " << +1 << " " << roi << " " << roi1 << std::endl;
            }

            // loop through back
            for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                SignalROI *roi1 = *it1;
                if (find(rois_u_loose.at(chid - 1).begin(), rois_u_loose.at(chid - 1).end(), roi1) ==
                    rois_u_loose.at(chid - 1).end())
                    std::cout << chid << " u " << -1 << " " << roi << " " << roi1 << std::endl;
            }
        }
    }

    for (int chid = 0; chid != nwire_v; chid++) {
        for (auto it = rois_v_loose.at(chid).begin(); it != rois_v_loose.at(chid).end(); it++) {
            SignalROI *roi = *it;
            // loop through front
            for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                SignalROI *roi1 = *it1;
                if (find(rois_v_loose.at(chid + 1).begin(), rois_v_loose.at(chid + 1).end(), roi1) ==
                    rois_v_loose.at(chid + 1).end())
                    std::cout << chid << " v " << +1 << " " << roi << " " << roi1 << std::endl;
            }

            // loop through back
            for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                SignalROI *roi1 = *it1;
                if (find(rois_v_loose.at(chid - 1).begin(), rois_v_loose.at(chid - 1).end(), roi1) ==
                    rois_v_loose.at(chid - 1).end())
                    std::cout << chid << " v " << -1 << " " << roi << " " << roi1 << std::endl;
            }
        }
    }

    for (int chid = 0; chid != nwire_w; chid++) {
        for (auto it = rois_w_tight.at(chid).begin(); it != rois_w_tight.at(chid).end(); it++) {
            SignalROI *roi = *it;
            // loop through front
            for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                SignalROI *roi1 = *it1;
                if (find(rois_w_tight.at(chid + 1).begin(), rois_w_tight.at(chid + 1).end(), roi1) ==
                    rois_w_tight.at(chid + 1).end())
                    std::cout << chid << " w " << +1 << " " << roi << " " << roi1 << std::endl;
            }

            // loop through back
            for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                SignalROI *roi1 = *it1;
                if (find(rois_w_tight.at(chid - 1).begin(), rois_w_tight.at(chid - 1).end(), roi1) ==
                    rois_w_tight.at(chid - 1).end())
                    std::cout << chid << " w " << -1 << " " << roi << " " << roi1 << std::endl;
            }
        }
    }
}

void ROI_refinement::ExtendROIs(int plane)
{
    bool flag = true;
    if (plane == 0) {
        // U plane ...
        for (int chid = 0; chid != nwire_u; chid++) {
            rois_u_loose.at(chid).sort(CompareRois());
            for (auto it = rois_u_loose.at(chid).begin(); it != rois_u_loose.at(chid).end(); it++) {
                SignalROI *roi = *it;
                // initialize the extended bins ...
                roi->set_ext_start_bin(roi->get_start_bin());
                roi->set_ext_end_bin(roi->get_end_bin());

                if (flag) {
                    // loop through front
                    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int ext_start_bin = roi->get_ext_start_bin();
                        int ext_end_bin = roi->get_ext_end_bin();
                        if (ext_start_bin > roi1->get_start_bin()) ext_start_bin = roi1->get_start_bin();
                        if (ext_end_bin < roi1->get_end_bin()) ext_end_bin = roi1->get_end_bin();
                        roi->set_ext_start_bin(ext_start_bin);
                        roi->set_ext_end_bin(ext_end_bin);
                        //	std::cout << roi->get_chid() << " " << roi1->get_chid() << " " << roi->get_start_bin()
                        //<< " " << roi->get_end_bin() << " " << roi1->get_start_bin() << " " << roi1->get_end_bin() <<
                        //std::endl;
                    }

                    // loop through back
                    for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int ext_start_bin = roi->get_ext_start_bin();
                        int ext_end_bin = roi->get_ext_end_bin();
                        if (ext_start_bin > roi1->get_start_bin()) ext_start_bin = roi1->get_start_bin();
                        if (ext_end_bin < roi1->get_end_bin()) ext_end_bin = roi1->get_end_bin();
                        roi->set_ext_start_bin(ext_start_bin);
                        roi->set_ext_end_bin(ext_end_bin);
                        //	std::cout << roi->get_chid() << " " << roi1->get_chid() << " " << roi->get_start_bin()
                        //<< " " << roi->get_end_bin() << " " << roi1->get_start_bin() << " " << roi1->get_end_bin() <<
                        //std::endl;
                    }
                    // std::cout << chid << " " << roi->get_start_bin() << " " << roi->get_end_bin() << std::endl;
                }
            }

            if (flag) {
                SignalROI *prev_roi = 0;
                for (auto it = rois_u_loose.at(chid).begin(); it != rois_u_loose.at(chid).end(); it++) {
                    SignalROI *roi = *it;
                    if (prev_roi != 0) {
                        if (prev_roi->get_ext_end_bin() > roi->get_ext_start_bin()) {
                            prev_roi->set_ext_end_bin(int(prev_roi->get_end_bin() * 0.5 + roi->get_start_bin() * 0.5));
                            roi->set_ext_start_bin(int(prev_roi->get_end_bin() * 0.5 + roi->get_start_bin() * 0.5));
                        }
                    }
                    prev_roi = roi;
                }
            }
            // for (auto it = rois_u_loose.at(chid).begin(); it!= rois_u_loose.at(chid).end();it++){
            //   SignalROI *roi =  *it;
            //   std::cout << chid << " " << roi->get_ext_start_bin() << " " << roi->get_ext_end_bin() << std::endl;
            // }
        }
    }
    else if (plane == 1) {
        // V plane
        for (int chid = 0; chid != nwire_v; chid++) {
            rois_v_loose.at(chid).sort(CompareRois());
            for (auto it = rois_v_loose.at(chid).begin(); it != rois_v_loose.at(chid).end(); it++) {
                SignalROI *roi = *it;
                // initialize the extended bins ...
                roi->set_ext_start_bin(roi->get_start_bin());
                roi->set_ext_end_bin(roi->get_end_bin());

                if (flag) {
                    // loop through front
                    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int ext_start_bin = roi->get_ext_start_bin();
                        int ext_end_bin = roi->get_ext_end_bin();
                        if (ext_start_bin > roi1->get_start_bin()) ext_start_bin = roi1->get_start_bin();
                        if (ext_end_bin < roi1->get_end_bin()) ext_end_bin = roi1->get_end_bin();
                        roi->set_ext_start_bin(ext_start_bin);
                        roi->set_ext_end_bin(ext_end_bin);
                        //	std::cout << roi->get_chid() << " " << roi1->get_chid() << " " << roi->get_start_bin()
                        //<< " " << roi->get_end_bin() << " " << roi1->get_start_bin() << " " << roi1->get_end_bin() <<
                        //std::endl;
                    }

                    // loop through back
                    for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int ext_start_bin = roi->get_ext_start_bin();
                        int ext_end_bin = roi->get_ext_end_bin();
                        if (ext_start_bin > roi1->get_start_bin()) ext_start_bin = roi1->get_start_bin();
                        if (ext_end_bin < roi1->get_end_bin()) ext_end_bin = roi1->get_end_bin();
                        roi->set_ext_start_bin(ext_start_bin);
                        roi->set_ext_end_bin(ext_end_bin);
                        //	std::cout << roi->get_chid() << " " << roi1->get_chid() << " " << roi->get_start_bin()
                        //<< " " << roi->get_end_bin() << " " << roi1->get_start_bin() << " " << roi1->get_end_bin() <<
                        //std::endl;
                    }
                    // std::cout << chid << " " << roi->get_start_bin() << " " << roi->get_end_bin() << std::endl;
                }
            }

            if (flag) {
                SignalROI *prev_roi = 0;
                for (auto it = rois_v_loose.at(chid).begin(); it != rois_v_loose.at(chid).end(); it++) {
                    SignalROI *roi = *it;
                    if (prev_roi != 0) {
                        if (prev_roi->get_ext_end_bin() > roi->get_ext_start_bin()) {
                            prev_roi->set_ext_end_bin(int(prev_roi->get_end_bin() * 0.5 + roi->get_start_bin() * 0.5));
                            roi->set_ext_start_bin(int(prev_roi->get_end_bin() * 0.5 + roi->get_start_bin() * 0.5));
                        }
                    }
                    prev_roi = roi;
                }
            }
            // for (auto it = rois_v_loose.at(chid).begin(); it!= rois_v_loose.at(chid).end();it++){
            //   SignalROI *roi =  *it;
            //   std::cout << chid << " " << roi->get_ext_start_bin() << " " << roi->get_ext_end_bin() << std::endl;
            // }
        }
    }
    else {
        // W plane
        for (int chid = 0; chid != nwire_w; chid++) {
            rois_w_tight.at(chid).sort(CompareRois());
            for (auto it = rois_w_tight.at(chid).begin(); it != rois_w_tight.at(chid).end(); it++) {
                SignalROI *roi = *it;
                // initialize the extended bins ...
                roi->set_ext_start_bin(roi->get_start_bin());
                roi->set_ext_end_bin(roi->get_end_bin());

                if (flag) {
                    // loop through front
                    for (auto it1 = front_rois[roi].begin(); it1 != front_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int ext_start_bin = roi->get_ext_start_bin();
                        int ext_end_bin = roi->get_ext_end_bin();
                        if (ext_start_bin > roi1->get_start_bin()) ext_start_bin = roi1->get_start_bin();
                        if (ext_end_bin < roi1->get_end_bin()) ext_end_bin = roi1->get_end_bin();
                        roi->set_ext_start_bin(ext_start_bin);
                        roi->set_ext_end_bin(ext_end_bin);
                        //	std::cout << roi->get_chid() << " " << roi1->get_chid() << " " << roi->get_start_bin()
                        //<< " " << roi->get_end_bin() << " " << roi1->get_start_bin() << " " << roi1->get_end_bin() <<
                        //std::endl;
                    }

                    // loop through back
                    for (auto it1 = back_rois[roi].begin(); it1 != back_rois[roi].end(); it1++) {
                        SignalROI *roi1 = *it1;
                        int ext_start_bin = roi->get_ext_start_bin();
                        int ext_end_bin = roi->get_ext_end_bin();
                        if (ext_start_bin > roi1->get_start_bin()) ext_start_bin = roi1->get_start_bin();
                        if (ext_end_bin < roi1->get_end_bin()) ext_end_bin = roi1->get_end_bin();
                        roi->set_ext_start_bin(ext_start_bin);
                        roi->set_ext_end_bin(ext_end_bin);
                        //	std::cout << roi->get_chid() << " " << roi1->get_chid() << " " << roi->get_start_bin()
                        //<< " " << roi->get_end_bin() << " " << roi1->get_start_bin() << " " << roi1->get_end_bin() <<
                        //std::endl;
                    }
                    // std::cout << chid << " " << roi->get_start_bin() << " " << roi->get_end_bin() << std::endl;
                }
            }
            if (flag) {
                SignalROI *prev_roi = 0;
                for (auto it = rois_w_tight.at(chid).begin(); it != rois_w_tight.at(chid).end(); it++) {
                    SignalROI *roi = *it;
                    if (prev_roi != 0) {
                        if (prev_roi->get_ext_end_bin() > roi->get_ext_start_bin()) {
                            prev_roi->set_ext_end_bin(int(prev_roi->get_end_bin() * 0.5 + roi->get_start_bin() * 0.5));
                            roi->set_ext_start_bin(int(prev_roi->get_end_bin() * 0.5 + roi->get_start_bin() * 0.5));
                        }
                    }
                    prev_roi = roi;
                }
            }
            // for (auto it = rois_w_tight.at(chid).begin(); it!= rois_w_tight.at(chid).end();it++){
            //   SignalROI *roi =  *it;
            //   std::cout << chid << " " << roi->get_ext_start_bin() << " " << roi->get_ext_end_bin() << std::endl;
            // }
        }
    }
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
