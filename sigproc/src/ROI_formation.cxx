
#include "ROI_formation.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::SigProc;

ROI_formation::ROI_formation(Waveform::ChannelMaskMap& cmm, int nwire_u, int nwire_v, int nwire_w, int nbins,
                             float th_factor_ind, float th_factor_col, int pad, float asy, int rebin, double l_factor,
                             double l_max_th, double l_factor1, int l_short_length, int l_jump_one_bin)
  : nwire_u(nwire_u)
  , nwire_v(nwire_v)
  , nwire_w(nwire_w)
  , nbins(nbins)
  , th_factor_ind(th_factor_ind)
  , th_factor_col(th_factor_col)
  , pad(pad)
  , asy(asy)
  , rebin(rebin)
  , l_factor(l_factor)
  , l_max_th(l_max_th)
  , l_factor1(l_factor1)
  , l_short_length(l_short_length)
  , l_jump_one_bin(l_jump_one_bin)
{
    self_rois_u.resize(nwire_u);
    self_rois_v.resize(nwire_v);
    self_rois_w.resize(nwire_w);

    loose_rois_u.resize(nwire_u);
    loose_rois_v.resize(nwire_v);
    loose_rois_w.resize(nwire_w);

    uplane_rms.resize(nwire_u);
    vplane_rms.resize(nwire_v);
    wplane_rms.resize(nwire_w);

    for (auto it = cmm["bad"].begin(); it != cmm["bad"].end(); it++) {
        int ch = it->first;
        std::vector<std::pair<int, int>> temps;
        bad_ch_map[ch] = temps;
        for (size_t ind = 0; ind < it->second.size(); ind++) {
            bad_ch_map[ch].push_back(std::make_pair(it->second.at(ind).first, it->second.at(ind).second));
            // std::cout << ch << " " <<  << std::endl;
        }
    }
    // std::cout << bad_ch_map.size() << std::endl;
}

void ROI_formation::apply_roi(int plane, Array::array_xxf& r_data, int flag)
{
    if (plane == 0) {
        for (int irow = 0; irow != r_data.rows(); irow++) {
            // refresh ...
            Waveform::realseq_t signal(r_data.cols(), 0);
            // loop ROIs and assign data

            std::vector<std::pair<int, int>> rois;
            if (flag == 1) {  // self ROI
                rois = self_rois_u.at(irow);
            }
            else if (flag == 2) {  // loose ROI
                rois = loose_rois_u.at(irow);
            }

            for (auto it = rois.begin(); it != rois.end(); it++) {
                int start_bin = it->first;
                int end_bin = it->second;

                float start_content = r_data(irow, start_bin);
                float end_content = r_data(irow, end_bin);
                for (int i = start_bin; i < end_bin + 1; i++) {
                    int content =
                        r_data(irow, i) -
                        ((end_content - start_content) * (i - start_bin) / (end_bin - start_bin) + start_content);
                    signal.at(i) = content;
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
            std::vector<std::pair<int, int>> rois;
            if (flag == 1) {  // self ROI
                rois = self_rois_v.at(irow);
            }
            else if (flag == 2) {  // loose ROI
                rois = loose_rois_v.at(irow);
            }

            for (auto it = rois.begin(); it != rois.end(); it++) {
                int start_bin = it->first;
                int end_bin = it->second;

                float start_content = r_data(irow, start_bin);
                float end_content = r_data(irow, end_bin);
                for (int i = start_bin; i < end_bin + 1; i++) {
                    int content =
                        r_data(irow, i) -
                        ((end_content - start_content) * (i - start_bin) / (end_bin - start_bin) + start_content);
                    signal.at(i) = content;
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
            std::vector<std::pair<int, int>> rois;
            if (flag == 1) {  // self ROI
                rois = self_rois_w.at(irow);
            }
            else if (flag == 2) {  // loose ROI
                rois = loose_rois_w.at(irow);
            }

            // loop ROIs and assign data
            for (auto it = rois.begin(); it != rois.end(); it++) {
                int start_bin = it->first;
                int end_bin = it->second;

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

ROI_formation::~ROI_formation() {}

void ROI_formation::Clear()
{
    self_rois_u.clear();
    self_rois_v.clear();
    self_rois_w.clear();

    loose_rois_u.clear();
    loose_rois_v.clear();
    loose_rois_w.clear();

    uplane_rms.clear();
    vplane_rms.clear();
    wplane_rms.clear();
}

void ROI_formation::extend_ROI_self(int plane)
{
    if (plane == 0) {
        for (size_t i = 0; i != self_rois_u.size(); i++) {
            std::vector<std::pair<int, int>> temp_rois;
            int temp_begin = 0, temp_end = 0;
            for (size_t j = 0; j != self_rois_u.at(i).size(); j++) {
                temp_begin = self_rois_u.at(i).at(j).first - pad;
                if (temp_begin < 0) temp_begin = 0;
                temp_end = self_rois_u.at(i).at(j).second + pad;
                if (temp_end >= nbins) temp_end = nbins - 1;
                // merge
                if (temp_rois.size() == 0) {
                    temp_rois.push_back(std::make_pair(temp_begin, temp_end));
                }
                else {
                    if (temp_begin < temp_rois.back().second) {
                        if (temp_end > temp_rois.back().second) {
                            temp_rois.back().second = temp_end;
                        }
                    }
                    else {
                        temp_rois.push_back(std::make_pair(temp_begin, temp_end));
                    }
                }
            }
            self_rois_u.at(i) = temp_rois;
        }
    }
    else if (plane == 1) {
        for (size_t i = 0; i != self_rois_v.size(); i++) {
            std::vector<std::pair<int, int>> temp_rois;
            int temp_begin = 0, temp_end = 0;
            for (size_t j = 0; j != self_rois_v.at(i).size(); j++) {
                temp_begin = self_rois_v.at(i).at(j).first - pad;
                if (temp_begin < 0) temp_begin = 0;
                temp_end = self_rois_v.at(i).at(j).second + pad;
                if (temp_end >= nbins) temp_end = nbins - 1;
                // merge
                if (temp_rois.size() == 0) {
                    temp_rois.push_back(std::make_pair(temp_begin, temp_end));
                }
                else {
                    if (temp_begin < temp_rois.back().second) {
                        if (temp_end > temp_rois.back().second) {
                            temp_rois.back().second = temp_end;
                        }
                    }
                    else {
                        temp_rois.push_back(std::make_pair(temp_begin, temp_end));
                    }
                }
            }
            self_rois_v.at(i) = temp_rois;
        }
    }
    else if (plane == 2) {
        for (size_t i = 0; i != self_rois_w.size(); i++) {
            std::vector<std::pair<int, int>> temp_rois;
            int temp_begin = 0, temp_end = 0;
            for (size_t j = 0; j != self_rois_w.at(i).size(); j++) {
                temp_begin = self_rois_w.at(i).at(j).first - pad;
                if (temp_begin < 0) temp_begin = 0;
                temp_end = self_rois_w.at(i).at(j).second + pad;
                if (temp_end >= nbins) temp_end = nbins - 1;
                // merge
                if (temp_rois.size() == 0) {
                    temp_rois.push_back(std::make_pair(temp_begin, temp_end));
                }
                else {
                    if (temp_begin < temp_rois.back().second) {
                        if (temp_end > temp_rois.back().second) {
                            temp_rois.back().second = temp_end;
                        }
                    }
                    else {
                        temp_rois.push_back(std::make_pair(temp_begin, temp_end));
                    }
                }
            }
            self_rois_w.at(i) = temp_rois;
        }
    }
}

void ROI_formation::create_ROI_connect_info(int plane)
{
    if (plane == 0) {
        // u
        for (int i = 0; i != nwire_u - 2; i++) {
            for (size_t j = 0; j != self_rois_u.at(i).size(); j++) {
                int start1 = self_rois_u.at(i).at(j).first;
                int end1 = self_rois_u.at(i).at(j).second;
                int length1 = end1 - start1 + 1;
                for (size_t k = 0; k != self_rois_u.at(i + 2).size(); k++) {
                    int start2 = self_rois_u.at(i + 2).at(k).first;
                    int end2 = self_rois_u.at(i + 2).at(k).second;
                    int length2 = end2 - start2 + 1;
                    if (std::abs(length2 - length1) < (length2 + length1) * asy) {
                        int start3 = (start1 + start2) / 2.;
                        int end3 = (end1 + end2) / 2.;
                        if (start3 < end3 && start3 <= end1 && start3 <= end2 && end3 >= start1 && end3 >= start2) {
                            // go through existing ones to make sure there is no overlap
                            int flag = 0;
                            for (size_t i1 = 0; i1 != self_rois_u.at(i + 1).size(); i1++) {
                                int max_start = start3;
                                if (self_rois_u.at(i + 1).at(i1).first > max_start)
                                    max_start = self_rois_u.at(i + 1).at(i1).first;
                                int min_end = end3;
                                if (self_rois_u.at(i + 1).at(i1).second < min_end)
                                    min_end = self_rois_u.at(i + 1).at(i1).second;
                                if (max_start < min_end) {
                                    flag = 1;
                                    break;
                                }
                            }
                            if (flag == 0) self_rois_u.at(i + 1).push_back(std::make_pair(start3, end3));
                        }
                    }
                }
            }
        }
    }
    else if (plane == 1) {
        // v
        for (int i = 0; i != nwire_v - 2; i++) {
            for (size_t j = 0; j != self_rois_v.at(i).size(); j++) {
                int start1 = self_rois_v.at(i).at(j).first;
                int end1 = self_rois_v.at(i).at(j).second;
                int length1 = end1 - start1 + 1;
                for (size_t k = 0; k != self_rois_v.at(i + 2).size(); k++) {
                    int start2 = self_rois_v.at(i + 2).at(k).first;
                    int end2 = self_rois_v.at(i + 2).at(k).second;
                    int length2 = end2 - start2 + 1;
                    if (std::abs(length2 - length1) < (length2 + length1) * asy) {
                        int start3 = (start1 + start2) / 2.;
                        int end3 = (end1 + end2) / 2.;
                        if (start3 < end3 && start3 <= end1 && start3 <= end2 && end3 >= start1 && end3 >= start2) {
                            // go through existing ones to make sure there is no overlap
                            int flag = 0;
                            for (size_t i1 = 0; i1 != self_rois_v.at(i + 1).size(); i1++) {
                                int max_start = start3;
                                if (self_rois_v.at(i + 1).at(i1).first > max_start)
                                    max_start = self_rois_v.at(i + 1).at(i1).first;
                                int min_end = end3;
                                if (self_rois_v.at(i + 1).at(i1).second < min_end)
                                    min_end = self_rois_v.at(i + 1).at(i1).second;
                                if (max_start < min_end) {
                                    flag = 1;
                                    break;
                                }
                            }
                            if (flag == 0) self_rois_v.at(i + 1).push_back(std::make_pair(start3, end3));
                        }
                    }
                }
            }
        }
    }
    else if (plane == 2) {
        // w?
        for (int i = 0; i != nwire_w - 2; i++) {
            for (size_t j = 0; j != self_rois_w.at(i).size(); j++) {
                int start1 = self_rois_w.at(i).at(j).first;
                int end1 = self_rois_w.at(i).at(j).second;
                int length1 = end1 - start1 + 1;
                for (size_t k = 0; k != self_rois_w.at(i + 2).size(); k++) {
                    int start2 = self_rois_w.at(i + 2).at(k).first;
                    int end2 = self_rois_w.at(i + 2).at(k).second;
                    int length2 = end2 - start2 + 1;
                    if (std::abs(length2 - length1) < (length2 + length1) * asy) {
                        int start3 = (start1 + start2) / 2.;
                        int end3 = (end1 + end2) / 2.;
                        if (start3 < end3 && start3 <= end1 && start3 <= end2 && end3 >= start1 && end3 >= start2) {
                            // go through existing ones to make sure there is no overlap
                            int flag = 0;
                            for (size_t i1 = 0; i1 != self_rois_w.at(i + 1).size(); i1++) {
                                int max_start = start3;
                                if (self_rois_w.at(i + 1).at(i1).first > max_start)
                                    max_start = self_rois_w.at(i + 1).at(i1).first;
                                int min_end = end3;
                                if (self_rois_w.at(i + 1).at(i1).second < min_end)
                                    min_end = self_rois_w.at(i + 1).at(i1).second;
                                if (max_start < min_end) {
                                    flag = 1;
                                    break;
                                }
                            }
                            if (flag == 0) self_rois_w.at(i + 1).push_back(std::make_pair(start3, end3));
                        }
                    }
                }
            }
        }
    }
}

double ROI_formation::cal_RMS(Waveform::realseq_t signal)
{
    double result = 0;
    if (signal.size() > 0) {
        // do quantile ...
        float par[3];
        par[0] = WireCell::Waveform::percentile(signal, 0.5 - 0.34);
        par[1] = WireCell::Waveform::percentile(signal, 0.5);
        par[2] = WireCell::Waveform::percentile(signal, 0.5 + 0.34);
        float rms = sqrt((pow(par[2] - par[1], 2) + pow(par[1] - par[0], 2)) / 2.);

        float rms2 = 0;
        float rms1 = 0;
        for (size_t i = 0; i != signal.size(); i++) {
            if (fabs(signal.at(i)) < 5.0 * rms) {
                rms1 += pow(signal.at(i), 2);
                rms2++;
            }
        }
        if (rms2 > 0) {
            result = sqrt(rms1 / rms2);
        }
    }

    return result;
}

void ROI_formation::find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data,
                                             const Array::array_xxf& r_data_tight)
{
    int offset = 0;
    if (plane == 0) {
        offset = 0;
    }
    else if (plane == 1) {
        offset = nwire_u;
    }
    else if (plane == 2) {
        offset = nwire_u + nwire_v;
    }

    for (int irow = 0; irow != r_data.rows(); irow++) {
        // calclulate rms for a row of r_data
        Waveform::realseq_t signal(nbins);
        Waveform::realseq_t signal1(nbins);
        Waveform::realseq_t signal2(nbins);

        if (bad_ch_map.find(irow + offset) != bad_ch_map.end()) {
            int ncount = 0;
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
                    signal.at(ncount) = r_data(irow, icol);
                    signal1.at(icol) = r_data(irow, icol);
                    signal2.at(icol) = r_data_tight(irow, icol);
                    ncount++;
                }
                else {
                    signal1.at(icol) = 0;
                    signal2.at(icol) = 0;
                }
            }
            signal.resize(ncount);
        }
        else {
            for (int icol = 0; icol != r_data.cols(); icol++) {
                signal.at(icol) = r_data(irow, icol);
                signal1.at(icol) = r_data(irow, icol);
                signal2.at(icol) = r_data_tight(irow, icol);
            }
        }
        // do threshold and fill rms
        double rms = cal_RMS(signal);
        double threshold = 0;
        if (plane == 0) {
            threshold = th_factor_ind * rms + 1;
            uplane_rms.at(irow) = rms;
        }
        else if (plane == 1) {
            threshold = th_factor_ind * rms + 1;
            vplane_rms.at(irow) = rms;
        }
        else if (plane == 2) {
            threshold = th_factor_col * rms + 1;
            wplane_rms.at(irow) = rms;
        }

        //  std::cout << plane << " " << signal.size() << " " << irow << " " << rms << std::endl;

        // create rois
        int roi_begin = -1;
        int roi_end = -1;

        std::vector<std::pair<int, int>> temp_rois;
        // now find ROI, above five sigma, and pad with +- six time ticks
        for (int j = 0; j < int(signal1.size()) - 1; j++) {
            double content = signal1.at(j);
            double content_tight = signal2.at(j);

            if (content > threshold || (content_tight > threshold)) {
                roi_begin = j;
                roi_end = j;
                for (int k = j + 1; k < int(signal1.size()); k++) {
                    if (signal1.at(k) > threshold || (signal2.at(k) > threshold)) {
                        roi_end = k;
                    }
                    else {
                        break;
                    }
                }
                int temp_roi_begin = roi_begin;  // filter_pad;
                if (temp_roi_begin < 0) temp_roi_begin = 0;
                int temp_roi_end = roi_end;  // filter_pad;
                if (temp_roi_end > int(signal1.size()) - 1) temp_roi_end = int(signal1.size()) - 1;

                //	if (abs(irow-1199)<=1&&plane==0) std::cout << "Tight: " << irow << " " << temp_roi_begin << " "
                //<< temp_roi_end << std::endl;

                if (temp_rois.size() == 0) {
                    temp_rois.push_back(std::make_pair(temp_roi_begin, temp_roi_end));
                }
                else {
                    if (temp_roi_begin <= temp_rois.back().second) {
                        temp_rois.back().second = temp_roi_end;
                    }
                    else {
                        temp_rois.push_back(std::make_pair(temp_roi_begin, temp_roi_end));
                    }
                }
                j = roi_end + 1;
            }
        }

        // if (plane==2 && irow == 69){
        //   std::cout << "Xin: " << irow << " " << rms << " " << temp_rois.size() << std::endl;
        //   for (size_t i=0;i!=temp_rois.size();i++){
        // 	std::cout << "Xin: " << temp_rois.at(i).first << " " << temp_rois.at(i).second << std::endl;
        //   }
        // }

        // fill rois ...
        if (plane == 0) {
            self_rois_u.at(irow) = temp_rois;
        }
        else if (plane == 1) {
            self_rois_v.at(irow) = temp_rois;
        }
        else {
            self_rois_w.at(irow) = temp_rois;
        }
        //    std::cout << plane << " " << irow << " " << temp_rois.size() << std::endl;
    }

    extend_ROI_self(plane);

    // for (size_t i=0;i!=self_rois_w.at(69).size();i++){
    //   std::cout << "Xin: " << self_rois_w.at(69).at(i).first << " " << self_rois_w.at(69).at(i).second << std::endl;
    // }

    create_ROI_connect_info(plane);
}

void ROI_formation::find_ROI_by_decon_itself(int plane, const Array::array_xxf& r_data)
{
    find_ROI_by_decon_itself(plane, r_data, r_data);
}

void ROI_formation::extend_ROI_loose(int plane)
{
    if (plane == 0) {
        // compare the loose one with tight one
        for (int i = 0; i != nwire_u; i++) {
            std::vector<std::pair<int, int>> temp_rois;
            for (size_t j = 0; j != loose_rois_u.at(i).size(); j++) {
                int start = loose_rois_u.at(i).at(j).first;
                int end = loose_rois_u.at(i).at(j).second;
                for (size_t k = 0; k != self_rois_u.at(i).size(); k++) {
                    int temp_start = self_rois_u.at(i).at(k).first;
                    int temp_end = self_rois_u.at(i).at(k).second;
                    if (start > temp_start && start < temp_end) start = temp_start;
                    // loop through all the tight one to examine start
                    if (end > temp_start && end < temp_end) end = temp_end;
                    // loop through all the tight one to examine the end
                }
                if (temp_rois.size() == 0) {
                    temp_rois.push_back(std::make_pair(start, end));
                }
                else {
                    if (start < temp_rois.back().second) {
                        temp_rois.back().second = end;
                    }
                    else {
                        temp_rois.push_back(std::make_pair(start, end));
                    }
                }
            }
            loose_rois_u.at(i) = temp_rois;
        }
    }
    else if (plane == 1) {
        for (int i = 0; i != nwire_v; i++) {
            std::vector<std::pair<int, int>> temp_rois;
            for (size_t j = 0; j != loose_rois_v.at(i).size(); j++) {
                int start = loose_rois_v.at(i).at(j).first;
                int end = loose_rois_v.at(i).at(j).second;
                for (size_t k = 0; k != self_rois_v.at(i).size(); k++) {
                    int temp_start = self_rois_v.at(i).at(k).first;
                    int temp_end = self_rois_v.at(i).at(k).second;
                    if (start > temp_start && start < temp_end) start = temp_start;
                    // loop through all the tight one to examine start
                    if (end > temp_start && end < temp_end) end = temp_end;
                    // loop through all the tight one to examine the end
                }
                if (temp_rois.size() == 0) {
                    temp_rois.push_back(std::make_pair(start, end));
                }
                else {
                    if (start < temp_rois.back().second) {
                        temp_rois.back().second = end;
                    }
                    else {
                        temp_rois.push_back(std::make_pair(start, end));
                    }
                }
            }
            loose_rois_v.at(i) = temp_rois;
        }
    }
}

double ROI_formation::local_ave(Waveform::realseq_t& signal, int bin, int width)
{
    double sum1 = 0;
    double sum2 = 0;

    for (int i = -width; i < width + 1; i++) {
        int current_bin = bin + i;

        while (current_bin < 0) current_bin += signal.size();
        while (current_bin >= int(signal.size())) current_bin -= signal.size();

        sum1 += signal.at(current_bin);
        sum2++;
    }

    if (sum2 > 0) {
        return sum1 / sum2;
    }
    else {
        return 0;
    }
}

int ROI_formation::find_ROI_end(Waveform::realseq_t& signal, int bin, double th, int jump_one_bin)
{
    int end = bin;
    double content = signal.at(end);
    while (content > th) {
        end++;
        if (end >= int(signal.size())) {
            content = signal.at(end - signal.size());
        }
        else {
            content = signal.at(end);
        }
        if (end >= int(signal.size()) - 1) {
            end = int(signal.size()) - 1;
            break;
        }
    }

    while (local_ave(signal, end + 1, 1) < local_ave(signal, end, 1) + 25 ||
           (jump_one_bin && local_ave(signal, end + 2, 1) < local_ave(signal, end, 1) + 25)) {  // ||
        //	(local_ave(signal,end+1,1) + local_ave(signal,end+2,1))*0.5 < local_ave(signal,end,1) ){
        end++;
        if (end >= int(signal.size()) - 1) {
            end = int(signal.size()) - 1;
            break;
        }
    }
    return end;
}
int ROI_formation::find_ROI_begin(Waveform::realseq_t& signal, int bin, double th, int jump_one_bin)
{
    // find the first one before bin and is below threshold ...
    int begin = bin;
    double content = signal.at(begin);
    while (content > th) {
        begin--;
        if (begin < 0) {
            content = signal.at(begin + int(signal.size()));
        }
        else {
            content = signal.at(begin);
        }
        if (begin <= 0) {
            begin = 0;
            break;
        }
    }

    // calculate the local average
    // keep going and find the minimum
    while (local_ave(signal, begin - 1, 1) < local_ave(signal, begin, 1) + 25 ||
           (jump_one_bin && local_ave(signal, begin - 2, 1) < local_ave(signal, begin, 1) + 25)) {  // ||
        // (local_ave(signal,begin-2,1) + local_ave(signal,begin-1,1))*0.5 < local_ave(signal,begin,1) ){
        begin--;
        if (begin <= 0) {
            begin = 0;
            break;
        }
    }

    return begin;
}

void ROI_formation::find_ROI_loose(int plane, const Array::array_xxf& r_data)
{
    int offset = 0;
    if (plane == 0) {
        offset = 0;
    }
    else if (plane == 1) {
        offset = nwire_u;
    }
    else if (plane == 2) {
        offset = nwire_u + nwire_v;
    }

    // form rebinned waveform ...
    for (int irow = 0; irow != r_data.rows(); irow++) {
        Waveform::realseq_t signal(nbins);                // remove bad ones
        Waveform::realseq_t signal1(nbins);               // all signal
        Waveform::realseq_t signal2(int(nbins / rebin));  // rebinned ones

        // std::cout << "xin1" << std::endl;

        if (bad_ch_map.find(irow + offset) != bad_ch_map.end()) {
            int ncount = 0;
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
                    signal.at(ncount) = r_data(irow, icol);
                    signal1.at(icol) = r_data(irow, icol);
                    ncount++;
                }
                else {
                    signal1.at(icol) = 0;
                }
            }
            signal.resize(ncount);
        }
        else {
            for (int icol = 0; icol != r_data.cols(); icol++) {
                signal.at(icol) = r_data(irow, icol);
                signal1.at(icol) = r_data(irow, icol);
            }
        }

        // std::cout << "xin2" << std::endl;

        // get rebinned waveform
        for (size_t i = 0; i != signal2.size(); i++) {
            double temp = 0;
            for (int j = 0; j != rebin; j++) {
                temp += signal1.at(rebin * i + j);
            }
            signal2.at(i) = temp;
        }

        // std::cout << "xin3" << " " << signal.size() << " " << signal2.size() << std::endl;

        // calculate rms ...
        float th = cal_RMS(signal) * rebin * l_factor;
        // if (irow==1240) std::cout << "a " << l_factor << " " << rebin << " " << th << " " << l_max_th << std::endl;
        if (th > l_max_th) th = l_max_th;

        std::vector<std::pair<int, int>> ROIs_1;
        std::vector<int> max_bins_1;
        int ntime = signal2.size();

        // if (irow == 1200 && plane==0){
        //   for (int j=0;j!=ntime;j++){
        // 	std::cout << j << " " << signal2.at(j) << " " << th << " " << l_factor1 << " " << std::endl;
        //   }
        // }

        for (int j = 1; j < ntime - 1; j++) {
            double content = signal2.at(j);
            double prev_content = signal2.at(j - 1);
            double next_content = signal2.at(j + 1);
            int flag_ROI = 0;
            int begin = 0;
            int end = 0;
            int max_bin = 0;
            if (content > th) {
                begin = find_ROI_begin(signal2, j, th * l_factor1, l_jump_one_bin);
                end = find_ROI_end(signal2, j, th * l_factor1, l_jump_one_bin);
                max_bin = begin;
                //	if (irow==1240) std::cout << "a: " << begin << " " << end << " " << j << std::endl;
                for (int k = begin; k <= end; k++) {
                    // std::cout << begin << " " << end << " " << max_bin << std::endl;
                    if (signal2.at(k) > signal2.at(max_bin)) {
                        max_bin = k;
                    }
                }
                flag_ROI = 1;
            }
            else {
                if (content > prev_content && content > next_content) {
                    begin = find_ROI_begin(signal2, j, prev_content, l_jump_one_bin);
                    end = find_ROI_end(signal2, j, next_content, l_jump_one_bin);
                    max_bin = begin;
                    for (int k = begin; k <= end; k++) {
                        if (signal2.at(k) > signal2.at(max_bin)) {
                            max_bin = k;
                        }
                    }
                    if (signal2.at(max_bin) - signal2.at(begin) + signal2.at(max_bin) - signal2.at(end) > th * 2) {
                        flag_ROI = 1;
                    }

                    int temp_begin = max_bin - l_short_length;
                    if (temp_begin < begin) temp_begin = begin;
                    int temp_end = max_bin + l_short_length;
                    if (temp_end > end) temp_end = end;
                    if ((signal2.at(max_bin) - signal2.at(temp_begin) > th * l_factor1 &&
                         signal2.at(max_bin) - signal2.at(temp_end) > th * l_factor1)) {
                        flag_ROI = 1;
                    }

                    // if (irow==1200 && plane==0)
                    //   std::cout << j << " " << begin << " " << end << " " << max_bin << " " << signal2.at(max_bin) *
                    //   2 - signal2.at(begin) - signal2.at(end) << " " << th*2 << " " << temp_begin << " " << temp_end
                    //   << " " << signal2.at(max_bin) - signal2.at(temp_begin) << " " << signal2.at(max_bin) -
                    //   signal2.at(temp_end) << " " << th*l_factor1 << " " << flag_ROI << std::endl;
                }
            }

            if (flag_ROI == 1) {
                // if (irow==1240) {
                //   std::cout << begin << " " << end << " " << j << " " << content << " " << th << " " << prev_content
                //   << " " << next_content << std::endl; for (int kk = begin; kk!=end; kk++){
                //     std::cout << kk << " a " << signal2.at(kk) << " " << th*l_factor1 << " " <<
                //     local_ave(signal2,kk,1) << std::endl; for (int kkk=0;kkk!=6;kkk++){
                //       std::cout << "b " << signal.at(kk*6+kkk) << std::endl;
                //     }
                //   }
                // }

                if (ROIs_1.size() > 0) {
                    if (begin <= ROIs_1.back().second) {
                        ROIs_1.back().second = end;
                        if (signal2.at(max_bin) > signal2.at(max_bins_1.back())) max_bins_1.back() = max_bin;
                    }
                    else {
                        ROIs_1.push_back(std::make_pair(begin, end));
                        max_bins_1.push_back(max_bin);
                    }
                }
                else {
                    ROIs_1.push_back(std::make_pair(begin, end));
                    max_bins_1.push_back(max_bin);
                }

                if (end < int(signal2.size())) {
                    j = end;
                }
                else {
                    j = signal2.size();
                }
            }
        }

        // std::cout << "xin4" << std::endl;

        // for (int j = 0; j!=int(ROIs_1.size());j++){
        //    int begin = ROIs_1.at(j).first * rebin;
        //    int end = ROIs_1.at(j).second *rebin + (rebin-1);
        //    if (irow ==1240) std::cout << "b " << begin << " " << end << std::endl;
        //  }

        if (ROIs_1.size() == 1) {
        }
        else if (ROIs_1.size() > 1) {
            int flag_repeat = 0;
            //  cout << "Xin1: " << ROIs_1.size() << endl;;
            while (flag_repeat) {
                flag_repeat = 1;
                for (int k = 0; k < int(ROIs_1.size() - 1); k++) {
                    int begin = ROIs_1.at(k).first;
                    int end = ROIs_1.at(k + 1).second;

                    double begin_content = signal2.at(begin);
                    double end_content = signal2.at(end);

                    int begin_1 = ROIs_1.at(k).second;
                    int end_1 = ROIs_1.at(k + 1).first;

                    int flag_merge = 1;
                    // Double_t sum1 = 0, sum2 = 0;
                    for (int j = begin_1; j <= end_1; j++) {
                        double current_content = signal2.at(j);
                        double content =
                            current_content -
                            ((end_content - begin_content) * (j * 1. - begin) / (end - begin * 1.) + begin_content);

                        if (content < th * l_factor1) {
                            flag_merge = 0;
                            break;
                        }
                        // sum1 += content;
                        // sum2 ++;
                        // cout << j << " " << content << endl;
                    }
                    // if (sum2 >0){
                    //   if (sum1/sum2 < th*factor1) flag_merge = 0;
                    // }

                    if (flag_merge == 1) {
                        ROIs_1.at(k).second = ROIs_1.at(k + 1).second;
                        ROIs_1.erase(ROIs_1.begin() + k + 1);
                        flag_repeat = 1;
                        break;
                    }
                }
                //	cout << "Xin2: " << ROIs_1.size() << endl;
            }
        }

        // std::cout << "xin5" << std::endl;
        // scale back ...
        for (int j = 0; j != int(ROIs_1.size()); j++) {
            int begin = ROIs_1.at(j).first * rebin;
            int end = ROIs_1.at(j).second * rebin + (rebin - 1);

            ROIs_1.at(j).first = begin;
            ROIs_1.at(j).second = end;

            // if (abs(irow-1199)<=1&& plane==0) std::cout << "Loose: "  << irow << " " << ROIs_1.at(j).first << " " <<
            // ROIs_1.at(j).second << std::endl;
        }

        if (plane == 0) {
            loose_rois_u.at(irow) = ROIs_1;
        }
        else if (plane == 1) {
            loose_rois_v.at(irow) = ROIs_1;
        }
        else {
            loose_rois_w.at(irow) = ROIs_1;
        }
        // std::cout << plane << " " << irow << " " << ROIs_1.size() << std::endl;
    }

    //  std::cout << "xin6" << std::endl;

    extend_ROI_loose(plane);
    // std::cout << "xin7" << std::endl;
}
