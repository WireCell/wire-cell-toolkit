#include "WireCellSigProc/L1SPFilter.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellIface/IFieldResponse.h"

#include "WireCellUtil/LassoModel.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Eigen.h"

#include <numeric>
#include <iostream>

WIRECELL_FACTORY(L1SPFilter, WireCell::SigProc::L1SPFilter, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace Eigen;
using namespace WireCell;
using namespace WireCell::SigProc;

using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;


L1SPFilter::L1SPFilter(double gain, double shaping, double postgain, double ADC_mV, double fine_time_offset,
                       double coarse_time_offset)
  : m_gain(gain)
  , m_shaping(shaping)
  , m_postgain(postgain)
  , m_ADC_mV(ADC_mV)
  , m_fine_time_offset(fine_time_offset)
  , m_coarse_time_offset(coarse_time_offset)
  , lin_V(0)
  , lin_W(0)
{
}

L1SPFilter::~L1SPFilter()
{
    delete lin_V;
    delete lin_W;
}

void L1SPFilter::init_resp()
{
    if (lin_V == 0 && lin_W == 0) {
        // get field response ...
        auto ifr = Factory::find_tn<IFieldResponse>(get<std::string>(m_cfg, "fields", "FieldResponse"));
        Response::Schema::FieldResponse fr = ifr->field_response();
        // Make a new data set which is the average FR, make an average for V and Y planes ...
        Response::Schema::FieldResponse fravg = Response::average_1D(fr);

        // get electronics response
        WireCell::Waveform::compseq_t elec;
        WireCell::Binning tbins(Response::as_array(fravg.planes[0]).cols(), 0,
                                Response::as_array(fravg.planes[0]).cols() * fravg.period);
        Response::ColdElec ce(m_gain, m_shaping);
        auto ewave = ce.generate(tbins);
        Waveform::scale(ewave, m_postgain * m_ADC_mV * (-1));  // ADC to electron ...
        elec = fwd_r2c(m_dft, ewave); 

        std::complex<float> fine_period(fravg.period, 0);

        // do a convolution here ...
        WireCell::Waveform::realseq_t resp_V = fravg.planes[1].paths[0].current;
        WireCell::Waveform::realseq_t resp_W = fravg.planes[2].paths[0].current;

        auto spectrum_V = fwd_r2c(m_dft, resp_V);
        auto spectrum_W = fwd_r2c(m_dft, resp_W);

        WireCell::Waveform::scale(spectrum_V, elec);
        WireCell::Waveform::scale(spectrum_W, elec);

        WireCell::Waveform::scale(spectrum_V, fine_period);
        WireCell::Waveform::scale(spectrum_W, fine_period);

        // Now this response is ADC for 1 electron .
        resp_V = inv_c2r(m_dft, spectrum_V);
        resp_W = inv_c2r(m_dft, spectrum_W);

        // convolute with V and Y average responses ...
        double intrinsic_time_offset = fravg.origin / fravg.speed;
        // std::cout << intrinsic_time_offset << " " << m_fine_time_offset << " " << m_coarse_time_offset << " " <<
        // m_gain << " " << 14.0 * units::mV/units::fC << " " << m_shaping << " " << fravg.period << std::endl;

        double x0 = (-intrinsic_time_offset - m_coarse_time_offset + m_fine_time_offset);
        double xstep = fravg.period;

        // boost::math::cubic_b_spline<double> spline_v(resp_V.begin(), resp_V.end(), x0, xstep);
        // boost::math::cubic_b_spline<double> spline_w(resp_W.begin(), resp_W.end(), x0, xstep);
        lin_V = new linterp<double>(resp_V.begin(), resp_V.end(), x0, xstep);
        lin_W = new linterp<double>(resp_W.begin(), resp_W.end(), x0, xstep);
    }
}

WireCell::Configuration L1SPFilter::default_configuration() const
{
    Configuration cfg;

    /// Name of component providing field responses
    cfg["fields"] = "FieldResponse";

    /// An array holding a waveform to use as the "smearing" filter.
    cfg["filter"] = Json::arrayValue;

    /// The tag identifying traces which represent "raw" (not
    /// deconvolved) ADC values.
    cfg["adctag"] = "raw";

    /// The tag identifying traces which represent "signal" processed
    /// (deconvolved) waveforms.
    cfg["sigtag"] = "gauss";

    /// The tag to place on the output waveforms
    cfg["outtag"] = "l1sp";

    // 4 sigma for raw waveform ROI identification
    cfg["raw_ROI_th_nsigma"] = 4;
    // 10 ADC for upper limit on ADC ...
    cfg["raw_ROI_th_adclimit"] = 10;
    // global offset
    cfg["overall_time_offset"] = 0;
    // need 3 us offset for collection plane relative to the induction plane ...
    cfg["collect_time_offset"] = 3.0;

    // ROI padding ticks ...
    cfg["roi_pad"] = 3;
    cfg["raw_pad"] = 15;

    // L1 fit parameters ...
    cfg["adc_l1_threshold"] = 6;
    cfg["adc_sum_threshold"] = 160;
    cfg["adc_sum_rescaling"] = 90.;
    cfg["adc_sum_rescaling_limit"] = 50.;
    cfg["adc_ratio_threshold"] = 0.2;

    cfg["l1_seg_length"] = 120;
    cfg["l1_scaling_factor"] = 500;
    cfg["l1_lambda"] = 5;
    cfg["l1_epsilon"] = 0.05;
    cfg["l1_niteration"] = 100000;
    cfg["l1_decon_limit"] = 100;  // 100 electrons

    cfg["l1_resp_scale"] = 0.5;
    cfg["l1_col_scale"] = 1.15;
    cfg["l1_ind_scale"] = 0.5;

    cfg["peak_threshold"] = 1000;
    cfg["mean_threshold"] = 500;

    cfg["gain"] = m_gain;
    cfg["shaping"] = m_shaping;
    cfg["postgain"] = m_postgain;
    cfg["ADC_mV"] = m_ADC_mV;

    cfg["fine_time_offset"] = m_fine_time_offset;
    cfg["coarse_time_offset"] = m_coarse_time_offset;

    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use

    return cfg;
}

void L1SPFilter::configure(const WireCell::Configuration& cfg)
{
    m_cfg = cfg;

    m_gain = get(cfg, "gain", m_gain);
    m_shaping = get(cfg, "shaping", m_shaping);
    m_postgain = get(cfg, "postgain", m_postgain);
    m_ADC_mV = get(cfg, "ADC_mV", m_ADC_mV);

    m_fine_time_offset = get(cfg, "fine_time_offset", m_fine_time_offset);
    m_coarse_time_offset = get(cfg, "coarse_time_offset", m_coarse_time_offset);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
}

bool L1SPFilter::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        return true;  // eos
    }

    std::string adctag = get<std::string>(m_cfg, "adctag");
    std::string sigtag = get<std::string>(m_cfg, "sigtag");
    std::string outtag = get<std::string>(m_cfg, "outtag");

    //    std::cout << smearing_vec.size() << std::endl;

    int roi_pad = 0;
    roi_pad = get(m_cfg, "roi_pad", roi_pad);
    int raw_pad = 0;
    raw_pad = get(m_cfg, "raw_pad", raw_pad);

    double raw_ROI_th_nsigma = get(m_cfg, "raw_ROI_th_nsigma", 4);
    double raw_ROI_th_adclimit = get(m_cfg, "raw_ROI_th_adclimit", 10);

    // std::cout << "Xin: " << raw_ROI_th_nsigma << " " << raw_ROI_th_adclimit << " " << overall_time_offset << " " <<
    // collect_time_offset << " " << roi_pad << " " << adc_l1_threshold << " " << adc_sum_threshold << " " <<
    // adc_sum_rescaling << " " << adc_sum_rescaling_limit << " " << l1_seg_length << " " << l1_scaling_factor << " " <<
    // l1_lambda << " " << l1_epsilon << " " << l1_niteration << " " << l1_decon_limit << " " << l1_resp_scale << " " <<
    // l1_col_scale << " " << l1_ind_scale << std::endl;
    init_resp();

    auto adctraces = Aux::tagged_traces(in, adctag);
    auto sigtraces = Aux::tagged_traces(in, sigtag);

    if (adctraces.empty() or sigtraces.empty() or adctraces.size() != sigtraces.size()) {
        std::cerr << "L1SPFilter got unexpected input: " << adctraces.size() << " ADC traces and " << sigtraces.size()
                  << " signal traces\n";
        THROW(RuntimeError() << errmsg{"L1SPFilter: unexpected input"});
    }

    m_period = in->tick();

    // std::cout << m_period/units::us << std::endl;

    /// here, use the ADC and signal traces to do L1SP
    ///  put result in out_traces
    ITrace::vector out_traces;

    std::map<int, std::set<int>> init_map;
    // do ROI from the decon signal
    for (auto trace : sigtraces) {
        int ch = trace->channel();
        int tbin = trace->tbin();
        auto const& charges = trace->charge();
        const int ntbins = charges.size();
        std::set<int> time_ticks;

        for (int qi = 0; qi < ntbins; qi++) {
            if (charges[qi] > 0) {
                time_ticks.insert(tbin + qi);
            }
        }

        init_map[ch] = time_ticks;
        // if (time_ticks.size()>0){
        // 	std::cout << ch << " " << time_ticks.size() << std::endl;
        // }
    }

    // do ROI from the raw signal
    int ntot_ticks = 0;

    std::map<int, std::shared_ptr<const WireCell::ITrace>> adctrace_ch_map;

    for (auto trace : adctraces) {
        int ch = trace->channel();

        adctrace_ch_map[ch] = trace;

        int tbin = trace->tbin();
        auto const& charges = trace->charge();
        const int ntbins = charges.size();
        std::set<int>& time_ticks = init_map[ch];

        if (ntot_ticks < ntbins) ntot_ticks = ntbins;

        Waveform::realseq_t tmp_charge = charges;
        double mean = Waveform::percentile(tmp_charge, 0.5);
        double mean_p1sig = Waveform::percentile(tmp_charge, 0.5 + 0.34);
        double mean_n1sig = Waveform::percentile(tmp_charge, 0.5 - 0.34);
        double cut = raw_ROI_th_nsigma * sqrt((pow(mean_p1sig - mean, 2) + pow(mean_n1sig - mean, 2)) / 2.);
        if (cut < raw_ROI_th_adclimit) cut = raw_ROI_th_adclimit;

        // if (ch==4090)
        // 	std::cout << cut << " " << raw_pad << std::endl;

        for (int qi = 0; qi < ntbins; qi++) {
            if (fabs(charges[qi]) > cut) {
                for (int qii = -raw_pad; qii != raw_pad + 1; qii++) {
                    if (tbin + qi + qii >= 0 && tbin + qi + qii < ntot_ticks) time_ticks.insert(tbin + qi + qii);
                }
            }
        }
        // if (time_ticks.size()>0){
        // 	std::cout << ch << " " << time_ticks.size() << std::endl;
        // }
    }

    // create ROIs ...
    std::map<int, std::vector<std::pair<int, int>>> map_ch_rois;

    for (auto it = init_map.begin(); it != init_map.end(); it++) {
        int wire_index = it->first;
        std::set<int>& time_slices_set = it->second;
        if (time_slices_set.size() == 0) continue;
        std::vector<int> time_slices;
        std::copy(time_slices_set.begin(), time_slices_set.end(), std::back_inserter(time_slices));

        std::vector<std::pair<int, int>> rois;
        std::vector<std::pair<int, int>> rois_save;

        rois.push_back(std::make_pair(time_slices.front(), time_slices.front()));
        for (size_t i = 1; i < time_slices.size(); i++) {
            if (time_slices.at(i) - rois.back().second <= roi_pad * 2) {
                rois.back().second = time_slices.at(i);
            }
            else {
                rois.push_back(std::make_pair(time_slices.at(i), time_slices.at(i)));
            }
        }

        // extend the rois to both side according to the bin content
        for (auto it = rois.begin(); it != rois.end(); it++) {
            int start_bin = it->first;
            int end_bin = it->second;
            start_bin = start_bin - roi_pad;
            end_bin = end_bin + roi_pad;
            if (start_bin < 0) start_bin = 0;
            if (end_bin > ntot_ticks - 1) end_bin = ntot_ticks - 1;
            it->first = start_bin;
            it->second = end_bin;
        }

        for (auto it = rois.begin(); it != rois.end(); it++) {
            if (rois_save.size() == 0) {
                rois_save.push_back(*it);
            }
            else if (it->first <= rois_save.back().second) {
                rois_save.back().second = it->second;
            }
            else {
                rois_save.push_back(*it);
            }
        }

        if (rois_save.size() > 0) map_ch_rois[wire_index] = rois_save;
        // for (auto it = rois_save.begin(); it!=rois_save.end(); it++){
        // std::cout << wire_index << " " << it->first << " " << it->second +1 << std::endl;
        // }
    }

    std::map<int, std::vector<int>> map_ch_flag_rois;
    // prepare for the output signal ...
    for (auto trace : sigtraces) {
        auto newtrace = std::make_shared<Aux::SimpleTrace>(trace->channel(), trace->tbin(), trace->charge());
        // How to access the sigtraces together ???
        if (map_ch_rois.find(trace->channel()) != map_ch_rois.end()) {
            std::vector<std::pair<int, int>>& rois_save = map_ch_rois[trace->channel()];
            std::vector<int> flag_rois;
            map_ch_flag_rois[trace->channel()] = flag_rois;

            bool flag_shorted = false;
            for (size_t i1 = 0; i1 != rois_save.size(); i1++) {
                //

                int temp_flag = L1_fit(newtrace, adctrace_ch_map[trace->channel()], rois_save.at(i1).first,
                                       rois_save.at(i1).second + 1, flag_shorted);
                map_ch_flag_rois[trace->channel()].push_back(temp_flag);

                for (int time_tick = rois_save.at(i1).first; time_tick != rois_save.at(i1).second + 1; time_tick++) {
                    // temporary hack to reset the data ...
                    if (newtrace->charge().at(time_tick - trace->tbin()) < 0)
                        newtrace->charge().at(time_tick - trace->tbin()) = 0;
                }
            }

            flag_shorted = false;
            int prev_time_tick = -2000;
            for (size_t i1 = 0; i1 != rois_save.size(); i1++) {
                if (rois_save.at(i1).first - prev_time_tick > 20) flag_shorted = false;

                if (map_ch_flag_rois[trace->channel()].at(i1) == 1) {
                    flag_shorted = true;
                }
                else if (map_ch_flag_rois[trace->channel()].at(i1) == 2) {
                    if (flag_shorted) {
                        L1_fit(newtrace, adctrace_ch_map[trace->channel()], rois_save.at(i1).first,
                               rois_save.at(i1).second + 1, flag_shorted);
                        map_ch_flag_rois[trace->channel()].at(i1) = 0;
                    }
                }
                else if (map_ch_flag_rois[trace->channel()].at(i1) == 0) {
                    flag_shorted = false;
                }
                prev_time_tick = rois_save.at(i1).second;
            }

            if (rois_save.size() > 0) {
                prev_time_tick = rois_save.back().second + 2000;

                for (size_t i1 = 0; i1 != rois_save.size(); i1++) {
                    if (prev_time_tick - rois_save.at(rois_save.size() - 1 - i1).second > 20) flag_shorted = false;

                    if (map_ch_flag_rois[trace->channel()].at(i1) == 1) {
                        flag_shorted = true;
                    }
                    else if (map_ch_flag_rois[trace->channel()].at(i1) == 2) {
                        if (flag_shorted) {
                            L1_fit(newtrace, adctrace_ch_map[trace->channel()], rois_save.at(i1).first,
                                   rois_save.at(i1).second + 1, flag_shorted);
                            map_ch_flag_rois[trace->channel()].at(i1) = 0;
                        }
                    }
                    else if (map_ch_flag_rois[trace->channel()].at(i1) == 0) {
                        flag_shorted = false;
                    }
                    prev_time_tick = rois_save.at(i1).first;
                }
            }
        }

        // std::cout << trace->channel() << std::endl;
        out_traces.push_back(newtrace);
    }

    for (size_t i2 = 0; i2 != out_traces.size(); i2++) {
        auto new_trace = std::make_shared<Aux::SimpleTrace>(out_traces.at(i2)->channel(), out_traces.at(i2)->tbin(),
                                                       out_traces.at(i2)->charge());
        int ch = out_traces.at(i2)->channel();
        if (map_ch_rois.find(ch) != map_ch_rois.end()) {
            std::vector<std::pair<int, int>>& rois_save = map_ch_rois[ch];
            std::vector<std::pair<int, int>> next_rois_save;
            if (map_ch_rois.find(ch + 1) != map_ch_rois.end()) next_rois_save = map_ch_rois[ch + 1];
            std::vector<std::pair<int, int>> prev_rois_save;
            if (map_ch_rois.find(ch - 1) != map_ch_rois.end()) prev_rois_save = map_ch_rois[ch - 1];

            for (size_t i1 = 0; i1 != rois_save.size(); i1++) {
                if (map_ch_flag_rois[ch].at(i1) == 2) {
                    bool flag_shorted = false;
                    // need to check now ...
                    for (size_t i3 = 0; i3 != next_rois_save.size(); i3++) {
                        if (map_ch_flag_rois[ch + 1].at(i3) == 1) {
                            if (rois_save.at(i1).first - 3 <= next_rois_save.at(i3).second + 3 &&
                                rois_save.at(i1).second + 3 >= next_rois_save.at(i3).first - 3) {
                                flag_shorted = true;
                                break;
                            }
                        }
                    }
                    if (!flag_shorted)
                        for (size_t i3 = 0; i3 != prev_rois_save.size(); i3++) {
                            if (map_ch_flag_rois[ch - 1].at(i3) == 1) {
                                if (rois_save.at(i1).first - 3 <= prev_rois_save.at(i3).second + 3 &&
                                    rois_save.at(i1).second + 3 >= prev_rois_save.at(i3).first - 3) {
                                    flag_shorted = true;
                                    break;
                                }
                            }
                        }
                    if (flag_shorted) {
                        for (int time_tick = rois_save.at(i1).first; time_tick != rois_save.at(i1).second + 1;
                             time_tick++) {
                            new_trace->charge().at(time_tick - out_traces.at(i2)->tbin()) = 0;
                        }
                    }
                }
            }
        }
        out_traces.at(i2) = new_trace;
    }

    std::cerr << "L1SPFilter: frame: " << in->ident() << " " << adctag << "[" << adctraces.size() << "] + " << sigtag
              << "[" << sigtraces.size() << "] --> " << outtag << "[" << out_traces.size() << "]\n";

    // Finally, we save the traces to an output frame with tags.

    IFrame::trace_list_t tl(out_traces.size());
    std::iota(tl.begin(), tl.end(), 0);

    auto sf = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick());
    sf->tag_traces(outtag, tl);
    out = IFrame::pointer(sf);

    return true;
}

int L1SPFilter::L1_fit(std::shared_ptr<Aux::SimpleTrace>& newtrace,
                       std::shared_ptr<const WireCell::ITrace>& adctrace, int start_tick, int end_tick,
                       bool flag_shorted)
{
    // fixme: this use of units is broken!!!
    double overall_time_offset = get(m_cfg, "overall_time_offset", 0.0) * units::us;
    double collect_time_offset = get(m_cfg, "collect_time_offset", 3.0) * units::us;

    double adc_l1_threshold = get(m_cfg, "adc_l1_threshold", 6);
    double adc_sum_threshold = get(m_cfg, "adc_sum_threshold", 160);
    double adc_sum_rescaling = get(m_cfg, "adc_sum_rescaling", 90.0);
    double adc_sum_rescaling_limit = get(m_cfg, "adc_sum_rescaling_limit", 50.0);
    double adc_ratio_threshold = get(m_cfg, "adc_ratio_threshold", 0.2);

    double l1_seg_length = get(m_cfg, "l1_seg_length", 120);
    double l1_scaling_factor = get(m_cfg, "l1_scaling_factor", 500);
    double l1_lambda = get(m_cfg, "l1_lambda", 5);
    double l1_epsilon = get(m_cfg, "l1_epsilon", 0.05);
    double l1_niteration = get(m_cfg, "l1_niteration", 100000);
    double l1_decon_limit = get(m_cfg, "l1_decon_limit", 100);

    double l1_resp_scale = get(m_cfg, "l1_resp_scale", 0.5);
    double l1_col_scale = get(m_cfg, "l1_col_scale", 1.15);
    double l1_ind_scale = get(m_cfg, "l1_ind_scale", 0.5);

    double peak_threshold = get(m_cfg, "peak_threshold", 10000);
    double mean_threshold = get(m_cfg, "mean_threshold", 500);

    std::vector<double> smearing_vec = get<std::vector<double>>(m_cfg, "filter");

    // algorithm
    const int nbin_fit = end_tick - start_tick;

    // fill the data ...
    VectorXd init_W = VectorXd::Zero(nbin_fit);
    VectorXd init_beta = VectorXd::Zero(nbin_fit);
    VectorXd final_beta = VectorXd::Zero(nbin_fit * 2);

    double temp_sum = 0;
    double temp1_sum = 0;
    double temp2_sum = 0;
    double max_val = -1;
    double min_val = 1;
    for (int i = 0; i != nbin_fit; i++) {
        init_W(i) = adctrace->charge().at(i + start_tick - newtrace->tbin());
        init_beta(i) = newtrace->charge().at(i + start_tick - newtrace->tbin());

        if (max_val < init_W(i)) max_val = init_W(i);
        if (min_val > init_W(i)) min_val = init_W(i);

        if (fabs(init_W(i)) > adc_l1_threshold) {
            temp_sum += init_W(i);
            temp1_sum += fabs(init_W(i));
            temp2_sum += fabs(init_beta(i));
        }
    }

    int flag_l1 = 0;  // do nothing
    // 1 do L1
    if (temp_sum / (temp1_sum * adc_sum_rescaling * 1.0 / nbin_fit) > adc_ratio_threshold &&
        temp1_sum > adc_sum_threshold) {
        flag_l1 = 1;  // do L1 ...
    }
    else if (temp1_sum * adc_sum_rescaling * 1.0 / nbin_fit < adc_sum_rescaling_limit) {
        flag_l1 = 2;  // remove signal ...
    }
    else if (temp2_sum > 30 * nbin_fit && temp1_sum < 2.0 * nbin_fit && max_val - min_val < 22) {
        flag_l1 = 2;
    }

    // if (adctrace->channel() == 4079){
    // std::cout << adctrace->channel() << " " << nbin_fit << " " << start_tick << " " << end_tick << " " << temp_sum <<
    // " " << temp1_sum << " " << temp2_sum << " " << max_val << " " << min_val << " " << flag_l1 << std::endl;
    // }

    // std::cout << temp_sum << " " << temp1_sum << " " << temp_sum/(temp1_sum*adc_sum_rescaling*1.0/nbin_fit) << " " <<
    // adc_ratio_threshold << " " << temp1_sum*adc_sum_rescaling*1.0/nbin_fit << " " << flag_l1 << std::endl;

    if (flag_l1 == 1) {
        // do L1 fit ...
        int n_section = std::round(nbin_fit / l1_seg_length * 1.0);
        if (n_section == 0) n_section = 1;
        std::vector<int> boundaries;
        for (int i = 0; i != n_section; i++) {
            boundaries.push_back(int(i * nbin_fit / n_section));
        }
        boundaries.push_back(nbin_fit);

        for (int i = 0; i != n_section; i++) {
            int temp_nbin_fit = boundaries.at(i + 1) - boundaries.at(i);
            VectorXd W = VectorXd::Zero(temp_nbin_fit);
            for (int j = 0; j != temp_nbin_fit; j++) {
                W(j) = init_W(boundaries.at(i) + j);
            }

            // for matrix G
            MatrixXd G = MatrixXd::Zero(temp_nbin_fit, temp_nbin_fit * 2);

            for (int i = 0; i != temp_nbin_fit; i++) {
                // X
                double t1 = i / 2. * units::us;  // us, measured time
                for (int j = 0; j != temp_nbin_fit; j++) {
                    // Y ...
                    double t2 = j / 2. * units::us;  // us, real signal time
                    double delta_t = t1 - t2;
                    if (delta_t > -15 * units::us - overall_time_offset &&
                        delta_t < 10 * units::us - overall_time_offset) {
                        G(i, j) = (*lin_W)(delta_t + overall_time_offset + collect_time_offset) * l1_scaling_factor *
                                  l1_resp_scale;
                        G(i, temp_nbin_fit + j) =
                            (*lin_V)(delta_t + overall_time_offset) * l1_scaling_factor * l1_resp_scale;
                    }
                }
            }

            double lambda = l1_lambda;  // 1/2.;
            WireCell::LassoModel m2(lambda, l1_niteration, l1_epsilon);
            m2.SetData(G, W);
            m2.Fit();
            VectorXd beta = m2.Getbeta();
            for (int j = 0; j != temp_nbin_fit; j++) {
                final_beta(boundaries.at(i) + j) = beta(j);
                final_beta(nbin_fit + boundaries.at(i) + j) = beta(temp_nbin_fit + j);
            }
        }

        double sum1 = 0;
        // double sum2 = 0;
        for (int i = 0; i != nbin_fit; i++) {
            sum1 += final_beta(i);
            // sum2 += final_beta(nbin_fit + i);
        }

        if (sum1 > adc_l1_threshold) {
            Waveform::realseq_t l1_signal(nbin_fit, 0);
            Waveform::realseq_t l2_signal(nbin_fit, 0);
            for (int j = 0; j != nbin_fit; j++) {
                l1_signal.at(j) = final_beta(j) * l1_col_scale + final_beta(nbin_fit + j) * l1_ind_scale;
            }
            int mid_bin = (smearing_vec.size() - 1) / 2;
            // std::cout << smearing_vec.size() << " " << mid_bin << std::endl;
            for (int j = 0; j != nbin_fit; j++) {
                double content = l1_signal.at(j);
                if (content > 0) {
                    for (size_t k = 0; k != smearing_vec.size(); k++) {
                        int bin = j + k - mid_bin;
                        if (bin >= 0 && bin < nbin_fit) l2_signal.at(bin) += content * smearing_vec.at(k);
                    }
                }
            }

            for (int j = 0; j != nbin_fit; j++) {
                if (l2_signal.at(j) < l1_decon_limit / l1_scaling_factor) {  // 50 electrons
                    l1_signal.at(j) = 0;
                }
                else {
                    l1_signal.at(j) = l2_signal.at(j) * l1_scaling_factor;
                }
            }

            {
                // go through the new data and clean the small peaks ...
                std::vector<int> nonzero_bins;
                std::vector<std::pair<int, int>> ROIs;
                for (int j = 0; j != nbin_fit; j++) {
                    if (l1_signal.at(j) > 0) nonzero_bins.push_back(j);
                }

                bool flag_start = false;
                int start_bin = -1, end_bin = -1;
                for (size_t j = 0; j != nonzero_bins.size(); j++) {
                    if (!flag_start) {
                        start_bin = nonzero_bins.at(j);
                        end_bin = nonzero_bins.at(j);
                        flag_start = true;
                    }
                    else {
                        if (nonzero_bins.at(j) - end_bin == 1) {
                            end_bin = nonzero_bins.at(j);
                        }
                        else {
                            ROIs.push_back(std::make_pair(start_bin, end_bin));
                            start_bin = nonzero_bins.at(j);
                            end_bin = nonzero_bins.at(j);
                        }
                    }
                }
                if (start_bin >= 0) ROIs.push_back(std::make_pair(start_bin, end_bin));

                for (size_t j = 0; j != ROIs.size(); j++) {
                    double max_val = -1;
                    double mean_val1 = 0;
                    double mean_val2 = 0;
                    for (int k = ROIs.at(j).first; k <= ROIs.at(j).second; k++) {
                        if (l1_signal.at(k) > max_val) max_val = l1_signal.at(k);
                        mean_val1 += l1_signal.at(k);
                        mean_val2++;
                    }
                    if (mean_val2 != 0) mean_val1 /= mean_val2;
                    if (max_val < peak_threshold && mean_val1 < mean_threshold) {
                        for (int k = ROIs.at(j).first; k <= ROIs.at(j).second; k++) {
                            l1_signal.at(k) = 0;
                        }
                    }
                    // std::cout << max_val << " " << mean_val1 << std::endl;
                }

                // std::cout << nonzero_bins.front() << " X " << nonzero_bins.back() << std::endl;
                // std::cout << ROIs.size() << std::endl;
                // for (size_t i=0;i!=ROIs.size();i++){
                //   std::cout << ROIs.at(i).first << " " << ROIs.at(i).second << std::endl;
                // }

                // finish cleaning ...
            }

            for (int time_tick = start_tick; time_tick != end_tick; time_tick++) {
                newtrace->charge().at(time_tick - newtrace->tbin()) = l1_signal.at(time_tick - start_tick);
            }
        }
    }
    else if (flag_l1 == 2) {
        if (flag_shorted) {
            for (int time_tick = start_tick; time_tick != end_tick; time_tick++) {
                // temporary hack to reset the data ...
                newtrace->charge().at(time_tick - newtrace->tbin()) = 0;
            }
        }
    }

    return flag_l1;
}
