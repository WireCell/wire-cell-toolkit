#include "WireCellSigProc/OmnibusPMTNoiseFilter.h"

#include "WireCellSigProc/Derivations.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/FrameTools.h"

WIRECELL_FACTORY(OmnibusPMTNoiseFilter, WireCell::SigProc::OmnibusPMTNoiseFilter, WireCell::IFrameFilter,
                 WireCell::IConfigurable)

#include "PMTNoiseROI.h"

using namespace WireCell;

using namespace WireCell::SigProc;

// m_pmt_rois moved from file-scope global to member variable m_pmt_rois
// to ensure thread safety and prevent state leakage across instances.

OmnibusPMTNoiseFilter::OmnibusPMTNoiseFilter(const std::string anode_tn, int pad_window, int min_window_length,
                                             int threshold, float rms_threshold, int sort_wires, float ind_th1,
                                             float ind_th2, int nwire_pmt_col_th)
  : m_intag("quiet")
  , m_outtag("raw")
  , m_anode_tn(anode_tn)
  , m_pad_window(pad_window)
  , m_min_window_length(min_window_length)
  , m_threshold(threshold)
  , m_rms_threshold(rms_threshold)
  , m_sort_wires(sort_wires)
  , m_ind_th1(ind_th1)
  , m_ind_th2(ind_th2)
  , m_nwire_pmt_col_th(nwire_pmt_col_th)
{
}
OmnibusPMTNoiseFilter::~OmnibusPMTNoiseFilter() {}

void OmnibusPMTNoiseFilter::configure(const WireCell::Configuration& config)
{
    m_anode_tn = get(config, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
    }

    m_pad_window = get(config, "pad_window", m_pad_window);
    m_min_window_length = get(config, "min_window_length", m_min_window_length);
    m_threshold = get(config, "threshold", m_threshold);
    m_rms_threshold = get(config, "rms_threshold", m_rms_threshold);
    m_sort_wires = get(config, "sort_wires", m_sort_wires);
    m_ind_th1 = get(config, "ind_th1", m_ind_th1);
    m_ind_th2 = get(config, "ind_th2", m_ind_th2);
    m_nwire_pmt_col_th = get(config, "nwire_pmt_col_th", m_nwire_pmt_col_th);
    m_intag = get(config, "intraces", m_intag);
    m_outtag = get(config, "outtraces", m_outtag);
}
WireCell::Configuration OmnibusPMTNoiseFilter::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["pad_window"] = m_pad_window;
    cfg["min_window_length"] = m_min_window_length;
    cfg["threshold"] = m_threshold;
    cfg["rms_threshold"] = m_rms_threshold;
    cfg["sort_wires"] = m_sort_wires;
    cfg["ind_th1"] = m_ind_th1;
    cfg["ind_th2"] = m_ind_th2;
    cfg["nwire_pmt_col_th"] = m_nwire_pmt_col_th;
    cfg["intraces"] = m_intag;
    cfg["outtraces"] = m_outtag;
    return cfg;
}

bool OmnibusPMTNoiseFilter::operator()(const input_pointer& in, output_pointer& out)
{
    if (!in) {  // eos processing
        out = nullptr;
        return true;
    }

    std::map<int, Waveform::realseq_t> bychan_coll;
    std::map<int, Waveform::realseq_t> bychan_indu;
    std::map<int, Waveform::realseq_t> bychan_indv;
    std::map<int, double> by_chan_rms;
    // go through all channels and calculate RMS as well as categorize them
    auto traces = Aux::tagged_traces(in, m_intag);
    for (const auto& trace : traces) {
        int ch = trace->channel();

        auto wpid = m_anode->resolve(ch);
        const int iplane = wpid.index();

        // Store waveform once, compute RMS from the stored copy
        if (iplane == 0) {
            bychan_indu[ch] = trace->charge();
            std::pair<double, double> results = Derivations::CalcRMS(bychan_indu[ch]);
            by_chan_rms[ch] = results.second;
        }
        else if (iplane == 1) {
            bychan_indv[ch] = trace->charge();
            std::pair<double, double> results = Derivations::CalcRMS(bychan_indv[ch]);
            by_chan_rms[ch] = results.second;
        }
        else {
            bychan_coll[ch] = trace->charge();
            std::pair<double, double> results = Derivations::CalcRMS(bychan_coll[ch]);
            by_chan_rms[ch] = results.second;
        }
    }

    // Remove PMT signal from Collection
    for (auto& cs : bychan_coll) {
        // ignore dead channels ...
        if (by_chan_rms[cs.first] > m_rms_threshold)
            IDPMTSignalCollection(cs.second, by_chan_rms[cs.first], cs.first);
    }

    // ID PMT signal in induction signal
    for (auto& cs : bychan_indu) {
        // ignore dead channels ...
        if (by_chan_rms[cs.first] > m_rms_threshold)
            IDPMTSignalInduction(cs.second, by_chan_rms[cs.first], cs.first, 0);
    }
    for (auto& cs : bychan_indv) {
        // ignore dead channels ...
        if (by_chan_rms[cs.first] > m_rms_threshold)
            IDPMTSignalInduction(cs.second, by_chan_rms[cs.first], cs.first, 1);
    }

    // Remove PMT signal from Induction ...

    for (int i = 0; i != int(m_pmt_rois.size()); i++) {
        m_pmt_rois.at(i)->sort_wires(m_sort_wires);
        // int flag_qx = 0;
        if (m_pmt_rois.at(i)->get_sorted_uwires().size() > 0 && m_pmt_rois.at(i)->get_sorted_vwires().size() > 0) {
            for (int j = 0; j != int(m_pmt_rois.at(i)->get_sorted_uwires().size()); j++) {
                if (m_pmt_rois.at(i)->get_average_uwires_peak_height(j) <
                        m_ind_th1 * m_pmt_rois.at(i)->get_average_wwires_peak_height() &&
                    m_pmt_rois.at(i)->get_max_uwires_peak_height(j) <
                        m_ind_th2 * m_pmt_rois.at(i)->get_max_wwires_peak_height() &&
                    m_pmt_rois.at(i)->get_sorted_uwires().at(j).size() <= m_pmt_rois.at(i)->get_sorted_wwires().size()) {
                    //	  flag_qx = 1;
                    for (int k = 0; k != int(m_pmt_rois.at(i)->get_sorted_uwires().at(j).size()); k++) {
                        RemovePMTSignal(bychan_indu[m_pmt_rois.at(i)->get_sorted_uwires().at(j).at(k)],
                                        m_pmt_rois.at(i)->get_start_bin(), m_pmt_rois.at(i)->get_end_bin());
                        // RemovePMTSignalInduction(hu[m_pmt_rois.at(i)->get_sorted_uwires().at(j).at(k)],m_pmt_rois.at(i)->get_start_bin(),m_pmt_rois.at(i)->get_end_bin());
                    }

                    // std::cout << i << " " <<  m_pmt_rois.at(i)->get_peaks().at(0) << " " <<
                    // m_pmt_rois.at(i)->get_peaks().size() << " " << m_pmt_rois.at(i)->get_sorted_uwires().at(j).size() <<
                    // " " << j << " " << m_pmt_rois.at(i)->get_average_uwires_peak_height(j) << " " <<
                    // m_pmt_rois.at(i)->get_average_wwires_peak_height() << " " <<
                    // m_pmt_rois.at(i)->get_max_uwires_peak_height(j) << " " <<
                    // m_pmt_rois.at(i)->get_max_wwires_peak_height() << " " <<
                    // m_pmt_rois.at(i)->get_sorted_uwires().at(j).size() << " " <<
                    // m_pmt_rois.at(i)->get_sorted_wwires().size() << std::endl;
                }
            }

            for (int j = 0; j != int(m_pmt_rois.at(i)->get_sorted_vwires().size()); j++) {
                if (m_pmt_rois.at(i)->get_average_vwires_peak_height(j) <
                        m_ind_th1 * m_pmt_rois.at(i)->get_average_wwires_peak_height() &&
                    m_pmt_rois.at(i)->get_max_vwires_peak_height(j) <
                        m_ind_th2 * m_pmt_rois.at(i)->get_max_wwires_peak_height() &&
                    m_pmt_rois.at(i)->get_sorted_vwires().at(j).size() <= m_pmt_rois.at(i)->get_sorted_wwires().size()) {
                    // flag_qx = 1;
                    for (int k = 0; k != int(m_pmt_rois.at(i)->get_sorted_vwires().at(j).size()); k++) {
                        RemovePMTSignal(bychan_indv[m_pmt_rois.at(i)->get_sorted_vwires().at(j).at(k)],
                                        m_pmt_rois.at(i)->get_start_bin(), m_pmt_rois.at(i)->get_end_bin());
                        // RemovePMTSignalInduction(hv[m_pmt_rois.at(i)->get_sorted_vwires().at(j).at(k)],m_pmt_rois.at(i)->get_start_bin(),m_pmt_rois.at(i)->get_end_bin());
                    }
                    //  std::cout << i << " " <<  m_pmt_rois.at(i)->get_peaks().at(0) << " " <<
                    //  m_pmt_rois.at(i)->get_peaks().size() << " " << m_pmt_rois.at(i)->get_sorted_vwires().at(j).size() <<
                    //  " " << j << " " << m_pmt_rois.at(i)->get_average_vwires_peak_height(j) << " " <<
                    //  m_pmt_rois.at(i)->get_average_wwires_peak_height() << " " <<
                    //  m_pmt_rois.at(i)->get_max_vwires_peak_height(j) << " " <<
                    //  m_pmt_rois.at(i)->get_max_wwires_peak_height() << " " <<
                    //  m_pmt_rois.at(i)->get_sorted_vwires().at(j).size() << " " <<
                    //  m_pmt_rois.at(i)->get_sorted_wwires().size() << std::endl;
                }
            }

            if (int(m_pmt_rois.at(i)->get_sorted_wwires().size()) >= m_nwire_pmt_col_th) {
                for (int j = 0; j != int(m_pmt_rois.at(i)->get_sorted_wwires().size()); j++) {
                    RemovePMTSignal(bychan_coll[m_pmt_rois.at(i)->get_sorted_wwires().at(j)],
                                    m_pmt_rois.at(i)->get_start_bin(), m_pmt_rois.at(i)->get_end_bin(), 1);
                }
            }

            // if (flag_qx == 1) std::cout << i << " " <<  m_pmt_rois.at(i)->get_peaks().at(0) << " " <<
            // m_pmt_rois.at(i)->get_peaks().size() << " " << m_pmt_rois.at(i)->get_uwires().size() << " " <<
            // m_pmt_rois.at(i)->get_vwires().size() << " " << m_pmt_rois.at(i)->get_sorted_uwires().size() << " " <<
            // m_pmt_rois.at(i)->get_sorted_vwires().size() << " " << m_pmt_rois.at(i)->get_sorted_wwires().size() << " " <<
            // m_pmt_rois.at(i)->get_average_uwires_peak_height() << " " <<
            // m_pmt_rois.at(i)->get_average_vwires_peak_height() << " " <<
            // m_pmt_rois.at(i)->get_average_wwires_peak_height() << " " << m_pmt_rois.at(i)->get_max_uwires_peak_height() <<
            // " " << m_pmt_rois.at(i)->get_max_vwires_peak_height() << " " << m_pmt_rois.at(i)->get_max_wwires_peak_height()
            // << std::endl;
        }
        delete m_pmt_rois.at(i);
    }
    m_pmt_rois.clear();

    // load results ...

    ITrace::vector itraces;
    for (auto cs : bychan_indu) {
        // fixme: that tbin though
        auto trace = new Aux::SimpleTrace(cs.first, 0, cs.second);
        itraces.push_back(ITrace::pointer(trace));
    }
    for (auto cs : bychan_indv) {
        // fixme: that tbin though
        auto trace = new Aux::SimpleTrace(cs.first, 0, cs.second);
        itraces.push_back(ITrace::pointer(trace));
    }
    for (auto cs : bychan_coll) {
        // fixme: that tbin though
        auto trace = new Aux::SimpleTrace(cs.first, 0, cs.second);
        itraces.push_back(ITrace::pointer(trace));
    }

    auto sframe = new Aux::SimpleFrame(in->ident(), in->time(), itraces, in->tick(), in->masks());
    if (! m_outtag.empty()) {
        IFrame::trace_list_t indices(itraces.size());
        for (size_t ind = 0; ind < itraces.size(); ++ind) {
            indices[ind] = ind;
        }
        sframe->tag_traces(m_outtag, indices);
    }

    out = IFrame::pointer(sframe);

    return true;
}

void OmnibusPMTNoiseFilter::IDPMTSignalCollection(Waveform::realseq_t& signal, double rms, int ch)
{
    int flag_start = 0;
    int start_bin = 0;
    int end_bin = 0;
    int peak_bin = 0;

    //  std::cout << m_threshold << " " << m_pad_window << " " << m_min_window_length << " " << rms << " " << ch <<
    //  std::endl;

    for (int i = 0; i != int(signal.size()); i++) {
        float content = signal.at(i);

        if (flag_start == 0) {
            if (content < -m_threshold * rms) {
                start_bin = i;
                flag_start = 1;
            }
        }
        else {
            if (content >= -m_threshold * rms) {
                end_bin = i - 1;
                if (end_bin > start_bin + m_min_window_length) {
                    float min = signal.at(start_bin);
                    peak_bin = start_bin;
                    for (int j = start_bin + 1; j != end_bin; j++) {
                        if (signal.at(j) < min) {
                            min = signal.at(j);
                            peak_bin = j;
                        }
                    }

                    PMTNoiseROI* ROI = new PMTNoiseROI(start_bin, end_bin, peak_bin, ch, signal.at(peak_bin));

                    // std::cout << start_bin << " " << end_bin << std::endl;
                    if (m_pmt_rois.size() == 0) {
                        m_pmt_rois.push_back(ROI);
                    }
                    else {
                        bool flag_merge = false;
                        for (int i = 0; i != int(m_pmt_rois.size()); i++) {
                            flag_merge = m_pmt_rois.at(i)->merge_ROI(*ROI);
                            if (flag_merge) {
                                delete ROI;
                                break;
                            }
                        }
                        if (!flag_merge) {
                            m_pmt_rois.push_back(ROI);
                        }
                    }
                    // std::cout << "h " << m_pmt_rois.size() << std::endl;

                    flag_start = 0;

                    // //adaptive baseline
                    //  float start_content =signal.at(start_bin);
                    //  for (int j=start_bin;j>=start_bin - m_pad_window;j--){
                    //    if (j<0) continue;
                    //    if (fabs(signal.at(j)) < fabs(start_content)){
                    // 	start_bin = j;
                    // 	start_content = signal.at(start_bin);
                    //    }
                    //  }
                    //  float end_content = signal.at(end_bin);
                    //  for (int j=end_bin; j<=end_bin+m_pad_window;j++){
                    //    if (j>=int(signal.size())) continue;
                    //    if (fabs(signal.at(j)) < fabs(end_content)){
                    // 	end_bin = j;
                    // 	end_content = signal.at(end_bin);
                    //    }
                    //  }

                    //  for (int j=start_bin;j<=end_bin;j++){
                    //    float content = start_content + (end_content - start_content) * (j - start_bin) / (end_bin -
                    //    start_bin*1.0); signal.at(j)=content;
                    //  }
                }
            }
        }
    }
}

void OmnibusPMTNoiseFilter::IDPMTSignalInduction(Waveform::realseq_t& signal, double rms, int ch, int plane)
{
    // std::cout <<  m_pmt_rois.size() << std::endl;

    for (int i = 0; i != int(m_pmt_rois.size()); i++) {
        PMTNoiseROI* ROI = m_pmt_rois.at(i);
        for (int j = 0; j != int(ROI->get_peaks().size()); j++) {
            int peak = ROI->get_peaks().at(j);
            int peak_m1 = peak - 1;
            if (peak_m1 < 0) peak_m1 = 0;
            int peak_m2 = peak - 2;
            if (peak_m2 < 0) peak_m2 = 0;
            int peak_m3 = peak - 3;
            if (peak_m3 < 0) peak_m3 = 0;
            int peak_p1 = peak + 1;
            if (peak_p1 >= int(signal.size())) peak_p1 = int(signal.size()) - 1;
            int peak_p2 = peak + 2;
            if (peak_p2 >= int(signal.size())) peak_p2 = int(signal.size()) - 1;
            int peak_p3 = peak + 3;
            if (peak_p3 >= int(signal.size())) peak_p3 = int(signal.size()) - 1;

            if (fabs(signal.at(peak)) > m_threshold * rms &&
                fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) >
                    fabs(signal.at(peak_m1)) + fabs(signal.at(peak_m2)) + fabs(signal.at(peak)) &&
                fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) >
                    fabs(signal.at(peak_p1)) + fabs(signal.at(peak_p2)) + fabs(signal.at(peak)) &&
                fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) >
                    fabs(signal.at(peak_m2)) + fabs(signal.at(peak_m3)) + fabs(signal.at(peak_m1)) &&
                fabs(signal.at(peak)) + fabs(signal.at(peak_m1)) + fabs(signal.at(peak_p1)) >
                    fabs(signal.at(peak_p2)) + fabs(signal.at(peak_p3)) + fabs(signal.at(peak_p1))) {
                //	    fabs(signal.at(peak))>= fabs(signal.at(peak_m2)) &&
                // fabs(signal.at(peak))>= fabs(signal.at(peak_p1)) &&
                // fabs(signal.at(peak))>= fabs(signal.at(peak_p2)) &&
                // fabs(signal.at(peak))>= fabs(signal.at(peak_p3)) &&
                // fabs(signal.at(peak))>= fabs(signal.at(peak_m3))
                // ){
                if (plane == 0) {
                    ROI->insert_uwires(ch, fabs(signal.at(peak)));
                    break;
                }
                else {
                    ROI->insert_vwires(ch, fabs(signal.at(peak)));
                    break;
                }
            }
        }
    }
}

void OmnibusPMTNoiseFilter::RemovePMTSignal(Waveform::realseq_t& signal, int start_bin, int end_bin, int flag)
{
    // int flag_start = 0;

    // adaptive baseline
    float start_content = signal.at(start_bin);
    for (int j = start_bin; j >= start_bin - m_pad_window; j--) {
        if (j < 0) continue;
        if (fabs(signal.at(j)) < fabs(start_content)) {
            start_bin = j;
            start_content = signal.at(start_bin);
        }
    }
    float end_content = signal.at(end_bin);
    for (int j = end_bin; j <= end_bin + m_pad_window; j++) {
        if (j >= int(signal.size())) continue;
        if (fabs(signal.at(j)) < fabs(end_content)) {
            end_bin = j;
            end_content = signal.at(end_bin);
        }
    }

    for (int j = start_bin; j <= end_bin; j++) {
        float content = (end_bin != start_bin)
            ? start_content + (end_content - start_content) * (j - start_bin) / (end_bin - start_bin * 1.0)
            : start_content;
        if (flag == 1) {
            if (signal.at(j) < 0) signal.at(j) = content;
        }
        else {
            signal.at(j) = content;
        }
    }
}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
