#ifndef WIRECELLSIGPROC_OMNIBUSPMTNOISEFILTER
#define WIRECELLSIGPROC_OMNIBUSPMTNOISEFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Waveform.h"

#include <vector>

class PMTNoiseROI;  // forward declaration

namespace WireCell {
    namespace SigProc {

        class OmnibusPMTNoiseFilter : public WireCell::IFrameFilter, public WireCell::IConfigurable {
           public:
            typedef std::vector<std::vector<int> > grouped_channels_t;

            /// Create an OmnibusPMTNoiseFilter.
            OmnibusPMTNoiseFilter(const std::string anode_tn = "AnodePlane", int pad_window = 5,
                                  int min_window_length = 4, int threshold = 4, float rms_threshold = 0.5,
                                  int sort_wires = 4, float ind_th1 = 2.0, float ind_th2 = 0.5,
                                  int nwire_pmt_col_th = 6);
            virtual ~OmnibusPMTNoiseFilter();

            /// IFrameFilter interface.
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            /// IConfigurable interface.
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            /// Explicitly inject required services
            void IDPMTSignalCollection(Waveform::realseq_t& signal, double rms, int ch);
            void IDPMTSignalInduction(Waveform::realseq_t& signal, double rms, int ch, int plane);
            void RemovePMTSignal(Waveform::realseq_t& signal, int start_bin, int end_bin, int flag = 0);

           private:
            std::string m_intag, m_outtag;
            std::string m_anode_tn;
            IAnodePlane::pointer m_anode;

            int m_pad_window;
            int m_min_window_length;
            int m_threshold;
            float m_rms_threshold;

            int m_sort_wires;
            float m_ind_th1;
            float m_ind_th2;
            int m_nwire_pmt_col_th;

            std::vector<PMTNoiseROI*> m_pmt_rois;
        };
    }  // namespace SigProc
}  // namespace WireCell

#endif
// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
