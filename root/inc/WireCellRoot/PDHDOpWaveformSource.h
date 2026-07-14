/** Read PDHD photon-detector waveform snippets from the temporary
 * light-data ROOT files (decoana TH1Ds) and emit one IFrame per
 * configured event.
 *
 * Each TH1D snippet (1024 ticks @ 16 ns) becomes one trace with
 * channel = offline OpChannel (== OpDet for PDHD).  Raw (ADC) and
 * LArSoft-deconvolved (SPE-normalized) snippets are emitted as
 * "raw" and "deconv" tagged trace sets.
 *
 * Time convention (see flash/docs/design.md §3.1): the decoana
 * histogram x-axis encodes the snippet start in 16 ns ticks relative
 * to the event's earliest snippet (t_first), whose absolute DTS time
 * is not recorded.  t_first is recovered as the mode of
 * (OpHit PeakTimeAbs - snippet start - deconv peak bin) over clear
 * deconv peaks.  The frame time() is (t_first - trigger) in WCT ns
 * and trace tbin counts ticks from t_first.
 */
#ifndef WIRECELLROOT_PDHDOPWAVEFORMSOURCE
#define WIRECELLROOT_PDHDOPWAVEFORMSOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Root {
        class PDHDOpWaveformSource : public Aux::Logger, public IFrameSource, public IConfigurable {
          public:
            PDHDOpWaveformSource();
            virtual ~PDHDOpWaveformSource();

            virtual bool operator()(IFrame::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_filename{""};
            int m_run{-1};
            int m_subrun{-1};  // <0: any
            int m_event{-1};
            double m_nominal_offset_us{250.0};
            std::string m_frame_tag{"pdhd_light"};
            // t_first recovery: only deconv snippets whose peak exceeds
            // this (PE/tick) vote, and only votes within the search
            // window around rd_timestamp are accepted.  Snippets are
            // observed up to ~3 ms before the TPC readout start, hence
            // the wide default.
            double m_min_peak{2.0};
            double m_t0_search_window_us{5000.0};

            int m_calls{0};
        };
    }  // namespace Root
}  // namespace WireCell

#endif
