/** Read PDVD photon-detector waveforms from the flat raw_waveform
 * TTree ROOT files and emit one IFrame per configured event.
 *
 * One TTree entry = one waveform: 1024-tick self-trigger snippets
 * (membrane X-ARAPUCA opch 20xx, PMT 30xx) or ~468k-tick cathode
 * full streams (opch 10xx), all 16 ns ticks, positive polarity.
 * Each becomes one trace with channel = offline DAPHNE OpChannel
 * (NOT the OpDet: cathode/membrane XAs have two channels per OpDet;
 * the ganging to OpDets happens in OpFlashFinder via its channel
 * map).
 *
 * Time convention: the `timestamp` branch is MICROSECONDS on one
 * common 16 ns-quantized clock for all record types.  A full stream
 * starts at its timestamp; a snippet's sample `trig_sample` is at
 * its timestamp.  The tick origin t0 is the earliest record start
 * of the event over ALL channels (independent of the opch selection,
 * so parallel per-population instances share one origin); trace tbin
 * counts ticks from t0 and frame time() is 0.  The light<->charge
 * offset is carried later by OpFlashFinder's offset_us metadata.
 *
 * Config `opch_lo`/`opch_hi` (inclusive, -1 = no bound) select the
 * population: cathode 1000-1999, membrane 2000-2999, PMT 3000-3999.
 */
#ifndef WIRECELLROOT_PDVDOPWAVEFORMSOURCE
#define WIRECELLROOT_PDVDOPWAVEFORMSOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Root {
        class PDVDOpWaveformSource : public Aux::Logger, public IFrameSource, public IConfigurable {
          public:
            PDVDOpWaveformSource();
            virtual ~PDVDOpWaveformSource();

            virtual bool operator()(IFrame::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_filename{""};
            int m_run{-1};
            int m_subrun{-1};  // <0: any
            int m_event{-1};
            std::string m_frame_tag{"pdvd_light"};
            // Inclusive OpChannel selection; -1 = unbounded.
            int m_opch_lo{-1};
            int m_opch_hi{-1};
            // A snippet's timestamp marks this sample (DAPHNE
            // self-trigger pre-trigger length); records with
            // nsamp <= snippet_nsamp are snippets, longer records are
            // full streams whose timestamp marks sample 0.
            int m_trig_sample{64};
            int m_snippet_nsamp{1024};

            int m_calls{0};
        };
    }  // namespace Root
}  // namespace WireCell

#endif
