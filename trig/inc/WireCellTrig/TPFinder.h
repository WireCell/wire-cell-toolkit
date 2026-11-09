/** Find trigger primitives in frames of ADC waveforms.

    Each trace is run through a per-channel TPG chain (TPGChannel) configured
    per wire plane.  One TP set is output per input frame, covering the
    frame's sample span.

    Configuration:

    - period :: the sample period the TPG assumes (default 500 ns).  A frame
      with another tick is an error.
    - anode :: optional IAnodePlane type:name, gives each channel its plane
      and the set its detector ident.  Without it all channels use one
      configuration and the detector is the "detector" value (default -1).
    - tag :: process only traces with this tag (default all).
    - tpg :: the TPG configuration (keys below), and planes :: an optional
      array of per-plane objects (by WirePlaneId::index()) overriding it.
      Keys: accum_limit, pedestal_init ("zero", "mode", "fixed"), pedestal,
      filter ("none", "running_sum", "abs_running_sum"), r_value,
      scale_factor, saturate_low, filter_pedestal, threshold,
      min_samples_over_threshold, settle (samples run through the filters
      before hit finding when a channel starts, dunetrigger's AbsRunningSum
      "stabilisation_period").
    - streaming :: false (default): every frame stands alone, the state is
      reset at its start and open hits are closed at its end, as dunetrigger
      does per waveform.  true: frames are pieces of one stream; a channel's
      state carries over when its next trace continues where the last ended,
      otherwise the channel is flushed and reset.  On end of stream, open hits
      are emitted in a last, empty-coverage set.  The latency is then the
      longest TP.
 */

#ifndef WIRECELLTRIG_TPFINDER
#define WIRECELLTRIG_TPFINDER

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITriggerPrimitiveFinder.h"
#include "WireCellTrig/TPGChannel.h"

#include <map>
#include <unordered_map>

namespace WireCell::Trig {

    class TPFinder : public Aux::Logger, public ITriggerPrimitiveFinder, public IConfigurable {
       public:
        TPFinder();
        virtual ~TPFinder();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_queue& outq);

       private:
        double m_period;
        std::string m_tag{""};
        bool m_streaming{false};
        int m_detector{-1};
        size_t m_settle{0};
        IAnodePlane::pointer m_anode;
        std::vector<TPGConfig> m_plane_cfgs;  // one, or per plane index

        // streaming state
        struct Channel {
            TPGChannel tpg;
            int64_t next{0};  // absolute sample count expected next
        };
        std::unordered_map<int, Channel> m_channels;
        double m_resid{0};  // frame time - (sample count * period)
        double m_last_end{0};
        bool m_have_last{false};
        int m_count{0};

        const TPGConfig& plane_config(int channel);
        TriggerPrimitive to_tp(int channel, const ThresholdHit& hit) const;
        void flush_all(ITriggerPrimitiveSet::primitives_t& tps);
    };

}  // namespace WireCell::Trig

#endif
