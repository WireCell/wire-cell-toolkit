/** Functions and classes for ProtoDUNE VD noise filter
 *  Copied from ProtoduneHD.h and extended with ShieldCouplingSub
 */

#ifndef WIRECELLSIGPROC_PDVD
#define WIRECELLSIGPROC_PDVD

#include "WireCellSigProc/Diagnostics.h"
#include "WireCellSigProc/CoherentNoiseDump.h"

#include "WireCellIface/IChannelFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDFT.h"
#include "WireCellIface/IFilterWaveform.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Bits.h"
#include "WireCellUtil/Logging.h"

#include <set>

namespace WireCell {
    namespace SigProc {
        namespace PDVD {

            bool SignalFilter(WireCell::Waveform::realseq_t& sig);
            float CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig);
            bool RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig);

            bool RemoveFilterFlags(WireCell::Waveform::realseq_t& sig);
            bool NoisyFilterAlg(WireCell::Waveform::realseq_t& spec, float min_rms, float max_rms);
            WireCell::Waveform::realseq_t CalcMedian_shieldCoupling_u(const WireCell::IChannelFilter::channel_signals_t& chansig);
            bool Signal_mask_top_u(WireCell::Waveform::realseq_t& sig);
            std::vector<std::vector<int> > SignalProtection(WireCell::Waveform::realseq_t& sig,
                                                            const WireCell::Waveform::compseq_t& respec,
                                                            const IDFT::pointer& dft,
                                                            int res_offset,
                                                            int pad_f, int pad_b,
                                                            const WireCell::Waveform::realseq_t& time_filter_wf,
                                                            const WireCell::Waveform::realseq_t& lf_filter_wf,
                                                            float upper_decon_limit = 0.02,
                                                            float upper_adc_limit = 15,
                                                            float protection_factor = 5.0, float min_adc_limit = 50,
                                                            WireCell::SigProc::CoherentNoiseDump* dump = nullptr);
            bool Subtract_WScaling(WireCell::IChannelFilter::channel_signals_t& chansig,
                                   const WireCell::Waveform::realseq_t& medians,
                                   const WireCell::Waveform::compseq_t& respec, int res_offset,
                                   std::vector<std::vector<int> >& rois,
                                   const IDFT::pointer& dft,
                                   const WireCell::Waveform::realseq_t& time_filter_wf,
                                   const WireCell::Waveform::realseq_t& lf_filter_wf,
                                   float upper_decon_limit1 = 0.08,
                                   float roi_min_max_ratio = 0.8, float rms_threshold = 0.,
                                   WireCell::SigProc::CoherentNoiseDump* dump = nullptr);

            float get_rms_and_rois(const WireCell::Waveform::realseq_t& signal, std::vector<std::vector<int> >& rois);
            bool Is_FEMB_noise(const WireCell::IChannelFilter::channel_signals_t& chansig, int& beg, int& end, float minwidth);

            class OneChannelNoise : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                OneChannelNoise(const std::string& anode_tn = "AnodePlane",
                                const std::string& noisedb = "OmniChannelNoiseDB");
                virtual ~OneChannelNoise();

                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

                void configure(const WireCell::Configuration& config);
                WireCell::Configuration default_configuration() const;

               private:
                std::string m_anode_tn, m_noisedb_tn;
                Diagnostics::Partial m_check_partial;
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                IDFT::pointer m_dft;
                bool m_adaptive_baseline{false};
                Log::logptr_t m_log;
            };

            class CoherentNoiseSub : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                CoherentNoiseSub(const std::string& anode = "AnodePlane",
                                 const std::string& noisedb = "OmniChannelNoiseDatabase", float rms_threshold = 0.);
                virtual ~CoherentNoiseSub();

                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

               private:
                std::string m_anode_tn, m_noisedb_tn;
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                IDFT::pointer m_dft;

                float m_rms_threshold;
                std::vector<std::string> m_time_filters;
                std::string m_lf_tighter_filter;
                std::string m_lf_loose_filter;

                // Opt-in coherent-NF debug dump. When m_dump_path is empty (default)
                // the apply() hot path is bit-identical to the pre-instrumentation
                // build. When set, one .npz per group is emitted under
                // <m_dump_path>/apa<N>/<plane>_g<gid>.npz.
                std::string m_dump_path;
                std::set<int> m_dump_groups;  // empty → all groups
            };

            class FEMBNoiseSub : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                FEMBNoiseSub(const std::string& anode = "AnodePlane", float width=50.0);
                virtual ~FEMBNoiseSub();
 
                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;
 
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;
 
               private:
                std::string m_anode_tn;
                IAnodePlane::pointer m_anode;
 
                float m_width;
            };

            /**
             * Shield Coupling Removal for ProtoDUNE VD
             * Removes capacitive coupling noise from shield readout groups
             */
            class ShieldCouplingSub : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                ShieldCouplingSub(const std::string& anode = "AnodePlane",
                                 const std::string& noisedb = "OmniChannelNoiseDB",
                                 float rms_threshold = 0.);
                virtual ~ShieldCouplingSub();

                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

               private:
                std::string m_anode_tn, m_noisedb_tn;
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                std::map<int, float> m_strip_lengths;
                IDFT::pointer m_dft;

                float m_rms_threshold; 
                // bool m_capa_weight;
                // bool m_calibrated;
                // int m_group_size;
                // int m_min_channels;
            };

        }  // namespace PDVD

    }  // namespace SigProc

}  // namespace WireCell

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End: