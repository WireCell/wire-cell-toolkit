/** Functions and classes for ProtoDUNE HD noise filter
 *  Modified from Protodune.h
 */

#ifndef WIRECELLSIGPROC_PDHD
#define WIRECELLSIGPROC_PDHD

#include "WireCellSigProc/Diagnostics.h"

#include "WireCellIface/IChannelFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDFT.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Bits.h"

namespace WireCell {
    namespace SigProc {
        namespace PDHD {

            bool SignalFilter(WireCell::Waveform::realseq_t& sig);
            float CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig);
            bool RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig);

            bool RemoveFilterFlags(WireCell::Waveform::realseq_t& sig);
            bool NoisyFilterAlg(WireCell::Waveform::realseq_t& spec, float min_rms, float max_rms);

            double filter_low(double freq, double cut_off = 0.08);
            double filter_time(double freq);
            double filter_low_loose(double freq);
            std::vector<std::vector<int> > SignalProtection(WireCell::Waveform::realseq_t& sig,
                                                            const WireCell::Waveform::compseq_t& respec,
                                                            const IDFT::pointer& dft,
                                                            int res_offset,
                                                            int pad_f, int pad_b, float upper_decon_limit = 0.02,
                                                            float decon_lf_cutoff = 0.08, float upper_adc_limit = 15,
                                                            float protection_factor = 5.0, float min_adc_limit = 50);
            bool Subtract_WScaling(WireCell::IChannelFilter::channel_signals_t& chansig,
                                   const WireCell::Waveform::realseq_t& medians,
                                   const WireCell::Waveform::compseq_t& respec, int res_offset,
                                   std::vector<std::vector<int> >& rois,
                                   const IDFT::pointer& dft,
                                   float upper_decon_limit1 = 0.08,
                                   float roi_min_max_ratio = 0.8, float rms_threshold = 0.);

            class OneChannelNoise : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                OneChannelNoise(const std::string& anode_tn = "AnodePlane",
                                const std::string& noisedb = "OmniChannelNoiseDB");
                virtual ~OneChannelNoise();

                // IChannelFilter interface

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

                void configure(const WireCell::Configuration& config);
                WireCell::Configuration default_configuration() const;

               private:
                std::string m_anode_tn, m_noisedb_tn;
                Diagnostics::Partial m_check_partial;  // at least need to expose them to configuration
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                IDFT::pointer m_dft;
            };

            /**
             * Coherent Noise Removal (CNR) for ProtoDUNE HD
             * Modified from Microboone
             */
            class CoherentNoiseSub : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                CoherentNoiseSub(const std::string& anode = "AnodePlane",
                                 const std::string& noisedb = "OmniChannelNoiseDB", float rms_threshold = 0.);
                virtual ~CoherentNoiseSub();

                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

                //// IChannelFilter interface

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

               private:
                std::string m_anode_tn, m_noisedb_tn;
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                IDFT::pointer m_dft;

                float m_rms_threshold;
            };

        }  // namespace PDHD

    }  // namespace SigProc

}  // namespace WireCell

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
