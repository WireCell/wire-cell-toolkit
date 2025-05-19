/** Inherently MicroBooNE-specific functions and classes
 */

#ifndef WIRECELLSIGPROC_MICROBOONE
#define WIRECELLSIGPROC_MICROBOONE

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
        namespace Microboone {

            bool Chirp_raise_baseline(WireCell::Waveform::realseq_t& sig, int bin1, int bin2);
            bool SignalFilter(WireCell::Waveform::realseq_t& sig);
            float CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig);
            bool RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig);

            bool RemoveFilterFlags(WireCell::Waveform::realseq_t& sig);
            bool NoisyFilterAlg(WireCell::Waveform::realseq_t& spec, float min_rms, float max_rms);

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
                                   float roi_min_max_ratio = 0.8, float rms_threshold = 0.,
                                   float correlation_threshold = 4.0, float default_scaling = 0.);

            // hold common config stuff
            class ConfigFilterBase : public WireCell::IConfigurable {
               public:
                ConfigFilterBase(const std::string& anode = "AnodePlane",
                                 const std::string& noisedb = "OmniChannelNoiseDB");
                virtual ~ConfigFilterBase();

                /// IConfigurable configuration interface
                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

                // FIXME: this method needs to die.
                void set_channel_noisedb(WireCell::IChannelNoiseDatabase::pointer ndb) { m_noisedb = ndb; }

               protected:
                std::string m_anode_tn, m_noisedb_tn;
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                IDFT::pointer m_dft;
            };

            /** Microboone style coherent noise subtraction.
             *
             * Fixme: in principle, this class could be general purpose
             * for other detectors.  However, it uses the functions above
             * which hard code microboone-isms.  If those
             * microboone-specific parameters can be pulled out to a
             * higher layer then this class can become generic and move
             * outside of this file.
             */
            class CoherentNoiseSub : public WireCell::IChannelFilter, public ConfigFilterBase {
               public:
                CoherentNoiseSub(const std::string& anode = "AnodePlane",
                                 const std::string& noisedb = "OmniChannelNoiseDB", 
                                 float rms_threshold = 0., 
                                 float correlation_threshold = 4.,
                                 float default_scaling = 0.);
                virtual ~CoherentNoiseSub();

                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

                //// IChannelFilter interface

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

               private:
                float m_rms_threshold;
                float m_correlation_threshold;
                float m_default_scaling;
            };

            /** Microboone style single channel noise subtraction.
             *
             * Fixme: in principle, this class could be general purpose
             * for other detectors.  However, it uses the functions above
             * which hard code microboone-isms.  If those
             * microboone-specific parameters can be pulled out to a
             * higher layer then this class can become generic and move
             * outside of this file.
             */

            class OneChannelNoise : public WireCell::IChannelFilter, public ConfigFilterBase {
               public:
                OneChannelNoise(const std::string& anode_tn = "AnodePlane",
                                const std::string& noisedb = "OmniChannelNoiseDB");
                virtual ~OneChannelNoise();

                //// IChannelFilter interface

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

               private:
                Diagnostics::Chirp m_check_chirp;      // fixme, these should be done via service interfaces
                Diagnostics::Partial m_check_partial;  // at least need to expose them to configuration
            };

            class OneChannelStatus : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                OneChannelStatus(const std::string anode_tn = "AnodePlane", double threshold = 3.5, int window = 5,
                                 int nbins = 250, double cut = 14);
                virtual ~OneChannelStatus();

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

                bool ID_lf_noisy(signal_t& sig) const;

               private:
                std::string m_anode_tn;
                IAnodePlane::pointer m_anode;
                IDFT::pointer m_dft;
                double m_threshold;
                int m_window;
                int m_nbins;
                double m_cut;
            };

            class ADCBitShift : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                ADCBitShift(int nbits = 12, int exam_nticks = 500, double threshold_sigma = 7.5,
                            double threshold_fix = 0.8);
                virtual ~ADCBitShift();

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

                virtual void configure(const WireCell::Configuration& config);
                virtual WireCell::Configuration default_configuration() const;

               private:
                // Number of bits 12
                int m_nbits;
                // How many ADC: 500
                int m_exam_nticks;
                // Threshold for ADC Bit shift: 7.5 sigma
                double m_threshold_sigma;
                // Threshold for correction: 80%
                double m_threshold_fix;
            };

        }  // namespace Microboone

    }  // namespace SigProc

}  // namespace WireCell

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
