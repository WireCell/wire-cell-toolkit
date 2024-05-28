/**
 * Functions and classes for correcting response functions in-situ
 */

#ifndef WIRECELLSIGPROC_RESPONSESHAPER
#define WIRECELLSIGPROC_RESPONSESHAPER

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
        namespace ResponseShaper {

            /**
             * Correct variations in per-channel response
             */
            class OneChannelResponse : public WireCell::IChannelFilter, public WireCell::IConfigurable {
               public:
                OneChannelResponse(const std::string& anode_tn = "AnodePlane",
                                const std::string& noisedb = "OmniChannelNoiseDB");
                virtual ~OneChannelResponse();

                // IChannelFilter interface

                /** Filter in place the signal `sig` from given `channel`. */
                virtual WireCell::Waveform::ChannelMaskMap apply(int channel, signal_t& sig) const;

                /** Filter in place a group of signals together. */
                virtual WireCell::Waveform::ChannelMaskMap apply(channel_signals_t& chansig) const;

                void configure(const WireCell::Configuration& config);
                WireCell::Configuration default_configuration() const;

               private:
                std::string m_elecresponse_tn{"ColdElec"};
                std::string m_anode_tn, m_noisedb_tn;
                IAnodePlane::pointer m_anode;
                IChannelNoiseDatabase::pointer m_noisedb;
                IDFT::pointer m_dft;
            };

        }  // namespace ResponseShaper

    }  // namespace SigProc

}  // namespace WireCell

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
