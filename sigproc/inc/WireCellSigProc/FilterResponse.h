/** This component provides per-wire filtering of
 * field responses based on a configuration data file. */

#ifndef WIRECELLSIGPROC_FILTERRESPONSE
#define WIRECELLSIGPROC_FILTERRESPONSE

#include "WireCellIface/IChannelResponse.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Units.h"

#include <string>
#include <unordered_map>

namespace WireCell {
    namespace SigProc {
        class FilterResponse : public IChannelResponse, public IConfigurable {
           public:
            FilterResponse(const char* filename = "", const int planeid = 0);

            virtual ~FilterResponse();

            // IChannelResponse
            virtual const Waveform::realseq_t& channel_response(int channel_ident) const;
            virtual Binning channel_response_binning() const;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           private:
            std::string m_filename;
            int m_planeid;
            std::unordered_map<int, Waveform::realseq_t> m_cr;
            Binning m_bins;
        };

    }  // namespace SigProc

}  // namespace WireCell
#endif
