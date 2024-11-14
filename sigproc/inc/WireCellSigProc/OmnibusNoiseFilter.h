/** Remove all possible noise from a microboone-like detector.
 *
 * This filter is a kitchen sink class and is a candidate for
 * factoring.
 */

#ifndef WIRECELLSIGPROC_OMNIBUSNOISEFILTER
#define WIRECELLSIGPROC_OMNIBUSNOISEFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IChannelFilter.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellAux/Logger.h"

#include <vector>
#include <map>
#include <string>

namespace WireCell {
    namespace SigProc {

        /**
           The Omnibus Noise Filter applies two series of IChannelFilter
           objects. The first series is applied on a per-channel basis and
           the second is applied on groups of channels as determined by
           its channel grouping.
        */
        class OmnibusNoiseFilter : public Aux::Logger,
                                   public WireCell::IFrameFilter,
                                   public WireCell::IConfigurable {
           public:
            typedef std::vector<std::vector<int> > grouped_channels_t;

            /// Create an OmnibusNoiseFilter.
            OmnibusNoiseFilter(std::string intag = "orig", std::string outtag = "raw");
            virtual ~OmnibusNoiseFilter();

            /// IFrameFilter interface.
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            /// IConfigurable interface.
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            // Local methods. Real code should not use these.  Expose for unit tests

            void set_channel_filters(std::vector<WireCell::IChannelFilter::pointer> filters) { m_perchan = filters; }
            void set_channel_status_filters(std::vector<WireCell::IChannelFilter::pointer> filters)
            {
                m_perchan_status = filters;
            }
            void set_grouped_filters(std::vector<WireCell::IChannelFilter::pointer> filters) { m_grouped = filters; }
            void set_channel_noisedb(WireCell::IChannelNoiseDatabase::pointer ndb) { m_noisedb = ndb; }

           private:
            // number of time ticks in the waveforms processed.  Set to 0 and first input trace sets it.
            size_t m_nticks;
            std::string m_intag, m_outtag;
            std::vector<WireCell::IChannelFilter::pointer> m_perchan, m_grouped, m_perchan_status;
            WireCell::IChannelNoiseDatabase::pointer m_noisedb;

            std::map<std::string, std::string> m_maskmap;

            size_t m_count{0};

            // This little struct, named MGCF (Multi-Group Channel Filters), is used to associate a set of
            // channels with a set of filters. See: https://github.com/WireCell/wire-cell-toolkit/issues/327
            struct MGCF {
                std::vector<IChannelNoiseDatabase::channel_group_t> channelgroups;
                std::vector<IChannelFilter::pointer> filters;
            };

            std::vector<MGCF>  m_multigroup_chanfilters;
        };

    }  // namespace SigProc
}  // namespace WireCell
#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
