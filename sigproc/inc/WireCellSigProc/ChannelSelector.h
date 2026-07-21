/**
   Make a new output frame with a set of traces selected from the
   input based on being in a set of channels and possibly tags.
 */

#ifndef WIRECELLSIGPROC_CHANNELSELECTOR
#define WIRECELLSIGPROC_CHANNELSELECTOR

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/TagRules.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <vector>
#include <unordered_set>

namespace WireCell {
    namespace SigProc {

        class ChannelSelector : public Aux::Logger,
                                public IFrameFilter,
                                public IConfigurable
        {
           public:
            ChannelSelector();
            virtual ~ChannelSelector();

            /// IFrameFilter interface.
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            /// IConfigurable interface.
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           protected:
            virtual void set_channels(const std::vector<int>& channels);
            void parse_channels(const Json::Value& jchannels);

           private:
            std::vector<std::string> m_tags;
            std::unordered_set<int> m_channels;
            int m_count{0};
            tagrules::Context m_ft;
            bool m_use_rules{false};
            bool m_select_all_channels{false};
        };
    }  // namespace SigProc
}  // namespace WireCell

#endif
