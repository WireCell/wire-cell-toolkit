/**
   Make a new output frame selecting just traces from the configured
   wire plane.
 */

#ifndef WIRECELLAUX_PLANESELECTOR
#define WIRECELLAUX_PLANESELECTOR

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/TagRules.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <vector>
#include <unordered_set>

namespace WireCell::Aux {

    class PlaneSelector : public Aux::Logger,
                          public IFrameFilter,
                          public IConfigurable
    {
      public:
        PlaneSelector();
        virtual ~PlaneSelector();

        /// IFrameFilter interface.
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        /// IConfigurable interface.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

      private:
        std::vector<std::string> m_tags;
        std::vector<std::string> m_summary_tags;
        tagrules::Context m_ft;
        int m_count{0};

        // The channel ident numbers we will select
        std::unordered_set<int> m_chids;
    };

}  // namespace WireCell::Aux

#endif
