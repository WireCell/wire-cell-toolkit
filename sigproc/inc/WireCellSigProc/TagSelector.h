/**
   Make a new output frame with a set of traces selected from the
   input based on trace tags.
 */

#ifndef WIRECELLSIGPROC_TAGSELECTOR
#define WIRECELLSIGPROC_TAGSELECTOR

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <vector>

namespace WireCell {
    namespace SigProc {

        class TagSelector : public Aux::Logger,
                            public IFrameFilter,
                            public IConfigurable
        {
           public:
            TagSelector();
            virtual ~TagSelector();

            /// IFrameFilter interface.
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            /// IConfigurable interface.
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           private:
            std::vector<std::string> m_tags;
            int m_count{0};
        };
    }  // namespace SigProc
}  // namespace WireCell

#endif

