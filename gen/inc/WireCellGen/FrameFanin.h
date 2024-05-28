/** Fanin N frames and output one frame.

    The FrameFanin has two operational modes.

    Fixed multiplicity mode is set by a non-zero "multiplicity" value.  This
    mode is appropriate for FrameFanin to operate in a DFP graph.

    Dynamic multiplicity mode is set by a zero "multiplicity" value.  This mode
    is relevant for using FrameFanin as a tool inside another DFP node which may
    generate a variable set of frames that need to be merged into one.  In this
    mode it will accept any size input vector of frames.  Arrays of trace tags
    and tag rules may still be supplied.  The arrays are mapped to the vector of
    input frames similarly to the case of fixed multiplicity.  However, any
    frames at indices in the vector equal to or larger than the size of an array
    will use the last element in that array.

    See also comments around the configuration parameters.

 */

#ifndef WIRECELL_GEN_FRAMEFANIN
#define WIRECELL_GEN_FRAMEFANIN

#include "WireCellIface/IFrameFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/TagRules.h"
#include "WireCellAux/Logger.h"

#include <vector>
#include <string>

namespace WireCell {
    namespace Gen {

        // Fan in N frames to one.
        class FrameFanin : Aux::Logger,
                           public IFrameFanin,
                           public IConfigurable
        {
          public:
            FrameFanin(size_t multiplicity = 2);
            virtual ~FrameFanin();

            virtual std::vector<std::string> input_types();

            virtual bool operator()(const input_vector& invec, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

          private:

            std::string tag(size_t port) const;

            size_t m_multiplicity{0};
            std::vector<std::string> m_tags;
            tagrules::Context m_ft;
            int m_count{0};
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
