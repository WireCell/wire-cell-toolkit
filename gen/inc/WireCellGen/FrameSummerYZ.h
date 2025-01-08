// A frame summer can "add" two frames together, possibly limiting
// them by tag.

#ifndef WIRECELL_GEN_FRAMESUMMERYZ
#define WIRECELL_GEN_FRAMESUMMERYZ

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFanin.h"

#include "WireCellAux/Logger.h"

#include <vector>
#include <string>

namespace WireCell {
    namespace Gen {
      class FrameSummerYZ : Aux::Logger, public IFrameFanin, public IConfigurable {
           public:
            FrameSummerYZ(size_t multiplicity = 2);
            virtual ~FrameSummerYZ();

	    virtual std::vector<std::string> input_types();

            // IJoinNode
	    virtual bool operator()(const input_vector& intup, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           private:
	    size_t m_multiplicity; 
	    int m_count{0};

        };
    }  // namespace Gen
}  // namespace WireCell

#endif
