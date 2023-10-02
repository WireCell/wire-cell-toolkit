// A frame summer can "add" two frames together, possibly limiting
// them by tag.

#ifndef WIRECELL_GEN_FRAMESUMMER
#define WIRECELL_GEN_FRAMESUMMER

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFanin.h"

#include "WireCellAux/Logger.h"

#include <vector>
#include <string>

namespace WireCell {
    namespace Gen {
      class FrameSummer : Aux::Logger, public IFrameFanin, public IConfigurable {
           public:
            FrameSummer(size_t multiplicity = 2);
            virtual ~FrameSummer();

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
