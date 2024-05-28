/**
   DeposOrBust - a pass through non-empty depo sets on port 0 else an
   empty frame on port 1.  Frame has same ident as the empty deposet
   but is otherwise completely devoid of any info.

   This provides a "branch" pattern which may be used to circumvent a
   pipeline that consumes depos and produces frames with a second path
   that immediately produces empty frames given empty depos.  See
   FrameSync for a "merge" pattern to match.
*/

#ifndef WIRECELLGEN_DEPOSORBUST
#define WIRECELLGEN_DEPOSORBUST 

#include "WireCellIface/IDepos2DeposOrFrame.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Gen {

    class DeposOrBust : Aux::Logger, public WireCell::IDepos2DeposOrFrame
    {
      public:

        DeposOrBust();
        virtual ~DeposOrBust();

        // hydra
        virtual bool operator()(input_queues& iqs, output_queues& oqs);

    };
    
}

#endif
