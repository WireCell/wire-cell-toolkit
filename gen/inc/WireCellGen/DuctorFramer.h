/**
   Apply a ductor across a depo set to produce a frame.

   This delegates to two sub DFP nodes:

   - an IDuctor is fed a stream of IDepos from the input IDepoSet and produces
     zero or more IFrames.

   - if the set of IFrame produced over the IDepoSet has more than one entry,
     the set is fed an IFrameFanin in order to merge them into a single IFrame.
     Note the IFrameFanin must be able to handle DYNAMIC MULTIPLICITY.  In the
     most likely case of using the concrete FrameFanin class for this, the
     multiplicity should be set to 0.
 */

#ifndef WIRECELLGEN_DUCTORFRAMER
#define WIRECELLGEN_DUCTORFRAMER

#include "WireCellIface/IDepoFramer.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFanin.h"
#include "WireCellIface/IDuctor.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Gen {

    class DuctorFramer : public Aux::Logger,
                         public IDepoFramer,
                         public IConfigurable
    {
      public:

        DuctorFramer();
        virtual ~DuctorFramer();

        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

      private:

        IDuctor::pointer m_ductor{nullptr};
        IFrameFanin::pointer m_fanin{nullptr};
        size_t m_count{0};
    };
}

#endif


