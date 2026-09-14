#ifndef WIRECELLEDEP_IONIZATION
#define WIRECELLEDEP_IONIZATION

#include "WireCellIface/ISimTruthSegments.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Edep {

    /** Ionization: ISimTruth -> ITrackSegmentSet.
     *
     *  A "stripper" that pulls the ionization track segments out of an
     *  ISimTruth (as produced by ParticleTracking) and re-emits them for the
     *  ionization chain (eg Gen::TrackSegmentSampler -> depos).  Because
     *  ISimTruth returns its segments as a shared pointer to the same object a
     *  standalone producer would emit, this is a pass-through with no copy.
     */
    class Ionization : public Aux::Logger, public ISimTruthSegments {
       public:
        Ionization();
        virtual ~Ionization();

        // ISimTruthSegments (IFunctionNode<ISimTruth, ITrackSegmentSet>)
        virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell::Edep

#endif
