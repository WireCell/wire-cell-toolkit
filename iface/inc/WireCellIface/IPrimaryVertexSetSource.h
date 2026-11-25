#ifndef WIRECELL_IPRIMARYVERTEXSETSOURCE
#define WIRECELL_IPRIMARYVERTEXSETSOURCE

#include "WireCellIface/ISourceNode.h"
#include "WireCellIface/IPrimaryVertexSet.h"

namespace WireCell {

    /** A source of IPrimaryVertexSet -- a producer of primary kinematics
     *  events to feed a tracking simulation (eg the edep ParticleTracking node).
     */
    class IPrimaryVertexSetSource : public ISourceNode<IPrimaryVertexSet> {
       public:
        typedef std::shared_ptr<IPrimaryVertexSetSource> pointer;

        virtual ~IPrimaryVertexSetSource();

        virtual std::string signature() { return typeid(IPrimaryVertexSetSource).name(); }

        // supply:
        // virtual bool operator()(IPrimaryVertexSet::pointer& out);
    };

}  // namespace WireCell

#endif
