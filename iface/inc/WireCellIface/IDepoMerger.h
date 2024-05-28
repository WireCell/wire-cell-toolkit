#ifndef WIRECELL_IDEPOMERGER
#define WIRECELL_IDEPOMERGER

#include "WireCellIface/IHydraNode.h"
#include "WireCellIface/IDepo.h"

namespace WireCell {

    class IDepoMerger : public IHydraNodeTT<std::tuple<IDepo, IDepo>, std::tuple<IDepo> > {
       public:

        using IHydraNodeTT<std::tuple<IDepo, IDepo>, std::tuple<IDepo> >::input_queues;
        using IHydraNodeTT<std::tuple<IDepo, IDepo>, std::tuple<IDepo> >::output_queues;

        typedef std::shared_ptr<IDepoMerger> pointer;
        virtual ~IDepoMerger();

        virtual std::string signature() { return typeid(IDepoMerger).name(); }

        // subclass supply:
        // virtual bool operator()(input_queues& inqs,
        //                         output_queues& outqs);
    };

}  // namespace WireCell
#endif
