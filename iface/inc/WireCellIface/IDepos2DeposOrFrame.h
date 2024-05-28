/** A node that takes in a queue of depo sets and outputs a queue of
 * depo sets on port 0 and an queue of frames on port 1.
 */

#ifndef WIRECELL_IDEPOS2DEPOSORFRAME
#define WIRECELL_IDEPOS2DEPOSORFRAME

#include "WireCellIface/IHydraNode.h"
#include "WireCellIface/IDepoSet.h"
#include "WireCellIface/IFrame.h"

namespace WireCell {
    
    class IDepos2DeposOrFrame : public IHydraNodeTT<std::tuple<IDepoSet>,
                                                    std::tuple<IDepoSet, IFrame>>
    {
      public:
        virtual ~IDepos2DeposOrFrame();
    
        virtual std::string signature() {
            return typeid(IDepos2DeposOrFrame).name();
        }
    };

}
#endif
