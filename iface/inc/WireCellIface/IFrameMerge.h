#ifndef WIRECELL_IFRAMEMERGE
#define WIRECELL_IFRAMEMERGE

#include "WireCellIface/IHydraNode.h"
#include "WireCellIface/IFrame.h"

namespace WireCell {

    class IFrameMerge : public IHydraNodeVT<IFrame, std::tuple<IFrame>, 0> {
      public:
        virtual ~IFrameMerge();
    
        virtual std::string signature() {
            return typeid(IFrameMerge).name();
        }
        virtual std::vector<std::string> input_types() = 0;
    };

}

#endif
