/** Shadow Ghosting

    This component is a cluster filter performing "deghosting".

    It will remove blobs which are likely to not contain true charge.

    It starts by forming geometric clusters and finding their mutual
    shadows in each view.  See aux/docs/cluster-shadow.org for entry
    point to documentation on this operation.

    It then classifies entire clusters as ghosts (or not) based on
    their shadowing.

 */

#ifndef WIRECELL_IMG_SHADOWGHOSTING
#define WIRECELL_IMG_SHADOWGHOSTING

#include "WireCellIface/IClusterFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Img {

    class ShadowGhosting : public Aux::Logger, public IClusterFilter, public IConfigurable {
      public:
        ShadowGhosting();
        virtual ~ShadowGhosting();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

      private:

        // Config: "shadow_type" is "c"/"channel" or "w"/"wire" and
        // determines the blob shadow type.
        std::string m_shadow_type{"channel"};

        int m_count{0};

    };

}

#endif
