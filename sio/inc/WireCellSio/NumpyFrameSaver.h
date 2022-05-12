/** Save frames to a Numpy file */

#ifndef WIRECELLSIO_NUMPYFRAMESAVER
#define WIRECELLSIO_NUMPYFRAMESAVER

#include "WireCellAux/Configurable.h"
#include "WireCellAux/Logger.h"

#include "WireCellIface/IFrameFilter.h"

#include "WireCellSio/Cfg/NumpyFrameSaver.hpp"

namespace WireCell {
    namespace Sio {

        using WireCell::Sio::Cfg::NumpyFrameSaver::Config;

        // This saver immediately saves each frame.
        class NumpyFrameSaver : public Aux::Logger,
                                public Aux::Configurable<Config>,
                                public virtual WireCell::IFrameFilter
        {
           public:
            NumpyFrameSaver();
            virtual ~NumpyFrameSaver();

            /// IFrameFilter
            virtual bool operator()(const WireCell::IFrame::pointer& inframe, WireCell::IFrame::pointer& outframe);

            /// IConfigurable
            // virtual WireCell::Configuration default_configuration() const;
            // virtual void configure(const WireCell::Configuration& config);
            virtual void configured();

           private:
            // using config_t = WireCellSio::Cfg::NumpyFrameSaver::Config;
            // config_t m_cfg;
            
            int m_save_count{0};  // count frames saved
            // Log::logptr_t l;
        };
    }  // namespace Sio
}  // namespace WireCell
#endif
