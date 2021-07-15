/** Some frames to a Numpy file */

#ifndef WIRECELLSIO_NUMPYFRAMESAVER
#define WIRECELLSIO_NUMPYFRAMESAVER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"

#include "WireCellSio/Cfg/NumpyFrameSaver/Structs.hpp"

namespace WireCell {
    namespace Sio {

        // This saver immediately saves each frame.
        class NumpyFrameSaver : public virtual WireCell::IFrameFilter, public WireCell::IConfigurable {
           public:
            NumpyFrameSaver();
            virtual ~NumpyFrameSaver();

            /// IFrameFilter
            virtual bool operator()(const WireCell::IFrame::pointer& inframe, WireCell::IFrame::pointer& outframe);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

           private:
            using config_t = WireCellSio::Cfg::NumpyFrameSaver::Config;
            config_t m_cfg;
            
            int m_save_count;  // count frames saved
            Log::logptr_t l;
        };
    }  // namespace Sio
}  // namespace WireCell
#endif
