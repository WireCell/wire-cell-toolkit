/** A pass-through for IFrame that saves IFrame to HDF
 */

#ifndef WIRECELLHDF5_HDF5FRAMETAP
#define WIRECELLHDF5_HDF5FRAMETAP

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

#include <hdf5.h>

namespace WireCell::Hio {

    class HDF5FrameTap : public Aux::Logger, public IFrameFilter, public IConfigurable, public ITerminal {
    public:
        HDF5FrameTap();
        virtual ~HDF5FrameTap();

        virtual bool operator()(const IFrame::pointer &inframe, IFrame::pointer &outframe);
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration &config);
        virtual void finalize();


    private:
        Configuration m_cfg;           /// copy of configuration
        IAnodePlane::pointer m_anode;  /// pointer to some APA, needed to associate chnnel ID to planes

        size_t m_calls{0};

        hid_t m_hfile;

    };

}  // namespace WireCell::Hio

#endif
