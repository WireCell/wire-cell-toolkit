/** A pass-through for IFrame that saves IFrame to HDF
 */

#ifndef WIRECELLHDF5_HDF5FRAMETAP
#define WIRECELLHDF5_HDF5FRAMETAP

//#include "WireCellIface/IAnodePlane.h"
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
        /** Configuration parameters.  All are optional.

            filename (string) - Name the HDF5 file to produce.  Default: "wct-frame.hdf5"

            trace_tags (list of string) - The tagged traces to save.  Default: whole frame

            digitize (bool) - If true, truncate to int16.

            baseline (float) - default 0, see below
            scale (float) - default 1, see below
            offset (float) - default 0, see below

            These three apply linear transform to trace samples prior to digitize:

                (sample + baseline)*scale + offset

            Beware that baseline and offset are both additive.

            tick0 (int) - force the output array to begin at this tick.  Default: not set.
            nticks (int) - force the otput array to have this number of columns. Default: not set.
            
            gzip (int, [0,...,9]) - The gzip compression level.  Zero is the
            default and compression is applied.

            chunk (int or 2-int array) - The chunk size along each dimension used for gzip.
         */
        float m_baseline{0.0}, m_scale{1.0}, m_offset{0.0};
        bool m_digitize{false};
        int m_gzip{0};
        std::vector<hsize_t> m_chunk = {256,256};

        Configuration m_cfg;           /// copy of configuration
        //IAnodePlane::pointer m_anode;  /// pointer to some APA, needed to associate chnnel ID to planes

        size_t m_calls{0};

        hid_t m_hfile;

    };

}  // namespace WireCell::Hio

#endif
