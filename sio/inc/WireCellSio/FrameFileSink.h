#ifndef WIRECELLSIO_FRAMEFILESINK
#define WIRECELLSIO_FRAMEFILESINK

#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <boost/iostreams/filtering_stream.hpp>

#include <string>
#include <vector>

namespace WireCell::Sio {

    class FrameFileSink : public Aux::Logger, public IFrameSink, public IConfigurable, public ITerminal {
    public:
        FrameFileSink();
        virtual ~FrameFileSink();

        virtual void finalize();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        virtual bool operator()(const IFrame::pointer& frame);

    private:

        /// Configuration:
        /// 
        /// Output stream name should be a file name with .tar .tar,
        /// .tar.bz2 or .tar.gz.
        ///
        /// In to it, individual files named following a fixed scheme
        /// will be streamed.
        ///
        /// Frames are written as Numpy .npy files.
        ///
        /// Future extensions to allow writing to a directory of
        /// individual files or to a zeromq stream are expected and
        /// their "outname" forms are reserved.
        std::string m_outname;

        /// The frame tag or trace tags to sink.  A tag first attempts
        /// to match a trace tag, else it attempts to match a frame
        /// tag.  If a tag is the empty string it matches untagged
        /// traces.  If no tags are given or the tag "*" is found the
        /// entire frame will be matched.  Where a trace is included
        /// in more than one match it will be saved redundantly.
        /// Note: the file retains no indication if a given tag is on
        /// a frame or trace scope.
        std::vector<std::string> m_tags;

        /// The channel mask map (CMM) can be saved to the stream as a
        /// .npy file named like "chanmask_<tag>_<ident>.npy" where
        /// <tag> is the map key, such as "bad".  Emitting the CMM may
        /// be suppressed with by setting "masks=false" in the
        /// configuration.  By default it is true.
        bool m_masks{true};


        /// The array may be transformed as:
        ///  (arr+baseline)*scale + offset
        /// This is mostly intended if digitizing/truncating.
        double m_baseline{0.0};
        double m_scale{1.0};
        double m_offset{0.0};

        /// After possible transform we may truncate from 32 bit float
        /// to 16 bit int.
        bool m_digitize{false};

        bool m_dense{false};
        int m_chbeg{0}, m_chend{0}, m_tbbeg{0}, m_tbend{0};

        // The output stream
        boost::iostreams::filtering_ostream m_out;

        void one_tag(const IFrame::pointer& frame,
                     const std::string& tag);
        void masks(const IFrame::pointer& frame);

        size_t m_count{0};
    };        
     
}

#endif
