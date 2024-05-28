#ifndef WIRECELLSIO_FRAMEFILESOURCE
#define WIRECELLSIO_FRAMEFILESOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/custard/pigenc.hpp"

#include <boost/iostreams/filtering_stream.hpp>

#include <string>
#include <vector>

namespace WireCell::Sio {

    class FrameFileSource : public Aux::Logger, public IFrameSource, public IConfigurable {
    public:
        FrameFileSource();
        virtual ~FrameFileSource();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        virtual bool operator()(IFrame::pointer& frame);
        
      private:
        
        /// Configuration:

        /** Config: "inname"
      
            Input stream name should be a file name with .tar .tar,
            .tar.bz2 or .tar.gz.
        
            From it, individual files named following a fixed scheme
            will be streamed.
        
            Frames are read from the tar stream as Numpy .npy files.
        */
        std::string m_inname{""};

        /** Config: "tags".

            A set of tags selecting "framelet" array to include as tagged
            traces.

            If the <tag> portion of the framelet array file name, eg

                frame_<tag>_<ident>

            is found in this set the framelet will provide tagged traces with
            the tag <tag>.

            Note, the use of a <tag> that is empty "" or the special "*" is
            intended for use only for the context of a frame file to indicate
            untagged traces.  In the context of an IFrame, it is against
            convention to explicitly tag traces with either label.

        */
        std::vector<std::string> m_tags;

        /** Config: "frame_tags".

            Apply these tags to the produced frame.
        */
        std::vector<std::string> m_frame_tags;        

        // The output stream
        boost::iostreams::filtering_istream m_in;

        IFrame::pointer load();

        // Classify a <tag> label from a framelet array name. 
        bool is_tagged(const std::string& tag);
        bool is_untagged(const std::string& tag);
        bool is_excluded(const std::string& tag);

        size_t m_count{0};
        bool m_eos_sent{false};

        // We must read-ahead one array to know that the sequence of
        // arrays for a given frame ident is complete.  This stashes
        // that extra read for use next time we are called.
        struct entry_t {
            entry_t() :fname(""), fsize(0), okay(false), ident(-1), type(""), tag(""), ext(""), pig() {}
            std::string fname{""};
            size_t fsize{0};
            bool okay{false};
            int ident{-1};
            std::string type{""}, tag{""}, ext{""};
            pigenc::File pig;
        };
        entry_t m_cur;

        void clear();
        bool read();

    };        
     
}

#endif
