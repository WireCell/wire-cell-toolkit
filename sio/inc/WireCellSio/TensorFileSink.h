#ifndef WIRECELL_SIO_TENSORFILESINK
#define WIRECELL_SIO_TENSORFILESINK

#include "WireCellIface/ITensorSetSink.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/iostreams/filtering_stream.hpp>
#pragma GCC diagnostic pop

namespace WireCell::Sio {

    class TensorFileSink : public Aux::Logger, public ITensorSetSink,
                           public WireCell::IConfigurable, public WireCell::ITerminal

    {
      public:
        TensorFileSink();
        virtual ~TensorFileSink();

        // IConfigurable
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& config);

        // ITerminal
        virtual void finalize();

        // ITensorSetSink
        virtual bool operator()(const ITensorSet::pointer &in);        

        /** Config: "outname"

            Name the output stream container.

            The stream container format is determined by the name
            suffix.  The usual container formats are accepted
            including tar with .tar suffix and with optional .gz/.bz2
            compression or zip with the .zip or the .npz suffix.

            The container stream will receive a mix of files in .npy
            format holding ITensor arrays and .json format holding
            ITensorSet and ITensor metadata objects.

            Note, when producing a .npz file and using numpy.load(),
            the content of the contained JSON files will be returned
            as raw bytes.  To get Python objects for these contained
            files use code like:

            >>> dat = json.loads(numpy.load("file.npz")["name.json"])

            Or:

            >>> dat = wirecell.util.ario.load("file.npz")["name"]

            If the name contains a printf conversion (a '%', eg
            "pctree-pr-evt%d.tar.gz") the sink writes ONE CONTAINER PER
            ITensorSet IDENT instead of one container for the whole stream: the
            current container is closed and a new one opened, named
            String::format(outname, ident), whenever the ident changes.  This
            lets a single process that streams many events write the same
            per-event files a one-event-per-process job writes.  A name with no
            '%' behaves exactly as before -- one container, opened in
            configure() -- so every existing job is unchanged.
        */
        std::string m_outname{"tensors.npz"};


        /** Config: "prefix"

            Name a prefix to be prepended to each contained file name
            sent out to the container stream.

            Three types of files are produced:

            <prefix>tensorset_<ident>_metadata.json 
            <prefix>tensor_<ident>_<index>_metadata.npy
            <prefix>tensor_<ident>_<index>_array.json

            Where <ident> is the value from ITensorSet::ident() and
            <index> is the index at which the ITensor is found in the
            ITensorSet.  Characters not surrounded by <>'s are literal.

            Note, no "_" is implicitly added between prefix string and
            the remainder.  Include it in the prefix if it is wanted.
        */
        std::string m_prefix{""};

        // if true, disable writing to file which makes it faster to debug
        bool m_dump_mode{false};

      private:
        
        using ostream_t = boost::iostreams::filtering_ostream;
        ostream_t m_out;
        size_t m_count{0};

        // Set when m_outname carries a printf conversion: one container per
        // ident (see the "outname" doc above).
        bool m_templated{false};
        // Whether m_out currently holds an open container, and for which
        // ident.  Only meaningful when m_templated.
        bool m_open{false};
        int m_open_ident{-1};

        /// Close any open container and open the one for this ident.
        void reopen(int ident);
        /// Close the open container, if any.
        void closeout();

        void numpyify(ITensor::pointer ten, const std::string& fname);
        void jsonify(const Configuration& cfg, const std::string& fname);


    };

}

#endif
