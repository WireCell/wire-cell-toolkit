#ifndef WIRECELL_SIO_TENSORFILESOURCE
#define WIRECELL_SIO_TENSORFILESOURCE

#include "WireCellIface/ITensorSetSource.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"


#include <boost/iostreams/filtering_stream.hpp>

namespace WireCell::Sio {

    class TensorFileSource : public Aux::Logger, public ITensorSetSource,
                             public WireCell::IConfigurable, public WireCell::ITerminal

    {
      public:
        TensorFileSource();
        virtual ~TensorFileSource();

        // IConfigurable
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& config);

        // ITerminal
        virtual void finalize();

        // ITensorSetSource
        virtual bool operator()(ITensorSet::pointer &out);

        /** Config: "inname"

            Name the input stream container.

            The stream container format is determined by the name
            suffix.  The usual container formats are accepted
            including tar with .tar suffix and with optional .gz/.bz2
            compression or zip with the .zip or the .npz suffix.

            The container stream will should contain a mix of files in
            .npy format holding ITensor arrays and .json format
            holding ITensorSet and ITensor metadata objects.

            Array .npy and metadata .json files may arrive in any
            order but all that correspond to a single ITensorSet must
            be contiguous in the stream.

        */
        std::string m_inname{""};

        /** Config: "prefix"

            Name a prefix to assume when matching names of files in
            the stream.  Any not matched are ignored.  

            Three patterns of file names are matched:

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

      private:
        
        using istream_t = boost::iostreams::filtering_istream;
        istream_t m_in;
        size_t m_count{0};
        bool m_eos_sent{false};

        ITensorSet::pointer load();
        bool read_head();
        void clear();

        // Must read the custard head to determine if we the stream
        // has gone past the last <ident> and thus cache it for next
        // time around.
        struct header_t {
            header_t() :fname(""), fsize(0) {}
            std::string fname{""};
            size_t fsize{0};
        };
        header_t m_cur;

    };

}

#endif
