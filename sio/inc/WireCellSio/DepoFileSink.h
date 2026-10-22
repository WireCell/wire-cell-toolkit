/**
   Sink stream depo sets to file.
 */

#ifndef WIRECELLSIO_DEPOFILESINK
#define WIRECELLSIO_DEPOFILESINK

#include "WireCellIface/IDepoSetSink.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/iostreams/filtering_stream.hpp>
#pragma GCC diagnostic pop

#include <string>
#include <vector>

namespace WireCell::Sio {

    class DepoFileSink : public Aux::Logger, public IDepoSetSink,
         public IConfigurable, public ITerminal {
    public:
        DepoFileSink();
        virtual ~DepoFileSink();

        // to close the stream
        virtual void finalize();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        // sink
        virtual bool operator()(const IDepoSet::pointer& deposet);

    private:

        /// Configuration:
        /// 
        /// Output stream name should be a file name with .tar .tar,
        /// .tar.bz2 or .tar.gz.
        ///
        /// In to it, individual files named following a fixed scheme
        /// will be streamed.
        ///
        /// Each depo set is written as two Numpy arrays each in their
        /// own .npy file in the output file stream.
        ///
        /// The arrays/files are name: depo_data_<ident>.npy and
        /// depo_info_<ident>.npy with <ident> coming from
        /// IDepoSet::ident().
        ///
        /// See @ref:WireCell::Aux::fill from DepoTools for array
        /// details.
        std::string m_outname{"wct-depos.tar.bz2"};

        /// If true, a third array depo_truth_<ident>.npy, double
        /// shaped (ndepos, 3) = (energy, orig_id, reserved), is
        /// written after depo_info (see Aux::fill_truth).  Default
        /// false: the output stream is byte-identical to before the
        /// knob existed.
        bool m_truth_sidecar{false};

        // The output stream
        boost::iostreams::filtering_ostream m_out;

        size_t m_count{0};
    };        
     
}

#endif
