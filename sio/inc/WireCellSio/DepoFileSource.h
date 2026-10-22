/** A file stream source for depos as depo sets.
 *
 * Files should be tar streams, possibly compressed, providing a
 * stream of a pair of files: depo_data_<ident>.npy and
 * depo_info_<ident>.npy.  
 */

#ifndef WIRECELLSIO_DEPOFILESOURCE
#define WIRECELLSIO_DEPOFILESOURCE

#include "WireCellIface/IDepoSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/iostreams/filtering_stream.hpp>
#pragma GCC diagnostic pop

#include <string>
#include <vector>
#include <memory>

namespace WireCell::Sio {

    /// One archive member read ahead of its use (see
    /// DepoFileSource::next()).  Defined in the .cxx.
    struct DepoFilePending;

    class DepoFileSource : public Aux::Logger,
                           public IDepoSetSource,
                           public IConfigurable
    {
      public:
        DepoFileSource();
        virtual ~DepoFileSource();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        virtual bool operator()(IDepoSet::pointer& ds);

      private:

        /// Configuration:
        /// 
        /// Input stream name should be a file name with .tar .tar,
        /// .tar.bz2 or .tar.gz.
        ///
        /// From it, individual files named following a fixed scheme
        /// will be streamed.
        ///
        /// Depos are read from the tar stream as Numpy .npy files.
        /// Each set is the pair depo_data_<ident>.npy and
        /// depo_info_<ident>.npy, optionally followed by the sidecar
        /// depo_truth_<ident>.npy (energy, orig_id, reserved) which
        /// DepoFileSink writes when its "truth_sidecar" is true.
        /// When the sidecar is present the energy is restored,
        /// otherwise it is the SimpleDepo default (1.0).
        std::string m_inname;

        // The output stream
        boost::iostreams::filtering_istream m_in;

        // One member read ahead while looking for a sidecar; it is
        // consumed first by the next call.
        std::unique_ptr<DepoFilePending> m_pending;

        IDepoSet::pointer next();

        double m_scale{1.0};
        size_t m_count{0};
        bool m_eos_sent{false};
    };        
     
}

#endif
