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

namespace WireCell::Sio {

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
        std::string m_inname;

        // The output stream
        boost::iostreams::filtering_istream m_in;

        IDepoSet::pointer next();

        double m_scale{1.0};
        size_t m_count{0};
        bool m_eos_sent{false};
    };        
     
}

#endif
