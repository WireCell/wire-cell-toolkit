#ifndef WIRECELLSIO_CLUSTERFILESOURCE
#define WIRECELLSIO_CLUSTERFILESOURCE

#include "WireCellUtil/Stream.h"
#include "WireCellIface/IClusterSource.h"
#include "WireCellIface/ITerminal.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/Logger.h"


#include <boost/iostreams/filtering_stream.hpp>

#include <string>

namespace WireCell::Sio {

    /** Source cluster related information from file.
     */
    class ClusterFileSource : public Aux::Logger,
         public IClusterSource, public ITerminal,
         public IConfigurable
    {
    public:
        ClusterFileSource();
        virtual ~ClusterFileSource();

        virtual void finalize();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
        virtual bool operator()(ICluster::pointer& cluster);

    private:
        /**
           Configuration: "inname".
           
           Name the input stream.
        */
        std::string m_inname{""};

        /**
           Configuration: "prefix" (default: "cluster")

           Prefix to match against file names in the input stream.
           See ClusterFileSink for details.
         */
        std::string m_prefix{"cluster"};

        /** Configuration: "anodes".

            Provide array of IAnodePlane instances for loading JSON.
        */

      private:

        using istream_t = boost::iostreams::filtering_istream;
        istream_t m_in;

        struct header_t {
            header_t() :fname(""), fsize(0) {}
            std::string fname{""};
            size_t fsize{0};
        };
        header_t m_cur;

        ICluster::pointer load();
        ICluster::pointer load_json(int ident);
        ICluster::pointer load_numpy(int ident);
        void clear();

        std::unique_ptr<Aux::ClusterLoader> m_json_loader;

        size_t m_count{0};
        bool m_eos_sent{false};
    };

}  // namespace WireCell
#endif
