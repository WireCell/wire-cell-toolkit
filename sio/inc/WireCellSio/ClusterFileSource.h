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

           Note, file name extension in the stream is used to
           determine loading mechanism.  JSON (.json) following
           cluster graph schema and Numpy (.npy) following cluster
           array schema are supported.
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
        bool load_filename();
        ICluster::pointer dispatch();
        void clear_load();

        ICluster::pointer load_json(int ident);
        ICluster::pointer load_numpy(int ident);
        void clear();

        std::unique_ptr<Aux::ClusterLoader> m_json_loader;

        size_t m_count{0};
        bool m_eos_sent{false};
        // Default true: use the imaging-time corners stored in the cluster file
        // for the dead-area bee patch instead of re-deriving them from the
        // reloaded shape (which drops RayGrid boundary layers and lets 2-view
        // dead blobs spill past the wire boundary).  Dead "corner" cloud only;
        // live/reco output is unaffected.  Set restore_corners:false to opt out.
        bool m_restore_corners{true};
        std::vector<IAnodePlane::pointer> m_anodes;
    };

}  // namespace WireCell
#endif
