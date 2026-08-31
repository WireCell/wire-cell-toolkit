#ifndef WIRECELL_CLUS_BEESINK
#define WIRECELL_CLUS_BEESINK

#include "WireCellClus/IBeeSink.h"

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/Bee.h"

#include <mutex>
#include <string>

namespace WireCell::Clus {

    /// Concrete shared Bee sink wrapping a single Bee::Sink (one open zip).
    /// Referenced by name from MultiAlgBlobClustering nodes via the "bee_sink"
    /// config.  See IBeeSink.h for the sharing/lifecycle contract.
    class BeeSink : public Aux::Logger, public IBeeSink, public IConfigurable {
       public:
        BeeSink();
        virtual ~BeeSink();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual void acquire();
        virtual void release();
        virtual size_t write(const Bee::Object& obj, size_t index, int run, int sub, int evt);

       private:
        Bee::Sink m_sink;
        int m_refs{0};
        std::mutex m_mtx;
        std::string m_outname{"mabc.zip"};
        int m_initial_index{0};
    };

}  // namespace WireCell::Clus

#endif
