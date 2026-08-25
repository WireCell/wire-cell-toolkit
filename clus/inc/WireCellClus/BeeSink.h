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
        /// Zip name.  A printf conversion in it (eg "mabc-pr-evt%d.zip") means
        /// ONE ZIP PER EVENT: the open zip is closed and a new one opened
        /// whenever write() is handed a different event number.  This is what
        /// lets a process that streams many events write the same per-event
        /// files a one-event-per-process job writes.  No '%' => one zip for the
        /// whole process, opened in configure(), exactly as before.
        std::string m_outname{"mabc.zip"};
        int m_initial_index{0};
        bool m_templated{false};
        /// Event number of the currently open zip; -1 when none is open.
        /// Only meaningful when m_templated.
        int m_open_evt{-1};
    };

}  // namespace WireCell::Clus

#endif
