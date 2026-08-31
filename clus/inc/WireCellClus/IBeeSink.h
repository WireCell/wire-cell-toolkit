/** A shared Bee sink: a single open Bee zip that multiple
    MultiAlgBlobClustering nodes write into so the whole chain produces one
    .zip instead of one per node.

    Writes are placed at an explicit event index (NOT Bee::Sink's
    name-collision auto-increment, which desyncs across multiple writers).  The
    zip is closed via reference counting: each client acquire()s in configure
    and release()s in finalize; the underlying store is closed on the last
    release(), independent of node finalize ordering.
 */
#ifndef WIRECELL_CLUS_IBEESINK
#define WIRECELL_CLUS_IBEESINK

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Bee.h"

namespace WireCell::Clus {

    class IBeeSink : public IComponent<IBeeSink> {
       public:
        virtual ~IBeeSink() {}

        /// Register/unregister a client.  The underlying zip is closed on the
        /// last release().
        virtual void acquire() = 0;
        virtual void release() = 0;

        /// Write one Bee object at an explicit event index with the given RSE.
        virtual size_t write(const Bee::Object& obj, size_t index, int run, int sub, int evt) = 0;
    };

}  // namespace WireCell::Clus

#endif
