/** Base API for interfacing with Bee.

    See WireCell::Aux::Bee namespace for higher level data type support.

    References:

    https://www.phy.bnl.gov/twister/bee/

    https://bnlif.github.io/wire-cell-docs/viz/uploads/

    For list of canonical Bee "detector" names see end of:

    https://github.com/WireCell/wire-cell-bee3/blob/main/events/static/js/bee/physics/experiment.js

    These are conventionally the same as how LArSoft defines its detector name.

    Not in docs:

    - "nq" key.  cluster number?
 
*/

#ifndef WIRECELLUTIL_BEE
#define WIRECELLUTIL_BEE

#include "WireCellUtil/Stream.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Point.h"

#include <boost/filesystem.hpp>

#include <unordered_set>

namespace WireCell::Bee {

    /// A Bee object maps to one Bee JSON file.
    class Object {

        // Note: exceptionally for Wire-Cell toolkit, any distance values
        // (x,y,z) held in m_data here are in Bee units (cm), not WCT units.
        Configuration m_data;

    public:

        /// Construct with metadata
        Object(const std::string& geom="uboone", // canonical detector geometry name
               const std::string& type="anonymous", // name of the "algorithm"
               int run=0, int sub=0, int evt=0);


        /// Drop arrays and reset the event number
        void reset(int evt);

        void detector(const std::string& geom);
        void algorithm(const std::string& type);
        std::string algorithm() const;

        /// Set the run/subrun/event numbers
        void rse(int run, int sub, int evt);

        std::vector<int> rse() const;

        // Add a WCT point.  Note, p is expected to be in usual WCT system-of-units.
        void append(const Point& p, double q=0);

        // Return self as a JSON string
        std::string json() const;

        void append(const Object& obj);

        size_t size() const;
        bool empty() const;

    public:

    };

    // todo: source

    /// A Bee Sink persists Objects to some store.
    class Sink {

        boost::iostreams::filtering_ostream m_out;
        size_t m_index{0};
        size_t m_rse{0};        // a hash
        std::unordered_set<std::string> m_seen;

    public:

        /// Construct a sink of Bee objects.  The "store" indicates some
        /// persistent resource to receive the objects.  This is typically a
        /// file name ending in ".zip".
        Sink();
        Sink(const std::string& store);

        /// Close current if exists and initialize a new store.
        void reset(const std::string& store);

        /// Destruct the sink.  This must be done in order to close the store.
        ~Sink();

        /// Write one Bee object to the sink store.
        ///
        /// An object is written to an "event id" or index.  The sink maintains
        /// this index and will advance it if the object has a r/s/e that
        /// differs from previous or has an "type" that has been written to the
        /// current index.  This index is returned.
        size_t write(const Object& obj);

        /// Flush and close the store. Subsequent writes will fail.
        void close();

    private:

        void index(const Object& obj);

    };

}

#endif
