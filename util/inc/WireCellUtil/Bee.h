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
#include "WireCellUtil/Units.h"

#include <boost/filesystem.hpp>

#include <unordered_set>

namespace WireCell::Bee {

    // Bee supports more than one type (schema) of JSON object.  The Object base
    // class handles the commonalities.
    class Object {

    protected: 

        // Note: exceptionally for Wire-Cell toolkit, any distance values
        // (x,y,z) held in m_data here are in Bee units (cm), not WCT units.
        Configuration m_data;

        // Subclass MUST set this for sink to be meaningful.
        std::string m_name{""};

        explicit Object(const std::string& name = "",
                        Configuration data = Json::nullValue)
            : m_data(data), m_name(name) {}

    public:

        /// Return self as a JSON string
        std::string json() const;

        /// Return digest of contents.
        size_t hash() const;

        /// The name must be suitable for use as part of a file name.  It
        /// identifies a Bee class of data (eg "clusters").
        std::string name() const { return m_name; }

        // Method to get JSON data
        Configuration& data() { return m_data; }
        const Configuration& data() const { return m_data; }
    };

    /// A Bee "object" represents a set of 3D points with attributes.  It maps
    /// to one Bee JSON file.  The "type" is used as the Object name.
    class Points : public Object {

    public:

        /// Default constructor sets geom and type to empty.  This will not make
        /// Bee happy but it is okay if this Points is an intermediate object
        /// that will later be .append()'d to another which provide geom and
        /// type.
        Points();

        /// Construct with metadata
        Points(const std::string& geom, // canonical detector geometry name, eg "uboone"
               const std::string& type, // name of the "algorithm"
               int run=0, int sub=0, int evt=0);


        /// Drop arrays and reset the e/s/r numbers.
        /// CAUTION: opposite order than may be expected.
        /// If sub or run are negative, any previous value is kept.
        void reset(int evt, int sub=-1, int run=-1);

        void detector(const std::string& geom);
        void algorithm(const std::string& type);
        std::string algorithm() const;

        /// Set the run/subrun/event numbers.  Pass negative to leave current
        /// value unchanged.
        void rse(int run, int sub, int evt);

        /// Return run/sub/evt as vector of int.
        std::vector<int> rse() const;

        /// Add a WCT point.  Note, p is expected to be in usual WCT system-of-units.
        void append(const Point& p, double q=0, int clid=0, int real_clid=0);

        /// Simply append obj's x,y,z,q,cluster_id arrays to this.
        void append(const Points& obj);

        size_t size() const;
        bool empty() const;

        // return the cluster id of the last appended point
        int back_cluster_id() const;

    };

    /// Represent a set of 2D areas or regions.  Each patch is described by an
    /// ordered set of 2D corner points.  The Patches will do the ordering.
    class Patches : public Object {
        std::vector<double> m_y, m_z; // a buffered patch in the making
        double m_tolerance{0};
        size_t m_minpts{3};        // at least a triangle
    public:

        /// Create a named patches.  A point is ignored if it is withing
        /// tolerance in y or z direction of a previously appended point.
        /// Tolerance must be provided in WCT system of units.  For a patch to
        /// be considered it must have at least minpts number of points after
        /// tolerance filtering is applied.
        explicit Patches(const std::string& name, double tolerance=0*units::mm,
                         size_t minpts=3);
        
        /// Append a single point to the current patch.
        void append(double y, double z);

        /// Flush buffered patch, saving it to the JSON data.
        void flush();

        /// Append all points of one patch and flush.
        template<typename It>
        void append(It ybeg, It yend, It zbeg, It zend) {
            m_y.insert(m_y.end(), ybeg, yend);
            m_z.insert(m_z.end(), zbeg, zend);
            flush();
        }

        // Number of patches stored 
        size_t size() const;

        // True if no patches are stored (note, could still have points buffered awaiting a close()
        bool empty() const;

        // Clear any stored JSON and buffered patch
        void clear();


    };


    // todo: source

    /// A Bee Sink persists Objects to some store.
    class Sink {

        boost::iostreams::filtering_ostream m_out;
        size_t m_index{0};
        std::unordered_set<std::string> m_seen; // names at current index

        int m_runNo{0};
        int m_subRunNo{0};
        int m_eventNo{0};
    public:

        /// Construct a sink of Bee objects.  The "store" indicates some
        /// persistent resource to receive the objects.  This is typically a
        /// file name ending in ".zip".
        Sink();
        Sink(const std::string& store);

        /// Construct with store path and initial index
        Sink(const std::string& store, size_t initial_index);

        /// Close current if exists and initialize a new store.
        void reset(const std::string& store);

        /// Reset store path with initial index
        void reset(const std::string& store, size_t initial_index);

        /// Set index directly
        void set_index(size_t index);
        
        /// Get current index
        size_t get_index() const;

        /// Destruct the sink.  This must be done in order to close the store.
        ~Sink();

        // Add method to set RSE
        void set_rse(int run, int subrun, int event) {
            m_runNo = run;
            m_subRunNo = subrun;
            m_eventNo = event;
        }    
        
    

        /// Write one Bee objects to the sink store.
        ///
        /// Bee objects are stored along two axis: "index" and "name".  Only one
        /// object of a given name can be stored in a given index.  The "index"
        /// increments monotonically from 0.  When a write() sees an object name
        /// on the existing index, the index is incremented.
        size_t write(const Object& obj);

        /// Increment to a new index.
        void flush();

        /// Close the store. Subsequent writes will fail.
        void close();

        /// Return path in store for an object at the current index.
        std::string store_path(const Object& obj) const;

    private:

        void index(const Object& obj);


    };

}

#endif
