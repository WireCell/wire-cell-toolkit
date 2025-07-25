#ifndef WIRECELL_CLUS_PR_COMMON
#define WIRECELL_CLUS_PR_COMMON

#include "WireCellUtil/Point.h"

namespace WireCell::Clus::PR {

    /// A WCPoint is a 3D point and corresponding wire indices and an index.
    //
    // FIXME: does this need any change given we now support wrapped wires and
    // multi APA/face detectors?
    struct WCPoint {
        WireCell::Point point;   // 3D point 
        int uvw[3] = {-1,-1,-1}; // wire indices
        int index{-1};           // point index in some container

        // FIXME: WCP had this, does WCT need it?
        // blob* b;

        
        // Return true if the point information has been filled.
        bool valid() const {
            if (index < 0) return false;
            return true;
        }
    };

    /** A Fit holds information predicted about a point by some "fit".
     *
     * A PR::Vertex has a scalar Fit object and PR::Segment has a vector<Fit>
     *
     * Note, WCP's ProtoSegment had struct-of-array instead of vector<Fit>
     */
    struct Fit {
        WireCell::Point point;
        double dQ{-1}, dx{0}, pu{-1}, pv{-1}, pw{-1}, pt{0}, reduced_chi2{-1};

        int index{-1};
        double range{-1};
        bool flag_fix{false};

        // Explicitly NOT defined:

        // bool flag_fit.  This seems never actually used in WCP.  If needed,
        // can we simply test on index or range?

        // Restore values to invalid
        void reset() {
            index = -1;
            flag_fix = false;
            range = -1;
        }

        double distance(const Point& p) {
            return (p - point).magnitude();
        }

        /** Return true if fit information has been filled */
        bool valid() {
            if (index < 0 || range < 0) return false;
            return true;
        }
    };

    
    /** Some mixin classes, eg used by Vertex and Segment.

        See also "Graphed" from PRGraph.h and "Flagged" from util.
     */

    /** A class that mixes in Identities caries its own ID number and that of an
       associated cluster.  A negative ID is considered invalid/undefined.
    */
    class Identities {
    public:

        /// Getters

        /// The ID number for this vertex.
        int ident() const { return m_ident; }

        /// The ID of the associated cluster
        int cluster_id() const { return m_cluster_id; }
         
        /// Chainable setters
        Identities& ident(int id) { m_ident = id; return *this; }
        Identities& cluster_id(int cid) { m_cluster_id = cid; return *this; }

    private:

        int m_ident{-1};
        int m_cluster_id{-1};

    };

}


#endif
