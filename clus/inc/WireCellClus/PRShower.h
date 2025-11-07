#ifndef WIRECELL_CLUS_PR_SHOWER
#define WIRECELL_CLUS_PR_SHOWER

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRTrajectoryView.h"

#include "WireCellUtil/Flagged.h"
#include "WireCellUtil/Point.h"

namespace WireCell::Clus::PR {

    /** The "flags" that may be set on a shower.
     */
    enum class ShowerFlags {
        /// The segment has no particular category.
        kUndefined = 0,
        /// The shower flag
        kShower = 1<<1,
        /// The kinematics flag
        kKinematics = 1<<2,
    };

    /** The data attributes of a PR::Shower
        
        The original WCShower has a huge number of attributes that are merely
        carried by the shower.  This struct factors out that data to facilitate
        writing stand-alone functions.

        Anything that is not part of the core shower-as-graph-view concept gets
        lumped into this ShowerData struct.

        Note, "flags" are set via Shower's Flagged base class.

    */
    struct ShowerData
    {
        int particle_type;
        double kenergy_range;
        double kenergy_dQdx;
        double kenergy_charge;
        double kenergy_best;
    
        WireCell::Point start_point;
        WireCell::Point end_point;
        WireCell::Vector init_dir;

        /// 1 for direct connection, 2 for indirection connection with a gap, 3
        /// for associations (not clear if this should be connected or not
        int start_connection_type;
    };


    /** Model a shower-like view of a trajectory.

        This is the WCT equivalent to a WCT WCShower.
     */
    class Shower
        : public TrajectoryView
        , public Flagged<ShowerFlags> // can set flags
        , public HasDPCs<Shower>      // has associated DynamicPointClouds.
    {
    public:

        Shower(Graph& graph);
        virtual ~Shower();

        // The bag of attributes is directly exposed to user.
        ShowerData data;


        // Getters

        /** Get the vertex that is considered the start of the shower.
         */
        VertexPtr start_vertex();

        /** Get the segment that is considered the start of the shower.
         */
        SegmentPtr start_segment();

        // Chainable setters

        /** Chainable setter of start vertex.

            The vertex must already be added to the underlying graph that this
            shower views.

            The vertex will be added to the view.

            The vertex will replace any existing start vertex and not remove the
            prior vertex from the shower's view.  Use `Shower::remove_vertex()`
            to explicitly remove from the view and `PR::remove_vertex()` to
            remove it from the underlying graph, if either are required.

            If the vertex lacks a valid descriptor, eg has yet to be added to
            the underlying graph, this function is a no-op and the stored
            start_vertex is nullified.
        */
        Shower& start_vertex(VertexPtr vtx);

        /** Chainable setter of start segment.

            This has the same semantics and caveats as the chainable setter:
            `start_vertex(VertexPtr)`.
        */
        Shower& start_segment(SegmentPtr seg);

    private:

        VertexPtr m_start_vertex;
        SegmentPtr m_start_segment;

    };

}
#endif
