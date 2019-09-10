#ifndef WIRECELLUTIL_PIMPOS
#define WIRECELLUTIL_PIMPOS

#include "WireCellUtil/Binning.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

namespace WireCell {


    /** \brief Pitch-Impact-Position.

	A Pimpos object encapsulates information and methods related
	to geometry and binning of a plane of parallel and equidistant
	wires and a further uniform sub division along their pitch
	directions.

        A wire region or wire bin is the locus *centered* on a wire
        and which extends +/- 1/2 pitch from the wire position.

        An impact region or impact bin is the smaller locus of pitch
        formed by uniformly splitting a wire region.  It is *bound* by
        two (not centered on one) impact positions.

        Impact positions are "in phase" with wires" such that there is
        one impact position coincident with a wire position.

        For 3D "lab" to "plane" coordinate transformations there is a
        plane coordinate system which has these three orthogonal axes:

	axis0) anti nominal drift direction (normal to the wire plane) and
	axis1) along the wire direction such that,
	axis2) the pitch direction is the cross product of axis0 X axis1.
    */

    class Pimpos {
    public:
	/** Create a Pimpos object for a particular plane.
	    
            \param nwires is the number of wires in the plane.

            \param minwirepitch is the location in the pitch
            coordinate of the first (zero index) wire.

            \param maxwirepitch is the location in the pitch
            coordinate of the last wire (ie, index=nwires-1).

	    \param wire is a Vector which sets the direction of the
	    wires in the plane.  If the underlying wires are not
	    exactly parallel, this should be some representative
	    average wire direction.

	    \param pitch is a Vector which sets the direction and of
            the pitch of the wires in the plane.  The pitch should be
            such that the cross product, wire (x) pitch, points in the
            anti-drift direction for the region this plane services.
            If the underlying wires are not exactly parallel, this
            pitcvector should be some representative average pitch
            direction.

	    \param origin is a Point which sets an origin for all
	    transforms.  In particular, the projection of this origin
	    point along the drift direction to the plane of wires sets
	    the origin for the pitch coordinate.

	    \param nbins gives the number of of impact bins covering
	    one wire region.

	    The pitch extents and the origin vector must be expressed
	    in the WCT system of (length) units.
	 */
	Pimpos(int nwires, double minwirepitch, double maxwirepitch,
               const Vector& wire = Vector(0,1,0),
               const Vector& pitch = Vector(0,0,1),
               const Point& origin = Point(0,0,0),
	       int nimpact_bins_per_wire_region=10);

        /// Trivial accessor
        int nimpbins_per_wire() const { return m_nimpbins_per_wire; }
        void set_nimpbins_per_wire(int nbins);

        //// Geometry related:

        /// Return given 3-point origin for plane pitch.
	const Point& origin() const { return m_origin; }

        /// Return an axis of the plan.  0=normal to plane (aka
        /// anti-drift), 1=wire direction, 2=pitch direction.
	const Vector& axis(int i) const { return m_axis[i]; }


	/// Return the vector from the origin to the given point.
	Vector relative(const Point& pt) const;

	/// Return the distance from origin to point along the given
	/// axis.  Default is pitch distance.
	double distance(const Point& pt, int axis=2) const;

	/// Transform the given point into the Pimpos coordinate
	/// system.
	Point transform(const Point& pt) const;


        //// Binning related:

        /// Return wire region binning.  Each bin is centered on the
        /// wire.  Bin edges extend a distance of 1/2 pitch between
        /// neighboring wires.
        const Binning& region_binning() const { return m_regionbins; }

        /// Return impact position binning.  Bin edges are at impact
        /// positions and are in-phase with wire centers.
        const Binning& impact_binning() const { return m_impactbins; }

        /// Return a pair of indices.  The first is the index of the
        /// wire closest to the given pitch.  The second is the
        /// *relative* index of the impact closest impact position
        /// (index=0 means the impact position coincident with the
        /// wire).
        std::pair<int, int> closest(double pitch) const;

        /// Return the impact position index coincident with the wire index.
        int wire_impact(int wireind) const;

        /// Return the impact position indices (bin edges) at either
        /// extreme of the wire region.  The smaller index for this
        /// wire is the larger index of wireind-1's values and vice
        /// versa.
        std::pair<int,int> wire_impacts(int wireind) const;

        /// Return the impact position index which is the reflection
        /// of the given impact position index through the given wire
        /// index.
        int reflect(int wireind, int impind) const;

	

    private:
        int m_nwires;
        int m_nimpbins_per_wire;
	Point m_origin;
	Vector m_axis[3];
        Binning m_regionbins, m_impactbins;
    };

}


#endif
