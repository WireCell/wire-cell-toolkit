#ifndef WIRECELL_ITRAJECTORY
#define WIRECELL_ITRAJECTORY

#include "WireCellIface/IData.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    /** An interface to one particle trajectory from a tracking simulation.
     *
     *  Trajectories form a tree: parent() gives the id() of the trajectory
     *  that spawned this one (negative for a primary), so the tree is
     *  reconstructed from the (id, parent) relation across an ITrajectorySet.
     *
     *  The interface models the tree-essential core (ids, PDG, initial point
     *  and initial momentum), balancing edep-sim's TG4Trajectory and Geant4's
     *  G4VTrajectory.  The particle name, the detailed trajectory points
     *  (position/time/momentum/process along the path) and other less-common
     *  attributes are carried, when provided, in metadata():
     *
     *    "name"   (string)  particle name
     *    "points" (array)   optional detailed trajectory points
     *
     *  Values are expressed in the WCT system of units.
     */
    class ITrajectory : public IData<ITrajectory> {
       public:
        virtual ~ITrajectory();

        /// The track id of this trajectory (its index within the set).
        virtual int id() const = 0;

        /// The track id of the parent trajectory, negative for a primary.
        virtual int parent() const = 0;

        /// The PDG particle code.
        virtual int pdg() const = 0;

        /// The initial position of the trajectory.
        virtual Point start() const = 0;

        /// The time at start().
        virtual double start_time() const = 0;

        /// The initial 3-momentum.
        virtual Vector momentum() const = 0;

        /// The initial total energy.
        virtual double energy() const = 0;

        /// Less-common attributes (name, detailed points, ...); empty if none.
        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
