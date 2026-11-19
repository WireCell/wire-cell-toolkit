#ifndef WIRECELL_ITRAJECTORYSET
#define WIRECELL_ITRAJECTORYSET

#include "WireCellIface/ITrajectory.h"

namespace WireCell {

    /** An interface to the set of particle trajectories of one event -- the
     *  trajectory tree.  Navigate the tree via each ITrajectory's id()/parent().
     */
    class ITrajectorySet : public IData<ITrajectorySet> {
       public:
        virtual ~ITrajectorySet();

        /// An identifier unique to this set (typically the event number).
        virtual int ident() const = 0;

        /// The trajectories in this set.
        virtual ITrajectory::shared_vector trajectories() const = 0;

        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
