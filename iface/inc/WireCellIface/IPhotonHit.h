#ifndef WIRECELL_IPHOTONHIT
#define WIRECELL_IPHOTONHIT

#include "WireCellIface/IData.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    /** An interface to one optical-photon hit from a tracking simulation.
     *
     *  Models edep-sim's TG4PhotonHit: a photon deposited energy in a photon
     *  sensitive detector, tracked from its start() (emission) to its stop()
     *  (detection).  The sensitive-detector name, the emission process and the
     *  wavelength are carried in metadata():
     *
     *    "sensitive_detector" (string)
     *    "process"            (int)
     *    "wavelength"         (double)
     *
     *  Values are expressed in the WCT system of units.
     */
    class IPhotonHit : public IData<IPhotonHit> {
       public:
        virtual ~IPhotonHit();

        /// Where the photon started (emission).
        virtual Point start() const = 0;

        /// Where the photon stopped (detection).
        virtual Point stop() const = 0;

        /// The time at start().
        virtual double start_time() const = 0;

        /// The time at stop().
        virtual double stop_time() const = 0;

        /// The deposited photon energy.
        virtual double energy() const = 0;

        /// The track id of the primary/particle that produced the photon.
        virtual int id() const = 0;

        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
