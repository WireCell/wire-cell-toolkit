#ifndef WIRECELL_GEN_DRIFTER
#define WIRECELL_GEN_DRIFTER

#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRandom.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"

#include <set>

namespace WireCell {

    // The xregions are a list of objects specifying locations on the
    //  X (drift) axis of the coordinate system in which depo
    //  positions are defined.  Three locations may be defined.
    //
    //

    namespace Gen {

        /** This component drifts depos bounded by planes
         * perpendicular to the X-axis.  The boundary planes are
         * specified with the "xregions" list.  Each list is an object
         * fully specified with a "cathode" an "anode" and a
         * "response" attribute giving X locations in the same
         * coordinate system as depos of three planes.
         *
         * - cathode :: a plane which bounds the maximum possible drift.
         * - anode :: a plane which bounds the minimum possible drift.
         * - response :: a plane to which all depositions are drifted.
         *
         * If "anode" is not given then its value is take to be that of
         * "response" and vice versa and at least one must be specified.
         * A "cathode" value must be specified.
         *
         * Any depo not falling between "anode" and "cathode" will be
         * dropped.
         *
         * Any depo falling between "response" and "cathode" will be
         * drifted to the "response" plane.
         *
         * Any depo falling between "anode" and "response" will be
         * ANTI-DRIFTED to the "response" plane.  Ie, it will be
         * "BACKED UP" in space an time as if it had be produced
         * earlier and at the response plane.
         *
         * Input depositions must be ordered in absolute time (their
         * current time) and output depositions are produced ordered
         * by their time after being drifted to the response plane.
         *
         * Diffusion and absorption effects and also, optionally,
         * fluctuations are applied.  Fano factor and Recombination
         * are not applied in this component (see IRecombinationModel
         * implementations).
         *
         * Typically a drifter is used just prior to a ductor and in
         * such cases the "response" plane should be made coincident
         * with the non-physical "response plane" which defines the
         * starting point for the field response functions.  The
         * location of the response plane *realtive* to the wire
         * planes can be found using:
         *
         * $ wriecell-sigproc response-info garfield-1d-3planes-21wires-6impacts-dune-v1.json.bz2
         * origin:10.00 cm, period:0.10 us, tstart:0.00 us, speed:1.60 mm/us, axis:(1.00,0.00,0.00)
         *    plane:0, location:9.4200mm, pitch:4.7100mm
         *    plane:1, location:4.7100mm, pitch:4.7100mm
         *    plane:2, location:0.0000mm, pitch:4.7100mm
         *
         * Here, "origin" gives the location of the response plane.
         * The location of the wire planes according to wire geometry
         * can be similarly dumped.
         *
         * $ wirecell-util wires-info protodune-wires-larsoft-v3.json.bz2
         * anode:0 face:0 X=[-3584.63,-3584.63]mm Y=[6066.70,6066.70]mm Z=[7.92,7.92]mm
         *     0: x=-3584.63mm dx=9.5250mm
         *     1: x=-3589.39mm dx=4.7620mm
         *     2: x=-3594.16mm dx=0.0000mm
         * ....
         * anode:5 face:1 X=[3584.63,3584.63]mm Y=[6066.70,6066.70]mm Z=[6940.01,6940.01]mm
         *     0: x=3584.63mm dx=-9.5250mm
         *     1: x=3589.39mm dx=-4.7620mm
         *     2: x=3594.16mm dx=0.0000mm
         *
         * Note, as can see, these two sources of information may not
         * be consistent w.r.t. the inter-plane separation distance
         * (4.71mm and 4.76mm, respectively).  This mismatch will
         * result in a relative shift in time between the planes for
         * various waveform features (eg induction zero crossings and
         * collection peak).
         *
         * For the example above, likely candidates for "anode" X
         * locations are:
         *
         *    x = -3594.16mm + 10cm
         *
         * and
         *
         *    x = +3594.16mm - 10cm
         */
        class Drifter : public Aux::Logger,
                        public IDrifter, public IConfigurable {
           public:
            Drifter();
            virtual ~Drifter();

            virtual void reset();
            virtual bool operator()(const input_pointer& depo, output_queue& outq);

            /// WireCell::IConfigurable interface.
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            // Implementation methods.

            // Do actual transport, producing a new depo
            IDepo::pointer transport(IDepo::pointer depo);

            // Return the "proper time" for a deposition
            double proper_time(IDepo::pointer depo);

            bool insert(const input_pointer& depo);
            void flush(output_queue& outq);
            void flush_ripe(output_queue& outq, double now);

            // Reset lifetime e.g. based on a larsoft database.
            // Detailed implementation in a subclass.
            virtual void set_lifetime(double lifetime_to_set) { m_lifetime = lifetime_to_set; };

           private:
            IRandom::pointer m_rng;
            std::string m_rng_tn;

            // Longitudinal and Transverse coefficients of diffusion
            // in units of [length^2]/[time].
            double m_DL, m_DT;

            // Electron absorption lifetime.
            double m_lifetime;

            // If true, fluctuate by number of absorbed electrons.
            bool m_fluctuate;

            double m_speed;    // drift speeds
            double m_toffset;  // time offset

	    double m_scale_factor;

            int n_dropped, n_drifted;

            // keep the depos sorted by time
            struct DepoTimeCompare {
                bool operator()(const IDepo::pointer& lhs, const IDepo::pointer& rhs) const;
            };

            // A little helper to carry the region extent and depo buffers.
            struct Xregion {
                Xregion(Configuration cfg);
                double anode, response, cathode;
                typedef std::set<IDepo::pointer, DepoTimeCompare> ordered_depos_t;
                ordered_depos_t depos;  // buffer depos

                bool inside_bulk(double x) const;
                bool inside_response(double x) const;
            };
            std::vector<Xregion> m_xregions;

            struct IsInsideBulk {
                const input_pointer& depo;
                IsInsideBulk(const input_pointer& depo)
                  : depo(depo)
                {
                }
                bool operator()(const Xregion& xr) const { return xr.inside_bulk(depo->pos().x()); }
            };
            struct IsInsideResp {
                const input_pointer& depo;
                IsInsideResp(const input_pointer& depo)
                  : depo(depo)
                {
                }
                bool operator()(const Xregion& xr) const { return xr.inside_response(depo->pos().x()); }
            };

        };  // Drifter

    }  // namespace Gen

}  // namespace WireCell

#endif
