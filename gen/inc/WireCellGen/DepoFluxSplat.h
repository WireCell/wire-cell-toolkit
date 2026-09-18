/** DepoFluxSplat is an approximate combination DepoTransform and OmnibusSigProc.

    DepoFluxSplat transforms an IDepoSet into an IFrame that represents a "true
    signal" comparable to the output of signal processing.  It does this by
    adding configurable amount of Gaussian smearing to the drifted depos to
    account for the smearing accrued by the convolution/deconvolution cycle of
    sim+sigproc and then bins the resulting distribution onto channels via their
    wires.

    DepoFluxSplat operates on the context of a single anode and will write
    either a sparse or dense frame.  When the frame is sparse, each depo results
    in a family of traces that will in general overlap in tick and channel with
    the family of traces from nearby depos.

    Differences between DepoFluxSplat and DepoSplat:
    - DepoFluxSplat is general to all detectors (no hard-wired magic numbers)
    - DepoFluxSplat consumes IDepoSet instead of IDepo

    Differences between DepoFluxSplat and DepoFluxWriter (from larwirecell)
    - DepoFluxSplat is "pure" Wire-Cell, no interace with art::Event.
    - DepoFluxSplat outputs IFrame instead of SimChannel
    - DepoFluxSplat operates in the context of a single AnodePlane.
    - DepoFluxSplat drops the misguided handling of "energy".

    DepoFluxSplat internals are largely copy-pasted from DepoFluxWriter.  It is
    hoped that DepoFluxWriter will be refactored to replace its internals with
    instances of DepoFluxSplat.

*/
#ifndef WIRECELLGEN_DEPOFLUXSPLAT
#define WIRECELLGEN_DEPOFLUXSPLAT

#include "WireCellAux/Logger.h"
#include "WireCellIface/IDepoFramer.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"

namespace WireCell::Gen {

    class DepoFluxSplat : public Aux::Logger,
                          public IDepoFramer,
                          public IConfigurable
    {
      public:
        DepoFluxSplat();
        ~DepoFluxSplat();

        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual void configure(const WireCell::Configuration& cfg);
        virtual Configuration default_configuration() const;

      private:

        /// Configuration:
        ///
        /// anode - the type/name for the IAnodePlane.
        ///
        IAnodePlane::pointer m_anode;

        /// field_response - name of an IFieldResponse from which drift_speed
        /// (speed) and response_plane (origin) is taken.  These two parameters
        /// may be explicitly overridden (see next).  If both are overridden
        /// then field_response parameter may be left unspecified.

        /// drift_speed, response_plane - give explicit drift speed and
        /// response plane distance.  Default is unspecified and the values
        /// provided by the field response are used.
        double m_speed{0}, m_origin{0};

        /// tick - the sample time over which to integrate depo flux
        /// into time bins.
        ///
        /// window_start - the start of an acceptance window in NOMINAL
        /// TIME.
        ///
        /// window_duration - the duration of the acceptance window.
        Binning m_tbins;

        /// sparse - bool, if true save a sparse IFrame, else dense.  Default is
        /// false.  Sparse frames can be very large an only useful for special
        /// studies that require the many-to-many mapping between depos and
        /// trace contributions.
        bool m_sparse{false};

        /// nsigma - number of sigma at which to truncate a depo Gaussian
        double m_nsigma{3.0};

        /// referenced_time - An arbitrary time that is SUBTRACTED from the
        /// NOMINAL TIME of the window start when setting the output frame time.
        double m_reftime{0};

        /// time_offsets - An arbitrary time that is ADDED to NOMINAL TIMES for
        /// each plane.  If set, this offset is reflected in the "tbin" value of
        /// traces for a given plane.
        std::vector<int> m_tick_offsets;


	std::vector<int> m_process_planes {0,1,2};

        /// smear_tran - Extra smearing applied to the depo in the
        /// transverse direction of each plane.  This given in units of
        /// pitch (ie, smear_tran=2.0 would additionally smear over 2
        /// pitch distances) and it is added in quadrature with the
        /// Gaussian sigma of the depo.  If zero (default) no extra
        /// smearing is applied.
        ///
        /// A single scalar number may be provided which will add a common
        /// smearing for all planes.  Or, a per-plane smearing may be provided.
        std::vector<double> m_smear_tran;

        /// smear_long - Extra smearing applied to the depo in the
        /// longitudinal direction.  This given in units of tick (ie,
        /// smear_long=2.0 would additionally smear across two ticks)
        /// and added in quadrature with the Gaussian sigma of the
        /// depo.  If zero (default) no extra smearing is applied.
        ///
        /// A single scalar number may be provided which will add a common
        /// smearing for all planes.  Or, a per-plane smearing may be provided.
        std::vector<double> m_smear_long;


        /// trio_file - if non-empty, the name of a .npz file to which the
        /// depo-to-channel association ("trios") is also written.  Empty
        /// (default) disables this entirely and the component behaves exactly
        /// as before.
        ///
        /// Each depo lands on one wire in each of the three planes, and those
        /// three channels at their ticks are a "trio": a cross-plane
        /// correspondence known to be real because one deposition made it.
        /// The frame this component produces sums every depo's contribution
        /// and so discards which depo went where; recovering that downstream
        /// would mean replaying Drifter (whose lifetime attenuation is
        /// stochastic) and this component.  Emitting it here is exact.
        ///
        /// Two arrays are appended per call, named with the frame ident:
        ///   trio_index_<ident> (M,6) int32   tbin_{u,v,w}, chan_{u,v,w}
        ///   trio_value_<ident> (M,3) double  charge, sigma_long, sigma_tran
        /// Depos landing on the same six indices are summed into one row; the
        /// sigmas are charge-weighted means over those depos.  Charge and
        /// sigmas are post-drift, as seen by this component.
        std::string m_trio_file{""};

        /// trio_min_charge - drop aggregated trios whose absolute summed
        /// charge is below this.  Default 0 keeps everything.
        double m_trio_min_charge{0};

        // internal, return nullptr if depo is not in anode plane.  o.w. return
        // anode face.
        IAnodeFace::pointer find_face(const WireCell::IDepo::pointer& depo);

        // Count executions for debug log.
        size_t m_count{0};

    };

}

#endif 
