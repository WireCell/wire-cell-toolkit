/** An RCResponse waveform is simple resistor-capacitor response.
 */

#ifndef WIRECELL_GEN_RCRESPONSE
#define WIRECELL_GEN_RCRESPONSE

#include "WireCellIface/IWaveform.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Response.h"

namespace WireCell::Gen {

    class RCResponse : public IWaveform, public IConfigurable {
      public:
        RCResponse(int nticks = 10000, double t0 = 0, double width = 1.0 * units::ms,
                   double tick = 0.5 * units::us);

        // IConfigurable interface
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // The starting point of the sampling
        virtual double waveform_start() const;
        // The sampling period aka bin width
        virtual double waveform_period() const;
        // The collection of samples
        virtual const sequence_type& waveform_samples() const;
        // The collection of samples rebinned
        virtual sequence_type waveform_samples(const WireCell::Binning& tbins) const;

      private:
        std::unique_ptr<Response::SimpleRC> m_rc;
        Configuration m_cfg;
        sequence_type m_wave;
    };

}  // namespace WireCell::Gen

#endif
