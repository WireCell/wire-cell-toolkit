/** Digitizer converts voltage waveforms to integer ADC ones.
 *
 * Resulting waveforms are still in floating-point form and should be
 * round()'ed and truncated to whatever integer representation is
 * wanted by some subsequent node.
 */

#ifndef WIRECELL_DIGITIZER
#define WIRECELL_DIGITIZER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"
#include <deque>

namespace WireCell {

    namespace Gen {

        class Digitizer : public Aux::Logger,
                          public IFrameFilter, public IConfigurable {
           public:
            Digitizer(const std::string& anode_tn = "AnodePlane",
                      int resolution = 12,  // bits of resolution
                      double gain = -1.0,   // input gain
                      std::vector<double> fullscale = {0.0, 2.0 * units::volt},
                      std::vector<double> baselines = {900 * units::mV, 900 * units::mV, 200 * units::mV});
            virtual ~Digitizer();

            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            // Voltage frame goes in, ADC frame comes out.
            virtual bool operator()(const input_pointer& inframe, output_pointer& outframe);

            // implementation method.  Return a "floating point ADC"
            // value for the given voltage assumed to have been lifted
            // to the baseline.
            double digitize(double voltage);

           private:
            IAnodePlane::pointer m_anode;
            std::string m_anode_tn;
            int m_resolution;
            double m_gain;
            std::vector<double> m_fullscale, m_baselines;
            std::string m_frame_tag;
            size_t m_count{0};
            bool m_round{false};  // apply round to FP ADC values
        };

    }  // namespace Gen
}  // namespace WireCell
#endif
