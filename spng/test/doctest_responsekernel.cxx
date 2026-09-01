/** Tests for SPNG ResponseKernel normalization.
 *
 * The kernel produced by ResponseKernel samples a fixed continuous detector
 * response V(t), in units of ADC count per electron, at the configured sample
 * period T.  Consequently the discrete kernel R[m] = V(m*T) must satisfy
 *
 *     sum_m R[m] * T  ==  integral V(t) dt
 *
 * and the right hand side does not depend on T.  So "DC * period" is an
 * invariant of the configured period.  This is the property that makes charge
 * come out right after deconvolution: a charge Q deposited in one tick yields
 * a measured ADC series Q*R[m], and dividing by the response recovers Q only
 * when R is the true per-sample response.
 *
 * Note the field response is sampled at its own (fine) period, 100 ns for the
 * DUNE garfield response, while the kernel is normally requested at the 500 ns
 * ADC tick.  When the two periods are equal, TorchLMN::resample_interval()
 * short circuits, so a kernel built at the fine period exercises no resampling
 * and serves as the anchor for the comparison.
 */

#include "WireCellSpng/ResponseKernel.h"
#include "WireCellSpng/Testing.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::SPNG;

// The field response sampling period.  Must match the "period" of the field
// response file named below.
static const double fr_period = 100 * units::ns;

// The conventional ADC tick, and a second, unrelated coarse period.  Neither
// divides into the other in a way that could hide a period-dependent scale.
static const double adc_tick = 500 * units::ns;
static const double other_tick = 1000 * units::ns;

static const std::string response_file = "dune-garfield-1d565.json.bz2";

// FieldResponse lives in WireCellSigProc and ColdElecResponse in WireCellGen.
// These are loaded at run time by the test binary; SPNG itself does not link
// against either plugin.
static int init_plugins()
{
    auto& pm = PluginManager::instance();
    pm.add("WireCellAux");
    pm.add("WireCellGen");
    pm.add("WireCellSigProc");
    return 0;
}
static int plugins_initialized = init_plugins();

// Configure the shared FR and ER components exactly once.
static void configure_responses()
{
    static bool done = false;
    if (done) return;
    done = true;

    {
        auto icfg = Factory::lookup<IConfigurable>("FieldResponse");
        auto cfg = icfg->default_configuration();
        cfg["filename"] = response_file;
        icfg->configure(cfg);
    }
    {
        auto icfg = Factory::lookup<IConfigurable>("ColdElecResponse");
        icfg->configure(icfg->default_configuration());
    }
}

// Return sum_m R[m] for a kernel built at the given plane and period.  The DC
// element of the 2D spectrum is the sum over the whole (wire, time) kernel, and
// pad_waveform() zero pads, so this is independent of the requested shape as
// long as the shape is large enough to hold the response.
static double kernel_dc(int plane_index, double period)
{
    configure_responses();

    ResponseKernel rk;
    Configuration cfg = rk.default_configuration();
    cfg["field_response"] = "FieldResponse";
    cfg["elec_response"] = "ColdElecResponse";
    cfg["elec_duration"] = 20 * units::us;
    cfg["period"] = period;
    cfg["plane_index"] = plane_index;
    cfg["scale"] = 1.0;
    rk.configure(cfg);

    // Generous shape so no part of the response is cropped at any period.
    ITorchSpectrum::shape_t shape = {128, 8192};
    auto spec = rk.spectrum(shape);
    return torch::real(spec[0][0]).item<double>();
}

TEST_SUITE("spng response kernel") {

    // The normalization must not depend on the sample period the kernel is
    // built at.  This currently FAILS: ResponseKernel.cxx multiplies the
    // fine-grid convolution by the coarse period m_cfg.period instead of by
    // the field response period, so the kernel is too large by the ratio of
    // the two periods and "DC * period" grows linearly with the period.
    TEST_CASE("spng response kernel normalization is period independent") {

        for (int plane_index : {0, 1, 2}) {

            const double ref = kernel_dc(plane_index, fr_period) * fr_period;
            REQUIRE(std::abs(ref) > 0.0);

            for (double period : {adc_tick, other_tick}) {
                const double got = kernel_dc(plane_index, period) * period;
                const double ratio = got / ref;

                INFO("plane=" << plane_index
                     << " period=" << period / units::ns << "ns"
                     << " ref=" << ref << " got=" << got
                     << " ratio=" << ratio
                     << " (expected 1, a ratio of period/" << fr_period / units::ns
                     << "ns indicates the coarse period is being applied to a"
                        " fine-grid convolution)");

                // 2% covers the small difference in rational padding between
                // periods; the defect under test is a factor of 5 and 10.
                CHECK(std::abs(ratio - 1.0) < 0.02);
            }
        }
    }
}
