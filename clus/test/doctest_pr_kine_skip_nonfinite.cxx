// pr_ml doc 26 sec 7.3: the dQ/dx fit can leave a whole segment at dQ = NaN
// (S1 batch 2 PDVD cosmic1176: 13 / 13 fit points of a shower-trajectory
// segment).  segment_cal_kine_dQdx then sums NaN -- its clamps (dE < 0,
// dE > 50 MeV/cm) never catch it -- segment_cal_4mom builds a NaN 4-momentum
// and Aux::ParticleInfo throws "energy-momentum relation violated", ending the
// job.  Knob PR::set_kine_skip_nonfinite (TaggerCheckNeutrino
// pr_kine_skip_nonfinite, default false): a fit point with a non-finite dQ or
// dx is a missing measurement and is skipped, like an invalid one.
//
// Fail-first record: before the knob existed this file did not build (no
// set_kine_skip_nonfinite); the legacy case below pins the NaN.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;

namespace {
    // Box model with the uBooNE constants, as doctest_cal_kine_dqdx_zero_dx.cxx
    // (a clus test must not depend on the gen plugin).
    struct BoxModel : public IRecombinationModel {
        double m_efield{0.273}, m_a{0.93}, m_b{0.212}, m_rho{1.38}, m_wi{23.6e-6};
        virtual ~BoxModel() {}
        double operator()(double dE, double dX = 0.0) override
        {
            const double coeff = m_b / (m_efield * m_rho);
            return std::log(m_a + coeff * dE / dX * units::cm * units::cm / units::MeV) / (coeff * m_wi) * dX / units::cm;
        }
        double dE(double dQ, double dX) override
        {
            const double coeff = m_b / (m_efield * m_rho);
            const double a_exp = std::exp(dQ / dX * units::cm * coeff * m_wi);
            const double numerator = (a_exp - m_a) * units::MeV / units::cm * dX;
            return numerator / coeff;
        }
    };

    PR::SegmentPtr make_seg(const std::vector<double>& dq)
    {
        auto seg = std::make_shared<PR::Segment>();
        std::vector<PR::Fit> fits;
        for (size_t i = 0; i < dq.size(); ++i) {
            PR::Fit f;
            f.point = Point(0, 0, 0.6 * units::cm * i);
            f.dQ = dq[i];
            f.dx = 0.6 * units::cm;
            f.index = (int) i;
            fits.push_back(f);
        }
        seg->fits(fits);
        return seg;
    }

    struct KnobReset {
        ~KnobReset() { PR::set_kine_skip_nonfinite(false); }
    };
}

TEST_CASE("pr kine skip_nonfinite: legacy (off) sums a NaN dQ into a NaN kinetic energy")
{
    KnobReset reset;
    PR::set_kine_skip_nonfinite(false);
    IRecombinationModel::pointer model = std::make_shared<BoxModel>();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(std::isnan(PR::segment_cal_kine_dQdx(make_seg({30e3, nan, 30e3}), model)));
    CHECK(std::isnan(PR::segment_cal_kine_dQdx(make_seg({nan, nan, nan}), model)));
}

TEST_CASE("pr kine skip_nonfinite: on skips NaN points, all-NaN gives 0, finite input unchanged")
{
    KnobReset reset;
    IRecombinationModel::pointer model = std::make_shared<BoxModel>();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    PR::set_kine_skip_nonfinite(false);
    const double clean_off = PR::segment_cal_kine_dQdx(make_seg({30e3, 30e3, 30e3}), model);
    PR::set_kine_skip_nonfinite(true);
    const double clean_on = PR::segment_cal_kine_dQdx(make_seg({30e3, 30e3, 30e3}), model);
    CHECK(clean_on == clean_off);                         // bit for bit without a NaN
    const double mixed = PR::segment_cal_kine_dQdx(make_seg({30e3, nan, 30e3}), model);
    CHECK(std::isfinite(mixed));
    CHECK(mixed > 0);
    CHECK(mixed < clean_on);                              // the NaN point contributes nothing
    CHECK(PR::segment_cal_kine_dQdx(make_seg({nan, nan, nan}), model) == 0.0);
}
