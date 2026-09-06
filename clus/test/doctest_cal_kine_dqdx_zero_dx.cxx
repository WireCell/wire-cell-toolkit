// doc pdvd/45 sec 5.4 follow-up: the vector form of cal_kine_dQdx sums
// recomb_model->dE(dQ, dx) over every fit point of a muon chain.  With dx == 0
// (a coincident pair of fit points) the box model computes dQ/dx = 0/0 and the
// NaN survives the sum into kenergy_dQdx -> kine_reco_Enu.  The prototype
// (ProtoSegment.cxx:1316-1339) divides by (dx + 1e-9) and then multiplies by
// dx, so a zero-length point contributes exactly 0; the port dropped the
// epsilon.  10 PDVD candidates hit this in production, 72 with excl_t0_frame on.
//
// Fail-first record (build of 2026-09-05): the legacy call returned -nan for the
// three-point vector below (CHECK(std::isfinite(e)) failed, -nan vs Approx(3.59883)).

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;

namespace {
    // Mirror of Gen::PracticalBoxRecombination::dE with the uBooNE constants
    // (gen/src/PracticalRecombinationModels.cxx:43-51).  Kept local because a
    // clus test must not depend on the gen plugin.
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
}

TEST_CASE("cal_kine_dQdx: legacy path pins the dx == 0 NaN (what the knob guards)")
{
    IRecombinationModel::pointer model = std::make_shared<BoxModel>();
    std::vector<double> dQ{50e3, 0.0, 50e3};
    std::vector<double> dx{1.0 * units::cm, 0.0, 1.0 * units::cm};
    // Default (skip_zero_dx = false) is the byte-identical production path.
    CHECK(std::isnan(Clus::PR::cal_kine_dQdx(dQ, dx, model)));
    // A positive-charge zero-length point is NaN too (exp(inf) * 0).
    std::vector<double> dQb{50e3, 50e3};
    std::vector<double> dxb{1.0 * units::cm, 0.0};
    CHECK(std::isnan(Clus::PR::cal_kine_dQdx(dQb, dxb, model)));
}

TEST_CASE("cal_kine_dQdx: skip_zero_dx drops dx <= 0 points and leaves the rest untouched")
{
    IRecombinationModel::pointer model = std::make_shared<BoxModel>();
    std::vector<double> dQ{50e3, 0.0, 50e3, 40e3};
    std::vector<double> dx{1.0 * units::cm, 0.0, 1.0 * units::cm, -0.5 * units::cm};
    const double e = Clus::PR::cal_kine_dQdx(dQ, dx, model, true);
    CHECK(std::isfinite(e));
    std::vector<double> dQ2{50e3, 50e3};
    std::vector<double> dx2{1.0 * units::cm, 1.0 * units::cm};
    CHECK(e == doctest::Approx(Clus::PR::cal_kine_dQdx(dQ2, dx2, model)));
    // Knob on with no degenerate point == legacy, bit for bit.
    CHECK(Clus::PR::cal_kine_dQdx(dQ2, dx2, model, true) == Clus::PR::cal_kine_dQdx(dQ2, dx2, model, false));
}
