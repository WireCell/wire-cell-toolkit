/** Aux::fill must write "child" as the row of the depo whose prior()
 *  the row is, so that a reader can rebuild the prior chain.
 *
 *  Defect (found 2026-10-21, pr_ml doc 19 sec 2.1): push_depo passed
 *  dp.size() taken after the emplace, so every gen>0 row pointed at
 *  itself and DepoFileSource made each prior its own prior.  After a
 *  DepoFileSink -> DepoFileSource round trip of drifted depos the
 *  active depos came back with prior()==nullptr and no warning.
 *
 *  Layout checked here (rows are pushed depth-first):
 *    row 0: depo A   gen 0  child 0
 *    row 1: A.prior  gen 1  child 0   (its child is row 0)
 *    row 2: A.prior.prior gen 2 child 1
 *    row 3: depo B   gen 0  child 0
 *    row 4: B.prior  gen 1  child 3
 */

#include "WireCellAux/DepoTools.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_CASE("depo fill child index points at the child row")
{
    auto a2 = std::make_shared<Aux::SimpleDepo>(2.0, Point(0, 0, 0), 1.0, nullptr, 0, 0, 7, 13, 0.5);
    auto a1 = std::make_shared<Aux::SimpleDepo>(1.0, Point(1, 0, 0), 1.0, a2, 0, 0, 7, 0, 1.0);
    auto a0 = std::make_shared<Aux::SimpleDepo>(0.0, Point(2, 0, 0), 1.0, a1, 0, 0, 7, 0, 1.0);
    auto b1 = std::make_shared<Aux::SimpleDepo>(3.0, Point(0, 1, 0), 2.0, nullptr, 0, 0, 8, 2212, 3.0);
    auto b0 = std::make_shared<Aux::SimpleDepo>(4.0, Point(0, 2, 0), 2.0, b1, 0, 0, 8, 0, 1.0);

    IDepo::vector depos = {a0, b0};

    Array::array_xxf data;
    Array::array_xxi info;
    Aux::fill(data, info, depos, true);

    REQUIRE(data.rows() == 5);
    REQUIRE(info.rows() == 5);
    CHECK(data.cols() == 7);
    CHECK(info.cols() == 4);

    // gen column
    CHECK(info(0, 2) == 0);
    CHECK(info(1, 2) == 1);
    CHECK(info(2, 2) == 2);
    CHECK(info(3, 2) == 0);
    CHECK(info(4, 2) == 1);

    // child column: the row whose prior() this row is
    CHECK(info(1, 3) == 0);
    CHECK(info(2, 3) == 1);
    CHECK(info(4, 3) == 3);

    // no prior row may point at itself
    for (int r = 0; r < info.rows(); ++r) {
        if (info(r, 2) > 0) {
            CHECK(info(r, 3) != r);
        }
    }

    // ids and pdg carried per row
    CHECK(info(2, 0) == 7);
    CHECK(info(2, 1) == 13);
    CHECK(info(4, 1) == 2212);
    CHECK(data(2, 0) == doctest::Approx(2.0));

    // truth sidecar is row-aligned and carries energy
    Eigen::ArrayXXd truth;
    Aux::fill_truth(truth, depos, true);
    REQUIRE(truth.rows() == 5);
    CHECK(truth.cols() == 3);
    CHECK(truth(2, 0) == doctest::Approx(0.5));
    CHECK(truth(4, 0) == doctest::Approx(3.0));
    CHECK(truth(0, 0) == doctest::Approx(1.0));
    CHECK(truth(0, 1) == 0.0);

    // without priors: active rows only, gen 0, child 0
    Aux::fill(data, info, depos, false);
    REQUIRE(info.rows() == 2);
    CHECK(info(0, 2) == 0);
    CHECK(info(1, 3) == 0);
}
