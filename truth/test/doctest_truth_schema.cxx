/** The truth schema is a contract: names are unique, orders are the
 *  documented ones, and the SBND attacher columns of truth_pf stay
 *  first and in order so the existing reader
 *  (root/src/SbndPrMagnifyTrackingVisitor.cxx write_truth_table)
 *  keeps working.
 */

#include "WireCellTruth/Schema.h"
#include "WireCellUtil/doctest.h"

#include <set>

using namespace WireCell;

namespace {
    bool unique(const std::vector<std::string>& v)
    {
        return std::set<std::string>(v.begin(), v.end()).size() == v.size();
    }
}

TEST_CASE("truth schema columns are unique and ordered")
{
    CHECK(Truth::schema_version == 1);
    CHECK(Truth::schema_version_string() == "wct-truth-schema-1");

    CHECK(unique(Truth::truth_nu_columns()));
    CHECK(unique(Truth::truth_pf_columns()));
    CHECK(unique(Truth::blobtruth_columns()));
    CHECK(unique(Truth::blobtruth_sum_columns()));
    CHECK(unique(Truth::truth_pc_keys()));
    CHECK(unique(Truth::depo_truth_columns()));

    // SBND attacher order (cfg/pgrapher/experiment/sbnd/docs/sbnd-2step-workflow.md:144-146)
    const std::vector<std::string> sbnd_nu = {
        "nu_idx", "pdg", "ccnc", "mode", "int_type", "flavor",
        "E", "vtx_x", "vtx_y", "vtx_z", "t", "edep"};
    CHECK(Truth::truth_nu_columns() == sbnd_nu);

    const std::vector<std::string> sbnd_pf = {
        "nu_row", "trackid", "parent_trackid", "mother_trackid", "pdg", "process",
        "E", "KE", "start_x", "start_y", "start_z", "start_t",
        "end_x", "end_y", "end_z", "end_t", "start_px", "start_py", "start_pz"};
    const auto& pf = Truth::truth_pf_columns();
    REQUIRE(pf.size() >= sbnd_pf.size());
    for (size_t i = 0; i < sbnd_pf.size(); ++i) {
        CHECK(pf[i] == sbnd_pf[i]);
    }

    CHECK(Truth::tid_ghost == -1);
    CHECK(Truth::tid_shared == -2);
    CHECK(Truth::cls_name(Truth::Cls::michel) == "michel");
    CHECK(Truth::cls_name(Truth::Cls::unknown) == "unknown");
    CHECK(std::string(Truth::truth_pc_name) == "truth");
    CHECK(Truth::truth_pc_keys().front() == "tid");
    CHECK(Truth::depo_truth_columns().size() == 3);
}
