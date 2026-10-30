// pr_ml doc 34 sec 8.18: Steiner::Grapher::filter_by_path_constraints read the wrong wire plane for a terminal and took an empty path tree for
// "close in 2-D".  FD-VD evt nc_b0_e1: 92 % of the terminals of CRM 101 deleted, skeleton recall 0.75.
//
// The historical arm is pinned too, because production configs without the knob must keep it.  Reverting the fix makes the own_wpid cases fail.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Units.h"
#include "WireCellClus/SteinerPathFilter.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;
using WireCell::Clus::Steiner::path_filter_removes;
using WireCell::Clus::Steiner::path_filter_terminal_wpid;

namespace {

    using fa_float_t = double;
    using fa_int_t = int;

    // A one-point blob in anode apa at (x, y, z).  Duplicated from the other clus doctests rather than factored into a shared fixture (CLAUDE.md M10).
    void add_blob(Cluster& cl, int apa, double x, double y, double z)
    {
        const int wpid = WirePlaneId(kAllLayers, 0, apa).ident();
        Blob& bl = cl.make_child();
        Dataset ds;
        ds.add("x", Array(std::vector<double>{x}));
        ds.add("y", Array(std::vector<double>{y}));
        ds.add("z", Array(std::vector<double>{z}));
        ds.add("uwire_index", Array(std::vector<int>{0}));
        ds.add("vwire_index", Array(std::vector<int>{0}));
        ds.add("wwire_index", Array(std::vector<int>{0}));
        ds.add("wpid", Array(std::vector<int>{wpid}));
        bl.value().local_pcs()["3d"] = ds;
        bl.value().local_pcs()["scalar"] = Dataset({
            {"charge", Array({(fa_float_t) 1.0})},
            {"center_x", Array({(fa_float_t) x})},
            {"center_y", Array({(fa_float_t) y})},
            {"center_z", Array({(fa_float_t) z})},
            {"wpid", Array({(fa_int_t) wpid})},
            {"npoints", Array({(fa_int_t) 1})},
            {"slice_index_min", Array({(fa_int_t) 0})},
            {"slice_index_max", Array({(fa_int_t) 4})},
            {"u_wire_index_min", Array({(fa_int_t) 0})},
            {"u_wire_index_max", Array({(fa_int_t) 1})},
            {"v_wire_index_min", Array({(fa_int_t) 0})},
            {"v_wire_index_max", Array({(fa_int_t) 1})},
            {"w_wire_index_min", Array({(fa_int_t) 0})},
            {"w_wire_index_max", Array({(fa_int_t) 1})},
            {"max_wire_interval", Array({(fa_int_t) 1})},
            {"min_wire_interval", Array({(fa_int_t) 1})},
            {"max_wire_type", Array({(fa_int_t) 0})},
            {"min_wire_type", Array({(fa_int_t) 0})},
        });
    }

    const double T2 = 1.8 * units::cm, T3 = 6.0 * units::cm;
    const double FAR3D = 10 * units::cm, NEAR3D = 1 * units::cm;

}  // namespace

TEST_CASE("steiner path filter: a terminal's wire plane")
{
    Points::node_t root;
    Grouping& grouping = *root.value.facade<Grouping>();
    Cluster& cl = grouping.make_child();
    // Point index 1 lies in anode 7, far from the x axis.  The cluster point nearest Point(1, 0, 0) is index 0, in anode 3.
    add_blob(cl, 3, 0.0, 0.0, 0.0);
    add_blob(cl, 7, 5000.0, 5000.0, 5000.0);
    REQUIRE(cl.npoints() == 2);

    SUBCASE("the terminal's own wire plane")
    {
        CHECK(path_filter_terminal_wpid(cl, 0, true).apa() == 3);
        CHECK(path_filter_terminal_wpid(cl, 1, true).apa() == 7);
    }
    SUBCASE("the historical lookup converts the index to a point (the defect, kept when the knob is off)")
    {
        CHECK(path_filter_terminal_wpid(cl, 1, false).apa() == 3);
    }
}

TEST_CASE("steiner path filter: the decision")
{
    const double empty[3] = {-1.0, -1.0, -1.0};
    const double close_all[3] = {0.1 * units::cm, 0.1 * units::cm, 0.1 * units::cm};
    const double close_two[3] = {0.1 * units::cm, 0.1 * units::cm, 5 * units::cm};
    const double close_one[3] = {0.1 * units::cm, 5 * units::cm, 5 * units::cm};

    SUBCASE("an empty path tree is not close when the knob is on")
    {
        CHECK_FALSE(path_filter_removes(FAR3D, empty, T2, T3, true));
    }
    SUBCASE("the historical arm reads the -1 as close (the defect, kept when the knob is off)")
    {
        CHECK(path_filter_removes(FAR3D, empty, T2, T3, false));
    }
    SUBCASE("the prototype's rule holds on real distances in both arms")
    {
        for (bool own : {false, true}) {
            CHECK(path_filter_removes(FAR3D, close_all, T2, T3, own));
            CHECK(path_filter_removes(FAR3D, close_two, T2, T3, own));
            CHECK_FALSE(path_filter_removes(FAR3D, close_one, T2, T3, own));     // one plane close: kept
            CHECK_FALSE(path_filter_removes(NEAR3D, close_all, T2, T3, own));    // near in 3-D: kept
            CHECK_FALSE(path_filter_removes(T3, close_all, T2, T3, own));        // exactly at the 3-D threshold: kept
        }
    }
    SUBCASE("near in 3-D is kept even with an empty tree, in both arms")
    {
        for (bool own : {false, true}) {
            CHECK_FALSE(path_filter_removes(NEAR3D, empty, T2, T3, own));
        }
    }
}
