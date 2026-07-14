#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellAux/LinterpFunction.h"

using namespace WireCell;
using namespace WireCell::Aux;

using spdlog::debug;

TEST_CASE("aux linterp function") {

    Configuration reg_cfg;
    reg_cfg["values"][0] = 0.0; // 0 = x
    reg_cfg["values"][1] = 2.0; // 1
    reg_cfg["values"][2] = 3.0; // 2
    reg_cfg["values"][3] = -1.0; // 3
    Configuration irr_cfg = reg_cfg;
    irr_cfg["coords"][0] = 0.0; // = x
    irr_cfg["coords"][1] = 10.0;
    irr_cfg["coords"][2] = 11.0;
    irr_cfg["coords"][3] = 20.;

    reg_cfg["start"] = 0.0;
    reg_cfg["step"] = 1.0;

    LinterpFunction reg, irr;

    {
        auto tmp = reg.default_configuration();
        reg.configure(update(tmp, reg_cfg));
    }
    {
        auto tmp = irr.default_configuration();
        irr.configure(update(tmp, irr_cfg));
    }

    // first point exact
    REQUIRE(0.0 == reg.scalar_function(0.0));
    REQUIRE(0.0 == irr.scalar_function(0.0));

    // second point exact
    REQUIRE(2.0 == reg.scalar_function(1.0));
    REQUIRE(2.0 == irr.scalar_function(10.0));

    // extrapolation
    REQUIRE(0.0 == reg.scalar_function(-1.0));
    REQUIRE(0.0 == irr.scalar_function(-1.0));
    REQUIRE(-1.0 == reg.scalar_function(+4.0));
    REQUIRE(-1.0 == irr.scalar_function(30.0));

    // midway interpolation
    REQUIRE(1.0 == reg.scalar_function(0.5));
    REQUIRE(1.0 == irr.scalar_function(5.0));

}
