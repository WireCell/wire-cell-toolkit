
#include "WireCellSpng/TorchContext.h"
#include "WireCellSpng/Testing.h"

using namespace WireCell::SPNG;

TEST_CASE("spng torch context") {
    TorchContext tc;
    CHECK(tc.devname() == "cpu");
    CHECK(tc.semname() == "Semaphore:torch-cpu");

    bool agm = torch::GradMode::is_enabled();
    REQUIRE(agm == true);         // we should not turn this off globally

    tc.enter();
    agm = torch::GradMode::is_enabled();
    REQUIRE(agm == false);

    tc.exit();
    agm = torch::GradMode::is_enabled();
    REQUIRE(agm == true);

    {
        TorchSemaphore ts(tc);  //  enter()
        agm = torch::GradMode::is_enabled();
        REQUIRE(agm == false);
    } // exit()

    agm = torch::GradMode::is_enabled();
    REQUIRE(agm == true);

}
