
#include "WireCellSpng/TorchContext.h"
#include "WireCellSpng/Testing.h"
#include "WireCellUtil/PluginManager.h"
using namespace WireCell::SPNG;


TEST_CASE("spng torch context") {
    auto& pm = WireCell::PluginManager::instance();
    pm.add("WireCellAux");

    {
        bool agm = torch::GradMode::is_enabled();
        CHECK(agm == true);         // on by default;
    }


    {
        TorchContext tc;
        CHECK(tc.devname() == "cpu");
        CHECK(tc.semname() == "Semaphore:torch-cpu");

        bool agm = torch::GradMode::is_enabled();
        CHECK(agm == true);         // not code-scoped
    }

    TorchContext tc;

    tc.enter();
    {
        bool agm = torch::GradMode::is_enabled();
        CHECK(agm == false);
    }

    tc.exit();
    {
        bool agm = torch::GradMode::is_enabled();
        CHECK(agm == true);
    }

    {
        TorchSemaphore ts(tc);  //  enter()
        bool agm = torch::GradMode::is_enabled();
        CHECK(agm == false);
    } // exit()

    {
        bool agm = torch::GradMode::is_enabled();
        CHECK(agm == true);
    }

}
