#include "WireCellAux/TensorDMpointtree.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetSource.h"

#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/PointTesting.h"
#include "WireCellUtil/ExecMon.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Logging.h"

#include <cassert>

using spdlog::debug;
using spdlog::info;
using spdlog::error;

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Tree;

static
Points::node_ptr make_simple_pctree(const size_t nnodes = 1000000)
{
    Points::node_ptr root = std::make_unique<Points::node_t>();


    for (size_t count=0; count<nnodes; ++count) {
        root->insert(Points({ {"3d", PointTesting::make_janky_track()} }));
    }
    return root;
}

int test_gen()
{
    ExecMon em("pc tree tensor DM round trip");

    auto root = make_simple_pctree();
    info(em("made pctree"));
    assert(root);

    const std::string datapath = "root";
    auto tens = as_tensors(*root.get(), datapath);
    info(em("convert to tensordm"));
    assert(tens.size() > 0);

    TensorIndex ti(tens);
    info(em("indexed tensors"));

    auto root2 = as_pctree(ti, datapath);
    info(em("convert to pctree"));
    assert(root2);
    
    return 0;
}

int test_file(const std::string& fname,
              const std::string& prefix,
              const std::string& datapath)
{
    ExecMon em("pc tree from tensor DM file");

    auto& pm [[maybe_unused]] = PluginManager::instance();
    assert(pm.add("WireCellSio"));

    auto cobj = Factory::lookup_tn<IConfigurable>("TensorFileSource");
    auto cfg = cobj->default_configuration();
    cfg["inname"] = fname;
    cfg["prefix"] = prefix;
    cobj->configure(cfg);
    auto in = Factory::find_tn<ITensorSetSource>("TensorFileSource");
    info(em("initialize components"));

    ITensorSet::pointer tensp;
    (*in)(tensp);
    info(em("exec TensorFileSource"));

    TensorIndex ti(*tensp->tensors());
    info(em("indexed tensors"));

    auto root = as_pctree(ti, datapath);
    info(em("as pctree"));

    return 0;
}

int main(int argc, char* argv[])
{
    Log::default_logging();

    if (argc == 1) {
        return test_gen();
    }

    if (argc == 4) {
        return test_file(argv[1], argv[2], argv[3]);
    }
    
    return 1;
}
