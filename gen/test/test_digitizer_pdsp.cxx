// Simple tests of the digitizer.
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Testing.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IFrameFilter.h" // digitizer

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

using namespace WireCell;
using namespace WireCell::Aux;

using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

int main(int argc, char* argv[])
{
    std::string outname = argv[0];
    outname += ".json";

    auto params = Persist::load("pgrapher/experiment/pdsp/params.jsonnet");

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");
    {
        auto icfg = Factory::lookup_tn<IConfigurable>("WireSchemaFile");
        auto cfg = icfg->default_configuration();
        cfg["filename"] = "protodune-wires-larsoft-v4.json.bz2";
        icfg->configure(cfg);
    }
    {
        auto icfg = Factory::lookup_tn<IConfigurable>("AnodePlane");
        auto cfg = icfg->default_configuration();
        cfg["ident"] = 0;
        cfg["wire_schema"] ="WireSchemaFile";
        cfg["faces"] = params["det"]["volumes"][0]["faces"];
        icfg->configure(cfg);
    }

    Configuration jcfg;
    {
        auto icfg = Factory::lookup_tn<IConfigurable>("Digitizer");
        jcfg = icfg->default_configuration();
        jcfg["gain"] = 1.0;
        jcfg["baselines"][0] = 0;
        jcfg["resolution"] = 12;
        jcfg["fullscale"][0] = 0;
        jcfg["fullscale"][1] = 2.0*units::volt;
        icfg->configure(jcfg);
    }
    auto digi = Factory::lookup_tn<IFrameFilter>("Digitizer");


    const double tick = 0.5*units::us;

    auto jfs = jcfg["fullscale"];
    const double vmin = jfs[0].asDouble();
    const double vmax = jfs[1].asDouble();
    const int res = jcfg["resolution"].asInt();
    const size_t adcmax = 1<<res;

    const double oversample = 10;
    const size_t nsamples = oversample * adcmax;
    const double dv = (vmax-vmin)/(oversample*adcmax);
    std::vector<float> ramp(nsamples);
    for (size_t ind=0; ind<nsamples; ++ind) {
        ramp[ind] = vmin + ind*dv;
    }
    ITrace::vector ramp_traces{std::make_shared<SimpleTrace>(0, 0, ramp)};

    auto iframe = std::make_shared<SimpleFrame>(0, 0, ramp_traces, tick);
    IFrame::pointer oframe;
    bool ok = (*digi)(iframe,oframe);
    AssertMsg(ok, "digitizer failed");

    auto oramp = oframe->traces()->at(0)->charge();

    // for (size_t ind=0; ind<adcmax; ++ind) {
    //     auto q = oramp[ind];
    //     std::cerr << ind << " " << ramp[ind] << " " << q << " " << round(q) << " " << floor(q) << "\n";
    // }

    Configuration jout;
    put(jout, "volts", ramp);
    put(jout, "adc", oramp);
    Persist::dump(outname, jout, true);
    std::cerr << outname << "\n";

    return 0;
}
