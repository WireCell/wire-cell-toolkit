// EmpiricalNoiseModel instances configured with the same spectra and
// sampling must share their spectrum caches.  A job makes one model per
// anode; when each kept its own caches, every anode that saw activity added
// tens of MB that were never released (full FD-VD: 320 anodes).

#include "WireCellAux/Testing.h"
#include "WireCellIface/IChannelSpectrum.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

static IChannelSpectrum::pointer make_model(const std::string& name, const std::string& anode_tn,
                                            const std::string& spectra_file)
{
    auto icfg = Factory::lookup<IConfigurable>("EmpiricalNoiseModel", name);
    auto cfg = icfg->default_configuration();
    cfg["spectra_file"] = spectra_file;
    cfg["nsamples"] = 2000;
    cfg["period"] = 0.5 * units::us;
    cfg["wire_length_scale"] = 1.0 * units::cm;
    cfg["anode"] = anode_tn;
    cfg["dft"] = "FftwDFT";
    cfg["chanstat"] = "";
    icfg->configure(cfg);
    return Factory::find<IChannelSpectrum>("EmpiricalNoiseModel", name);
}

TEST_CASE("gen empiricalnoisemodel shares spectra between instances")
{
    auto anodes = Testing::anodes("pdsp");
    REQUIRE(anodes.size() > 1);
    Testing::get_dft();

    const std::string spectra = "protodune-noise-spectra-v1.json.bz2";
    auto a = make_model("ensa", "AnodePlane:0", spectra);
    auto b = make_model("ensb", "AnodePlane:0", spectra);

    const int chid = anodes[0]->channels().front();
    const auto& sa = a->channel_spectrum(chid);
    const auto& sb = b->channel_spectrum(chid);
    REQUIRE(sa.size() > 0);
    CHECK(sa == sb);
    CHECK(&sa == &sb);  // one copy, held by the shared cache

    // A different anode with like channels reuses the same entries.
    auto c = make_model("ensc", "AnodePlane:1", spectra);
    const int chid1 = anodes[1]->channels().front();
    const auto& sc = c->channel_spectrum(chid1);
    const auto& sa1 = a->channel_spectrum(chid);
    CHECK(&sa1 == &sa);  // repeated lookups stay stable
    if (sc == sa) {
        CHECK(&sc == &sa);
    }
}
