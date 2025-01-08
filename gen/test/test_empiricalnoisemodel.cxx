#include "WireCellAux/Testing.h"

#include "WireCellIface/IChannelSpectrum.h"
#include "WireCellIface/IChannelStatus.h"

#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Persist.h"

#include <boost/crc.hpp>

#include <unordered_set>

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::String;

static
std::ostream& dump_spectra(std::ostream& os, const std::string& fname)
{
    auto jspec = Persist::load(fname);

    std::vector<int> planes, nsamples;
    std::vector<float> periods, gains, shapings, wirelens, consts;
    const int nspecs = jspec.size();
    std::cerr << nspecs << " in " << fname << "\n";
    for (int ind=0; ind<nspecs; ++ind) {
        auto j = jspec[ind];
        planes.push_back(j["plane"].asInt());
        nsamples.push_back(j["nsamples"].asInt());
        periods.push_back(j["period"].asFloat());
        gains.push_back(j["gain"].asFloat());
        shapings.push_back(j["shapings"].asFloat());
        wirelens.push_back(j["wirelens"].asFloat());
        consts.push_back(j["consts"].asFloat());
        Stream::write(os, format("freqs_%03d.npy", ind), convert<std::vector<float>>(j["freqs"]));
        Stream::write(os, format( "amps_%03d.npy", ind), convert<std::vector<float>>(j["amps"]));
        std::cerr << "\t nsamples:" << nsamples.back() << " size=" << j["freqs"].size() << "\n";
    }
    std::cerr << "\t nplanes:" << planes.size() << "\n";
    Stream::write(os, "planes.npy", planes);
    Stream::write(os, "nsamples.npy", nsamples);
    Stream::write(os, "periods.npy", periods);
    Stream::write(os, "gains.npy", gains);
    Stream::write(os, "shapings.npy", shapings);
    Stream::write(os, "wirelens.npy", wirelens);
    Stream::write(os, "consts.npy", consts);
    return os;
}

int main(int argc, char* argv[])
{
    // Testing::loginit(argv[0]);

    auto rng = Testing::get_random();
    auto dft = Testing::get_dft();

    Stream::filtering_ostream out;
    // fixme: need PR #163 merged to save as .npz
    std::string outfile = argc > 1 ? argv[1] : std::string(argv[0]) + ".tar";
    Stream::output_filters(out, outfile);
    
    std::string spectra_file = "protodune-noise-spectra-v1.json.bz2";
    dump_spectra(out, spectra_file);

    size_t nsamples = 2000;     // explicitly NOT an FFT "best" length

    auto anodes = Testing::anodes("pdsp");
    Testing::get_default<IChannelStatus>("StaticChannelStatus");
    {
        auto ienm = Factory::lookup<IConfigurable>("EmpiricalNoiseModel");
        auto cfg = ienm->default_configuration();
        cfg["spectra_file"] = spectra_file;
        cfg["nsamples"] = (unsigned int)nsamples;
        cfg["period"] = 0.5*units::us;
        cfg["wire_length_scale"] = units::cm,
        cfg["anode"] = "AnodePlane:0";
        cfg["dft"] = "FftwDFT";
        cfg["chanstat"] = "StaticChannelStatus";
        ienm->configure(cfg);
    }
    auto enm = Factory::find<IChannelSpectrum>("EmpiricalNoiseModel");

    std::map<uint32_t, IChannelSpectrum::amplitude_t> amps;

    size_t bad_sizes = 0;
    for (auto chid : anodes[0]->channels()) {
        const auto& spec = enm->channel_spectrum(chid);
        if (spec.size() != nsamples) {
            if (! bad_sizes) {
                std::cerr << "EmpiricalNoiseModel is broken: I asked for " << nsamples
                          << " sampled but got " << spec.size()
                          << " on channel ID " << chid
                          << " will repeat this error\n";
            }
            ++bad_sizes;
        }
        //assert(spec.size() == nsamples);

        // We do this craziness because SOMEONE violated the
        // IChannelSpectrum API and return always the same vector.
        boost::crc_32_type crc;
        for (const auto& amp : spec) {
            crc.process_bytes(&amp, sizeof(float));
        }
        auto key = crc.checksum();
        amps[key] = spec;
    }
    if (bad_sizes) {
        std::cerr << "Got " << bad_sizes << " channels that give spectrum sizes that I did not request\n";
    }

    std::cerr << "Got " << amps.size() << " unique spectra\n";
    int count=0;
    for (const auto& [key,amp] : amps) {
        std::cerr << "Writing spec of size: " << amp.size() << "\n";
        Stream::write(out, format("chspec_%03d.npy", count), amp);
        ++count;
    }

    out.flush();
    out.pop();
    std::cerr << "May now wish to run:\n";
    std::cerr << "wirecell-test plot -n empiricalnoisemodel " << outfile << " test_empiricalnoisemodel.pdf\n";
    return 0;
}
