// Tests of the low-energy depo source and the TA truth node, using the FD HD
// APA 0 geometry if its wires file is on WIRECELL_PATH.

#include "fdhd_anode.h"

#include "WireCellTrig/LowEDepoSource.h"
#include "WireCellTrig/TATruth.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <map>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {
    Configuration ar39_source()
    {
        Configuration src;
        src["name"] = "Ar39";
        src["kind"] = "decay";
        src["activity_bq_per_kg"] = 1.01;
        src["branches"][0]["beta"]["q"] = 0.565 * units::MeV;
        src["branches"][0]["beta"]["z"] = 19;
        src["branches"][0]["beta"]["shape"] = "unique_first";
        return src;
    }

    std::vector<double> column(const ITensor::pointer& ten, const std::string& name)
    {
        const auto md = ten->metadata();
        const size_t ncols = md["columns"].size();
        size_t icol = ncols;
        for (size_t ind = 0; ind < ncols; ++ind)
            if (md["columns"][(int) ind].asString() == name) icol = ind;
        REQUIRE(icol < ncols);
        const double* data = (const double*) ten->data();
        const size_t nrows = ten->shape()[0];
        std::vector<double> col(nrows);
        for (size_t row = 0; row < nrows; ++row) col[row] = data[row * ncols + icol];
        return col;
    }
}  // namespace

TEST_CASE("trig lowe depo source")
{
    auto anode = trig_test::fdhd_anode();
    if (!anode) {
        MESSAGE("FD HD wires file not on WIRECELL_PATH, skipping");
        return;
    }
    LowEDepoSource src;
    auto cfg = src.default_configuration();
    cfg["anode"] = "AnodePlane:trigtest";
    cfg["recombination"] = trig_test::box_recombination();
    cfg["nframes"] = 2;
    cfg["readout_time"] = 0.2 * units::ms;
    cfg["pre_time"] = 0.1 * units::ms;
    cfg["sources"][0] = ar39_source();
    src.configure(cfg);
    // two faces of about 3.6 x 6 x 2.3 m3 of argon: about 140 t
    CHECK(src.mass() / units::kg > 1.2e5);
    CHECK(src.mass() / units::kg < 1.6e5);

    const double expected = 1.01 * (src.mass() / units::kg) * 0.3e-3;  // decays per frame
    for (int frame = 0; frame < 2; ++frame) {
        IDepoSet::pointer ds;
        REQUIRE(src(ds));
        REQUIRE(ds);
        CHECK(ds->ident() == frame);
        std::map<int, int> markers, ndepos;
        double last = -1e100;
        for (const auto& depo : *ds->depos()) {
            CHECK(depo->time() >= last);
            last = depo->time();
            CHECK(depo->id() % LowEDepoSource::max_sources == 0);
            CHECK(depo->time() >= -0.1 * units::ms);
            CHECK(depo->time() < 0.2 * units::ms);
            if (depo->charge() == 0) {
                ++markers[depo->id()];
                CHECK(depo->energy() > 0);
                CHECK(depo->energy() <= 0.565 * units::MeV);
            }
            else {
                ++ndepos[depo->id()];
                CHECK(depo->charge() < 0);  // electrons
            }
        }
        const double n = markers.size();
        CHECK(std::abs(n - expected) < 5 * std::sqrt(expected));
        for (const auto& [id, count] : markers) CHECK(count == 1);
        for (const auto& [id, count] : ndepos) CHECK(markers.count(id) == 1);
    }
    IDepoSet::pointer ds;
    REQUIRE(src(ds));
    CHECK(ds == nullptr);
    CHECK(!src(ds));
}

TEST_CASE("trig ta truth")
{
    auto anode = trig_test::fdhd_anode();
    if (!anode) {
        MESSAGE("FD HD wires file not on WIRECELL_PATH, skipping");
        return;
    }
    const double speed = 1.6 * units::mm / units::us;
    const double offset_w = 80 * units::us;

    TATruth truth;
    auto cfg = truth.default_configuration();
    cfg["anode"] = "AnodePlane:trigtest";
    cfg["drift_speed"] = speed;
    cfg["plane_offsets"][2] = offset_w;
    cfg["sources"][0]["name"] = "snb";
    cfg["sources"][0]["kind"] = "marley";
    cfg["sources"][1] = ar39_source();
    cfg["metadata"]["boost"] = 100;
    truth.configure(cfg);

    // A signal interaction (source 0) and a background one (source 1) at
    // the same point in face 0, at different times.
    auto face = anode->faces()[0];
    auto wplane = face->planes()[2];
    const double xresp = wplane->pimpos()->origin().x();
    const Point pos(1 * units::m, -3 * units::m, 1 * units::m);
    REQUIRE(face->sensitive().inside(pos));
    const int wind = wplane->pimpos()->region_binning().bin(wplane->pimpos()->distance(pos, 2));
    const int channel = wplane->wires()[wind]->channel();
    const int sig_id = 64 * 7 + 0, bkg_id = 64 * 8 + 1;
    const double tsig = 0.1 * units::ms, tbkg = 1.0 * units::ms;
    IDepo::vector depos{
        std::make_shared<Aux::SimpleDepo>(tsig, pos, 0.0, nullptr, 0, 0, sig_id, 0, 15 * units::MeV),
        std::make_shared<Aux::SimpleDepo>(tsig, pos, -5000.0, nullptr, 0, 0, sig_id, 11, 0.2 * units::MeV),
        std::make_shared<Aux::SimpleDepo>(tbkg, pos, 0.0, nullptr, 0, 0, bkg_id, 0, 0.3 * units::MeV),
        std::make_shared<Aux::SimpleDepo>(tbkg, pos, -2000.0, nullptr, 0, 0, bkg_id, 11, 0.1 * units::MeV),
    };
    auto ds = std::make_shared<Aux::SimpleDepoSet>(3, depos);

    const double arrive = (pos.x() - xresp) / speed + offset_w;
    auto make_ta = [&](double t, int ch) {
        TriggerActivity ta;
        TriggerPrimitive tp;
        tp.channel = ch;
        tp.time_start = t;
        tp.samples_over_threshold = 10;
        ta.inputs.push_back(tp);
        ta.ninputs = 1;
        ta.time_start = t;
        ta.time_end = t + 5 * units::us;
        ta.channel_start = ta.channel_end = ch;
        ta.plane = wplane->planeid();
        return ta;
    };
    ITriggerActivitySet::activities_t tas{make_ta(tsig + arrive - 1 * units::us, channel),
                                          make_ta(tbkg + arrive + 1 * units::us, channel + 1),
                                          make_ta(tbkg + arrive, channel + 100)};
    auto taset = std::make_shared<Aux::SimpleTriggerActivitySet>(3, 0, 0, 3 * units::ms, tas);

    ITriggerActivityTruth::output_pointer out;
    REQUIRE(truth(std::make_tuple(ds, taset), out));
    REQUIRE(out);
    CHECK(out->ident() == 3);
    CHECK(out->metadata()["nchannels"].asInt() == 2560);
    CHECK(out->metadata()["mass_kg"].asDouble() > 1.2e5);
    CHECK(out->metadata()["plane_channels"][2].size() == 960);
    CHECK(out->metadata()["plane_channels"][0].size() == 800);
    CHECK(out->metadata()["mass_kg"].asDouble() < 1.6e5);
    CHECK(out->metadata()["metadata"]["boost"].asInt() == 100);
    CHECK(out->metadata()["sources"][0]["signal"].asBool());
    CHECK(!out->metadata()["sources"][1]["signal"].asBool());
    auto tensors = out->tensors();
    REQUIRE(tensors->size() == 3);
    auto tat = tensors->at(0), intt = tensors->at(1), lnk = tensors->at(2);
    CHECK(tat->metadata()["name"].asString() == "tas");

    auto qtot = column(tat, "q_total"), qsig = column(tat, "q_signal"), dom = column(tat, "dominant_id");
    auto off = column(tat, "offset_est");
    REQUIRE(qtot.size() == 3);
    CHECK(qtot[0] == doctest::Approx(5000));
    CHECK(qsig[0] == doctest::Approx(5000));
    CHECK(dom[0] == sig_id);
    CHECK(off[0] == doctest::Approx(offset_w - 1 * units::us));
    CHECK(qtot[1] == doctest::Approx(2000));  // within the channel tolerance
    CHECK(qsig[1] == 0);
    CHECK(dom[1] == bkg_id);
    CHECK(qtot[2] == 0);
    CHECK(dom[2] == -1);

    auto ids = column(intt, "id"), signal = column(intt, "signal"), energy = column(intt, "energy");
    auto qcol = column(intt, "q_collection");
    REQUIRE(ids.size() == 2);
    CHECK(ids[0] == sig_id);
    CHECK(signal[0] == 1);
    CHECK(energy[0] == doctest::Approx(15 * units::MeV));
    CHECK(qcol[0] == doctest::Approx(5000));
    CHECK(column(intt, "t_collection")[0] == doctest::Approx(tsig + arrive));
    CHECK(signal[1] == 0);
    CHECK(column(lnk, "q").size() == 2);

    // end of stream
    ITriggerActivityTruth::output_pointer eos;
    REQUIRE(truth(std::make_tuple(IDepoSet::pointer(), ITriggerActivitySet::pointer()), eos));
    CHECK(eos == nullptr);
}

TEST_CASE("trig lowe depo source daughters")
{
    auto anode = trig_test::fdhd_anode();
    if (!anode) {
        MESSAGE("FD HD wires file not on WIRECELL_PATH, skipping");
        return;
    }
    LowEDepoSource src;
    auto cfg = src.default_configuration();
    cfg["anode"] = "AnodePlane:trigtest";
    cfg["recombination"] = trig_test::box_recombination();
    cfg["readout_time"] = 10 * units::ms;
    cfg["pre_time"] = 0;
    const double mean_life = 237 * units::us;
    Configuration bi;
    bi["name"] = "Bi214";
    bi["kind"] = "decay";
    bi["rate_hz"] = 2000;
    bi["branches"][0]["beta"]["q"] = 3.27 * units::MeV;
    bi["branches"][0]["beta"]["z"] = 84;
    bi["daughter"]["name"] = "Po214";
    bi["daughter"]["mean_life"] = mean_life;
    Configuration po;
    po["name"] = "Po214";
    po["kind"] = "decay";
    po["rate_hz"] = 0;
    po["branches"][0]["alpha"] = 7.6869 * units::MeV;
    cfg["sources"][0] = bi;
    cfg["sources"][1] = po;
    src.configure(cfg);

    IDepoSet::pointer ds;
    REQUIRE(src(ds));
    std::map<int, IDepo::pointer> markers;
    for (const auto& depo : *ds->depos())
        if (depo->charge() == 0) markers[depo->id()] = depo;
    std::vector<IDepo::pointer> parents, daughters;
    for (const auto& [id, m] : markers) (id % 64 == 0 ? parents : daughters).push_back(m);
    REQUIRE(parents.size() > 10);
    // all but the parents near the window end have a daughter
    CHECK(daughters.size() <= parents.size());
    CHECK(daughters.size() >= parents.size() * 0.85);
    // each daughter sits on a parent, later
    double dtsum = 0;
    for (const auto& d : daughters) {
        bool found = false;
        for (const auto& p : parents) {
            if ((d->pos() - p->pos()).magnitude() < 1e-6 && d->time() >= p->time()) {
                found = true;
                dtsum += d->time() - p->time();
                break;
            }
        }
        CHECK(found);
    }
    CHECK(dtsum / daughters.size() == doctest::Approx(mean_life).epsilon(0.3));

    cfg["sources"][0]["daughter"]["name"] = "nope";
    CHECK_THROWS(src.configure(cfg));
}
