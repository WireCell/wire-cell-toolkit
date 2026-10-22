/** DepoFileSink -> DepoFileSource must round-trip the prior chain,
 *  id and pdg, and (with the truth sidecar) the energy.
 *
 *  Before 2026-10-21 the writer stored, on every gen>0 row, the row's
 *  own index as its "child" (Aux::fill), and the reader therefore
 *  called set_prior(row, row): active depos came back with
 *  prior()==nullptr, pdg 0 and energy 1.0, silently.  Three cases:
 *
 *   1. new writer, no sidecar: priors restored, energy is the default;
 *   2. new writer, truth_sidecar on: energy restored too, and a
 *      following set without a sidecar is still read correctly (the
 *      read-ahead member is handed to the next call);
 *   3. a legacy archive (child == own row) still round-trips through
 *      the reader's fallback.
 */

#include "WireCellSio/DepoFileSink.h"
#include "WireCellSio/DepoFileSource.h"
#include "WireCellAux/DepoTools.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/doctest.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace WireCell;

namespace {

    std::string tmpname(const std::string& base)
    {
        const char* dir = std::getenv("TMPDIR");
        return std::string(dir ? dir : "/tmp") + "/" + base;
    }

    // Two active depos: A with a two-deep prior chain carrying id 7,
    // pdg 13 and energy 0.5 on the undrifted end; B with one prior.
    IDepo::vector make_depos()
    {
        auto a2 = std::make_shared<Aux::SimpleDepo>(2.0, Point(0, 0, 0), 1.0, nullptr, 0.1, 0.2, 7, 13, 0.5);
        auto a1 = std::make_shared<Aux::SimpleDepo>(1.0, Point(1, 0, 0), 1.0, a2, 0.1, 0.2, 7, 0, 1.0);
        auto a0 = std::make_shared<Aux::SimpleDepo>(0.0, Point(2, 0, 0), 1.0, a1, 0.1, 0.2, 7, 0, 1.0);
        auto b1 = std::make_shared<Aux::SimpleDepo>(3.0, Point(0, 1, 0), 2.0, nullptr, 0.3, 0.4, 8, 2212, 3.0);
        auto b0 = std::make_shared<Aux::SimpleDepo>(4.0, Point(0, 2, 0), 2.0, b1, 0.3, 0.4, 8, 0, 1.0);
        return {a0, b0};
    }

    void write_sets(const std::string& fname, bool sidecar, const std::vector<int>& idents)
    {
        std::remove(fname.c_str());
        Sio::DepoFileSink sink;
        auto cfg = sink.default_configuration();
        cfg["outname"] = fname;
        cfg["truth_sidecar"] = sidecar;
        sink.configure(cfg);
        for (int ident : idents) {
            auto ds = std::make_shared<Aux::SimpleDepoSet>(ident, make_depos());
            REQUIRE(sink(ds));
        }
        sink.finalize();
    }

    std::vector<IDepoSet::pointer> read_sets(const std::string& fname)
    {
        Sio::DepoFileSource src;
        auto cfg = src.default_configuration();
        cfg["inname"] = fname;
        src.configure(cfg);
        std::vector<IDepoSet::pointer> out;
        while (true) {
            IDepoSet::pointer ds;
            src(ds);
            if (!ds) break;
            out.push_back(ds);
        }
        return out;
    }

    void check_chain(const IDepoSet::pointer& ds, bool expect_energy)
    {
        REQUIRE(ds);
        const auto& depos = *ds->depos();
        REQUIRE(depos.size() == 2);

        auto a0 = depos[0];
        CHECK(a0->id() == 7);
        REQUIRE(a0->prior());
        auto a1 = a0->prior();
        CHECK(a1->time() == doctest::Approx(1.0));
        REQUIRE(a1->prior());
        auto a2 = a1->prior();
        CHECK(a2->time() == doctest::Approx(2.0));
        CHECK(a2->id() == 7);
        CHECK(a2->pdg() == 13);
        CHECK(a2->prior() == nullptr);
        CHECK(a2->energy() == doctest::Approx(expect_energy ? 0.5 : 1.0));

        auto b0 = depos[1];
        CHECK(b0->id() == 8);
        REQUIRE(b0->prior());
        auto b1 = b0->prior();
        CHECK(b1->pdg() == 2212);
        CHECK(b1->prior() == nullptr);
        CHECK(b1->energy() == doctest::Approx(expect_energy ? 3.0 : 1.0));
        CHECK(b1->extent_long() == doctest::Approx(0.3));
    }
}

TEST_CASE("depofile round trip restores priors without sidecar")
{
    const std::string fname = tmpname("wct-depofile-roundtrip-plain.tar");
    write_sets(fname, false, {11, 22});
    auto sets = read_sets(fname);
    REQUIRE(sets.size() == 2);
    CHECK(sets[0]->ident() == 11);
    CHECK(sets[1]->ident() == 22);
    check_chain(sets[0], false);
    check_chain(sets[1], false);
    std::remove(fname.c_str());
}

TEST_CASE("depofile round trip restores energy with sidecar")
{
    const std::string fname = tmpname("wct-depofile-roundtrip-sidecar.tar");
    write_sets(fname, true, {5, 6, 7});
    auto sets = read_sets(fname);
    REQUIRE(sets.size() == 3);
    for (size_t i = 0; i < 3; ++i) {
        CHECK(sets[i]->ident() == 5 + (int)i);
        check_chain(sets[i], true);
    }
    std::remove(fname.c_str());
}

TEST_CASE("depofile reader hands a read-ahead member to the next set")
{
    // set 1 with sidecar, set 2 without: the reader's look-ahead for
    // set 2's sidecar must not swallow set 3's data member.
    const std::string fname = tmpname("wct-depofile-roundtrip-mixed.tar");
    std::remove(fname.c_str());
    {
        boost::iostreams::filtering_ostream out;
        Stream::output_filters(out, fname);
        auto depos = make_depos();
        Array::array_xxf data;
        Array::array_xxi info;
        Eigen::ArrayXXd truth;
        Aux::fill(data, info, depos, true);
        Aux::fill_truth(truth, depos, true);
        Stream::write(out, "depo_data_1.npy", data);
        Stream::write(out, "depo_info_1.npy", info);
        Stream::write(out, "depo_truth_1.npy", truth);
        Stream::write(out, "depo_data_2.npy", data);
        Stream::write(out, "depo_info_2.npy", info);
        Stream::write(out, "depo_data_3.npy", data);
        Stream::write(out, "depo_info_3.npy", info);
        Stream::write(out, "depo_truth_3.npy", truth);
        out.flush();
        out.pop();
    }
    auto sets = read_sets(fname);
    REQUIRE(sets.size() == 3);
    CHECK(sets[0]->ident() == 1);
    CHECK(sets[1]->ident() == 2);
    CHECK(sets[2]->ident() == 3);
    check_chain(sets[0], true);
    check_chain(sets[1], false);
    check_chain(sets[2], true);
    std::remove(fname.c_str());
}

TEST_CASE("depofile reader falls back on legacy self-referencing child rows")
{
    const std::string fname = tmpname("wct-depofile-roundtrip-legacy.tar");
    std::remove(fname.c_str());
    {
        auto depos = make_depos();
        Array::array_xxf data;
        Array::array_xxi info;
        Aux::fill(data, info, depos, true);
        // reproduce the pre-2026-10-21 writer: child == own row on prior rows
        for (int r = 0; r < info.rows(); ++r) {
            if (info(r, 2) > 0) info(r, 3) = r;
        }
        boost::iostreams::filtering_ostream out;
        Stream::output_filters(out, fname);
        Stream::write(out, "depo_data_9.npy", data);
        Stream::write(out, "depo_info_9.npy", info);
        out.flush();
        out.pop();
    }
    auto sets = read_sets(fname);
    REQUIRE(sets.size() == 1);
    CHECK(sets[0]->ident() == 9);
    check_chain(sets[0], false);
    std::remove(fname.c_str());
}
