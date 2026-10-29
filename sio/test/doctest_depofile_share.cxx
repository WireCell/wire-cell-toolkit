/** DepoFileSource "share_loaded" (pr_ml doc 33 perf round 2): sources
 *  that read the same file with the same scale hand out ONE loaded copy
 *  of its depo sets instead of one copy each (the clustering / PR jobs
 *  run one source per anode face on the event's drifted depos: 8-18
 *  copies of ~1.8 M depos held at once, 68 % of the PDVD PR job's heap
 *  peak).  Default off: every source loads its own copy, as before.
 *
 *   1. off: two sources give distinct set objects with equal content;
 *   2. on: two sources give the same set object, equal in content to an
 *      unshared read, set by set and then EOS;
 *   3. on: once every source sent EOS and the sets are released, a new
 *      source loads the file again (the cache holds no strong reference);
 *   4. on, another scale: not shared with the first scale.
 */

#include "WireCellSio/DepoFileSink.h"
#include "WireCellSio/DepoFileSource.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellUtil/doctest.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace WireCell;

namespace {

    std::string tmpname(const std::string& base)
    {
        const char* dir = std::getenv("TMPDIR");
        return std::string(dir ? dir : "/tmp") + "/" + base;
    }

    IDepo::vector make_depos(int n)
    {
        IDepo::vector out;
        for (int i = 0; i < n; ++i) {
            auto p = std::make_shared<Aux::SimpleDepo>(2.0 * i, Point(i, 0, 0), -1.0 * (i + 1), nullptr, 0.1, 0.2, 100 + i, 13, 0.5);
            out.push_back(std::make_shared<Aux::SimpleDepo>(1.0 * i, Point(i, 1, 0), -1.0 * (i + 1), p, 0.3, 0.4, 100 + i, 0, 1.0));
        }
        return out;
    }

    void write_sets(const std::string& fname, const std::vector<int>& idents)
    {
        std::remove(fname.c_str());
        Sio::DepoFileSink sink;
        auto cfg = sink.default_configuration();
        cfg["outname"] = fname;
        cfg["truth_sidecar"] = true;
        sink.configure(cfg);
        for (int ident : idents) REQUIRE(sink(std::make_shared<Aux::SimpleDepoSet>(ident, make_depos(3 + ident))));
        sink.finalize();
    }

    struct Src {
        Sio::DepoFileSource src;
        Src(const std::string& fname, bool share, double scale = 1.0)
        {
            auto cfg = src.default_configuration();
            cfg["inname"] = fname;
            cfg["scale"] = scale;
            if (share) cfg["share_loaded"] = true;
            src.configure(cfg);
        }
        IDepoSet::pointer get()
        {
            IDepoSet::pointer ds;
            src(ds);
            return ds;
        }
    };

    void same_content(const IDepoSet::pointer& a, const IDepoSet::pointer& b)
    {
        REQUIRE(a); REQUIRE(b);
        CHECK(a->ident() == b->ident());
        const auto& da = *a->depos();
        const auto& db = *b->depos();
        REQUIRE(da.size() == db.size());
        for (size_t i = 0; i < da.size(); ++i) {
            CHECK(da[i]->time() == db[i]->time());
            CHECK(da[i]->charge() == db[i]->charge());
            CHECK(da[i]->pos() == db[i]->pos());
            CHECK(da[i]->id() == db[i]->id());
            CHECK(da[i]->energy() == db[i]->energy());
            REQUIRE(da[i]->prior()); REQUIRE(db[i]->prior());
            CHECK(da[i]->prior()->pdg() == db[i]->prior()->pdg());
            CHECK(da[i]->prior()->time() == db[i]->prior()->time());
        }
    }
}

TEST_CASE("depofile share_loaded")
{
    const std::string fname = tmpname("wct-depofile-share.tar.bz2");
    write_sets(fname, {1, 2});

    // the reference: one unshared source
    std::vector<IDepoSet::pointer> ref;
    {
        Src r(fname, false);
        for (auto ds = r.get(); ds; ds = r.get()) ref.push_back(ds);
    }
    REQUIRE(ref.size() == 2);

    SUBCASE("off: one copy per source")
    {
        Src a(fname, false), b(fname, false);
        auto a1 = a.get(), b1 = b.get();
        CHECK(a1.get() != b1.get());
        same_content(a1, b1);
        same_content(a1, ref[0]);
    }
    SUBCASE("on: one copy for all sources, set by set, then EOS")
    {
        Src a(fname, true), b(fname, true);
        auto a1 = a.get(), b1 = b.get();
        CHECK(a1.get() == b1.get());
        same_content(a1, ref[0]);
        auto a2 = a.get(), b2 = b.get();
        CHECK(a2.get() == b2.get());
        same_content(a2, ref[1]);
        CHECK(a.get() == nullptr);
        CHECK(b.get() == nullptr);
    }
    SUBCASE("on: released sets are not held; a later source loads again")
    {
        const IDepoSet* first = nullptr;
        {
            Src a(fname, true);
            auto a1 = a.get();
            first = a1.get();
            same_content(a1, ref[0]);
            a.get(); a.get();   // set 2, EOS
        }
        Src c(fname, true);
        auto c1 = c.get();
        same_content(c1, ref[0]);
        CHECK(c1->depos()->size() == ref[0]->depos()->size());
        (void) first;
    }
    SUBCASE("on: a different scale is a different copy")
    {
        Src a(fname, true, 1.0), b(fname, true, -1.0);
        auto a1 = a.get(), b1 = b.get();
        CHECK(a1.get() != b1.get());
        CHECK(b1->depos()->at(0)->charge() == doctest::Approx(-a1->depos()->at(0)->charge()));
    }
    std::remove(fname.c_str());
}
