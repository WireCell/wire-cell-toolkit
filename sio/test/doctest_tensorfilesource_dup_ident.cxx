/** TensorFileSource: one set per tensorset metadata member, even when idents repeat.
 *
 * A tensor set's ident is the art event number in the lar chain, which restarts at 1 in
 * every subrun; a lar job over several subruns (the gen2 BNB 1M production, one flat CAF
 * = 20 subruns) writes consecutive sets with EQUAL idents into one container.  The
 * reader used to group members by ident alone, so the later set silently replaced the
 * earlier one (sbnd-wirecell-production issue 2, 2026-10-21: 12 of 528 events lost).
 */
#include "WireCellSio/TensorFileSink.h"
#include "WireCellSio/TensorFileSource.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/doctest.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace WireCell;

namespace {
    ITensorSet::pointer make_set(int ident, int tag)
    {
        std::vector<size_t> shape{2};
        std::vector<float> vals{1.0f * tag, 2.0f * tag};
        auto ten = std::make_shared<Aux::SimpleTensor>(shape, vals.data());
        auto tens = std::make_shared<ITensor::vector>();
        tens->push_back(ten);
        Configuration md;
        md["tag"] = tag;
        return std::make_shared<Aux::SimpleTensorSet>(ident, md, tens);
    }
    std::string tmpname(const std::string& base)
    {
        const char* dir = std::getenv("TMPDIR");
        return std::string(dir ? dir : "/tmp") + "/" + base;
    }
    std::vector<int> read_tags(const std::string& fname, std::vector<int>& idents)
    {
        Sio::TensorFileSource src;
        Configuration cfg = src.default_configuration();
        cfg["inname"] = fname;
        cfg["prefix"] = "clustering_";
        src.configure(cfg);
        std::vector<int> tags;
        while (true) {
            ITensorSet::pointer out;
            if (!src(out) || !out) break;   // EOS
            idents.push_back(out->ident());
            tags.push_back(out->metadata()["tag"].asInt());
        }
        return tags;
    }
}

TEST_CASE("tensorfilesource consecutive sets with the same ident stay separate")
{
    const std::string fname = tmpname("wct-tfs-dup-ident.tar");
    std::remove(fname.c_str());
    {
        Sio::TensorFileSink sink;
        Configuration cfg = sink.default_configuration();
        cfg["outname"] = fname;
        cfg["prefix"] = "clustering_";
        sink.configure(cfg);
        // subrun A event 8, subrun B event 8, subrun B event 9: idents 8, 8, 9
        CHECK(sink(make_set(8, 1)));
        CHECK(sink(make_set(8, 2)));
        CHECK(sink(make_set(9, 3)));
        sink.finalize();
    }
    std::vector<int> idents;
    auto tags = read_tags(fname, idents);
    CHECK(tags == std::vector<int>{1, 2, 3});
    CHECK(idents == std::vector<int>{8, 8, 9});
    std::remove(fname.c_str());
}
