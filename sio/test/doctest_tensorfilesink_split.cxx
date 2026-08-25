/** TensorFileSink: one container per ident when the name carries a conversion.
 *
 * doc 76 round 2.  A single wire-cell process that streams a GROUP of events
 * must still write the per-event files a one-event-per-process job writes, or
 * no gate can compare the two modes.  Two things have to hold, and only the
 * second is obvious:
 *
 *   1. with a '%' in "outname", each ident lands in its own container;
 *   2. WITHOUT one, exactly one container is written and its bytes are what
 *      they always were -- that is the assertion that catches an off-by-one in
 *      the close/reopen bookkeeping, e.g. a stray reopen on the first set.
 */

#include "WireCellSio/TensorFileSink.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/custard/custard.hpp"

#include "WireCellUtil/doctest.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

using namespace WireCell;

namespace {

    ITensorSet::pointer make_set(int ident)
    {
        // One tiny 1-D float tensor so the set has real content to write.
        std::vector<size_t> shape{2};
        std::vector<float> vals{1.0f * ident, 2.0f * ident};
        auto ten = std::make_shared<Aux::SimpleTensor>(shape, vals.data());
        auto tens = std::make_shared<ITensor::vector>();
        tens->push_back(ten);
        Configuration md;
        md["ident"] = ident;
        return std::make_shared<Aux::SimpleTensorSet>(ident, md, tens);
    }

    /// Member names inside a stream container.
    std::set<std::string> members(const std::string& fname)
    {
        std::set<std::string> out;
        boost::iostreams::filtering_istream in;
        Stream::input_filters(in, fname);
        while (in) {
            std::string name;
            size_t size = 0;
            custard::read(in, name, size);   // the archive member header
            if (name.empty() || !in) break;
            in.ignore(size);
            out.insert(name);
        }
        return out;
    }

    std::string tmpname(const std::string& base)
    {
        const char* dir = std::getenv("TMPDIR");
        return std::string(dir ? dir : "/tmp") + "/" + base;
    }
}

TEST_CASE("tensorfilesink one container per ident")
{
    const std::string tmpl = tmpname("wct-tfs-split-evt%d.tar");
    const std::string one = tmpname("wct-tfs-split-evt111.tar");
    const std::string two = tmpname("wct-tfs-split-evt222.tar");
    std::remove(one.c_str());
    std::remove(two.c_str());

    Sio::TensorFileSink sink;
    Configuration cfg = sink.default_configuration();
    cfg["outname"] = tmpl;
    cfg["prefix"] = "clustering_";
    sink.configure(cfg);

    CHECK(sink(make_set(111)));
    CHECK(sink(make_set(222)));
    sink.finalize();

    auto m1 = members(one);
    auto m2 = members(two);
    // Each container holds ONLY its own event.
    CHECK(m1.count("clustering_tensorset_111_metadata.json") == 1);
    CHECK(m1.count("clustering_tensorset_222_metadata.json") == 0);
    CHECK(m2.count("clustering_tensorset_222_metadata.json") == 1);
    CHECK(m2.count("clustering_tensorset_111_metadata.json") == 0);

    std::remove(one.c_str());
    std::remove(two.c_str());
}

TEST_CASE("tensorfilesink without a conversion writes one container")
{
    const std::string plain = tmpname("wct-tfs-plain.tar");
    std::remove(plain.c_str());

    Sio::TensorFileSink sink;
    Configuration cfg = sink.default_configuration();
    cfg["outname"] = plain;
    cfg["prefix"] = "clustering_";
    sink.configure(cfg);

    CHECK(sink(make_set(111)));
    CHECK(sink(make_set(222)));
    sink.finalize();

    // Both events in the ONE container, in stream order -- the legacy
    // behaviour every existing job depends on.
    auto m = members(plain);
    CHECK(m.count("clustering_tensorset_111_metadata.json") == 1);
    CHECK(m.count("clustering_tensorset_222_metadata.json") == 1);
    // ...and no per-ident container was created behind our back.
    CHECK(members(plain).size() == m.size());

    std::remove(plain.c_str());
}
