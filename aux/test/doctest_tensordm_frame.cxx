#include "WireCellAux/TensorDMframe.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/doctest.h"

#include <vector>
#include <iostream>

using namespace WireCell;
using namespace WireCell::String;
namespace PC = WireCell::PointCloud;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;
using spdlog::debug;

template<typename T>
void dump_array(const PC::Dataset& ds, const std::string& name, const std::string& prefix = "\t\t")
{
    REQUIRE(ds.has(name));
    // if (! ds.has(name)) {
    //     debug("tensordm frame: {} no array \"{}\"", prefix, name);
    //     return;
    // }

    auto arr = ds.get(name);
    REQUIRE(arr->is_type<T>());

    std::stringstream ss;
    ss << "doctest_tensordm_frame: array: " << prefix << "\n";
    ss << "\tarray \"" << name << "\" size major: " << arr->size_major() << " dtype: " << arr->dtype() << "\n";
    ss << "\tmetadata: " << arr->metadata() << "\n";
    auto vec = arr->elements<T>();
    ss << vec.size() << " " << name;
    for (auto one : vec) ss << " " << one;
    debug(ss.str());

}

#if 0
static void dump_tensor(ITensor::pointer iten, const std::string& prefix = "")
{
    std::stringstream ss;
    ss << "doctest_tensordm_frame: tensor: " << prefix << "\n";

    const std::vector<std::string> keys = {"datatype","datapath","tag","arrays"};

    auto md = iten->metadata();
    for (const auto& key : keys) {
        auto val = md[key];
        if (val.isNull()) continue;
        ss << "\t" << key << " = " << val << "\n";
    }
    ss << "\t"
       << "array: dtype=" << iten->dtype()
       << " shape=[";
    for (auto s : iten->shape()) ss << " " << s;
    ss << " ]\n";
    debug(ss.str());
}
#endif

void test_tensor(IFrame::pointer frame, FrameTensorMode mode, bool truncate,
                 const size_t ntraces, const size_t nchans, const std::string& frame_tag)
{
    std::stringstream ss;
    ss << "test frame/tensor round trip with"
       << " mode=" << (int)mode
       << " truncate=" << truncate
       << " ntraces=" << ntraces
       << " nchans " << nchans
       << " frame_tag=" << frame_tag;
    debug(ss.str());

    const std::string datapath = "frames/0";

    auto itenvec = as_tensors(frame, datapath, mode, truncate);
    CHECK(itenvec.size() > 0);
    {
        auto itenset = as_tensorset(itenvec);
        CHECK(itenset);
    }

    auto ften = itenvec[0];
    auto fmd = ften->metadata();
    CHECK(fmd["tags"][0].asString() == frame_tag);
    CHECK(fmd["datatype"].asString() == "frame");
    CHECK(fmd["datapath"].asString() == "frames/0");

    TensorIndex ti(itenvec);
    auto ften2 = ti.at_of("frame");
    CHECK(ften == ften2);

    CHECK(ti.at(datapath) == ften);

    CHECK(fmd["traces"].size() == ntraces);
    size_t count=0;
    for (auto jtrdp : fmd["traces"]) {
        debug("traces [{} / {}]", count++, ntraces);

        auto iten = ti.at(jtrdp.asString());
        CHECK(iten);

        // dump_tensor(iten,"\t");
        auto md = iten->metadata();
        CHECK(startswith(md["datapath"].asString(), datapath));
        CHECK(md["datatype"].asString() == "trace");
    }

    count = 0;
    for (auto jtrdp : fmd["tracedata"]) {
        debug("tracedata [{} / {}]", count++, fmd["tracedata"].size());
        auto iten = ti.at(jtrdp.asString());
        CHECK(iten);

        // dump_tensor(iten, "\t");
        auto md = iten->metadata();
        auto tag = get<std::string>(md, "tag", "");

        auto dp = md["datapath"].asString();
        CHECK(startswith(dp, datapath));
        CHECK(md["datatype"].asString() == "pcdataset");

        auto ds = as_dataset(ti, dp);
        auto keys = ds.keys();
        if (tag.empty()) {
            // convention to save whole frame, must have chids
            REQUIRE(ds.has("chid"));
            auto chid = ds.get("chid");
            REQUIRE(chid);
            REQUIRE(chid->size_major() > 0);
            REQUIRE(chid->is_type<int>());
            dump_array<int>(ds, "chid","\t");
        }
        else {
            REQUIRE(ds.has("index"));
            auto arr = ds.get("index");
            REQUIRE(arr);
            REQUIRE(arr->size_major() > 0);
            REQUIRE(arr->is_type<size_t>());
            dump_array<size_t>(ds, "index","\t");
            if (ds.has("chid")) {
                // may have chids
                dump_array<int>(ds, "chid","\t");
            }
        }
    }

    // round trip back to frame
    auto frame2 = as_frame(ti, datapath);
    CHECK(frame2->ident() == frame->ident());
    CHECK(frame2->time() == frame->time());
    CHECK(frame2->tick() == frame->tick());
    {
        // check traces
        auto ts1 = frame->traces();
        auto ts2 = frame2->traces();
        CHECK(ts1->size() == ts2->size());
        std::set<int> chids1, chids2;
        for (size_t ind=0; ind<ts1->size(); ++ind) {
            auto t1 = ts1->at(ind);
            auto t2 = ts2->at(ind);
            chids1.insert(t1->channel());
            chids2.insert(t2->channel());

            std::stringstream ss;
            ss << "trace: " << ind << "/" << ts1->size()
               << " chan: " << t1->channel() << " =?= " << t2->channel()
               << ", tbin: " << t1->tbin() << " =?= " << t2->tbin();
            debug(ss.str());

            if (mode == FrameTensorMode::sparse) {
                // Only get exact representation round-trip in sparse mode
                CHECK(t1->channel() == t2->channel());
                CHECK(t1->tbin() == t2->tbin());

                const auto& qs1 = t1->charge();
                const auto& qs2 = t2->charge();
                // std::cerr << "trace: " << ind << "/" << ts1->size()
                //           << " sizes: " << qs1.size() << " =?= " << qs2.size() << "\n";
                CHECK(qs1.size() == qs2.size());

                for (size_t iq=0; iq < qs1.size(); ++iq) {
                    float q1 = qs1[iq];
                    float q2 = qs2[iq];
                    if (truncate) {
                        // I/O truncated to int16_t, so must do likewise in this check
                        q1 = (int16_t)q1;
                    }

                    // if (q1 != q2) {
                    //     std::cerr << "charge mismatch: " <<  q1 << " != " <<  q2 << "\n";
                    // }
                    CHECK(q1 == q2);
                }
            }
        }
        CHECK(chids1 == chids2);
    }
    {
        // check channel mask map
        auto m1 = frame->masks();
        auto m2 = frame2->masks();
        CHECK(m1.size() == m2.size());
        for (const auto& [tag, cms] : m1) {
            auto it2 = m2.find(tag);
            CHECK(it2 != m2.end());
            CHECK(it2->second.size() == cms.size());
        }
        for (const auto& [tag, cms] : m2) {
            auto it1 = m1.find(tag);
            CHECK(it1 != m1.end());
            CHECK(it1->second.size() == cms.size());
        }
    }
    // only get exact representation round trip in sparse mode
    if (mode == FrameTensorMode::sparse) {
        // check tagged traces and summaries
        const auto tts1 = frame->trace_tags();
        const auto tts2 = frame2->trace_tags();
        CHECK(tts1.size() == tts2.size());
        //std::set ttss1(tts1.begin(), tts1.end());
        //std::set ttss2(tts2.begin(), tts2.end());
	//	std::set ttss1(tts1);
	//std::set ttss2(tts2);
        //CHECK(ttss1 == ttss2);

        for (const auto& tag : frame->trace_tags()) {
            const auto& tt1 = frame->tagged_traces(tag);
            const auto& tt2 = frame2->tagged_traces(tag);
            // std::cerr << "tagged traces: " << tt1.size() << " =?= " << tt2.size() << "\n";
            CHECK(tt1.size() == tt2.size());
            for (size_t ind=0; ind<tt1.size(); ++ind) {
                CHECK(tt1[ind] == tt2[ind]);
            }

            const auto& ts1 = frame->trace_summary(tag);
            const auto& ts2 = frame2->trace_summary(tag);
            CHECK(ts1.size() == ts2.size());
            for (size_t ind=0; ind<ts1.size(); ++ind) {
                CHECK(ts1[ind] == ts2[ind]);
            }
        }
    }

    // And one more back to tensor
    auto itenvec2 = as_tensors(frame2, datapath, mode, truncate);
    CHECK(itenvec2.size() > 0);
    if (mode == FrameTensorMode::sparse) {
        auto& tens1 = itenvec;
        auto& tens2 = itenvec2;
        CHECK(tens1.size() == tens2.size());
        for (size_t ind=0; ind<tens1.size(); ++ind) {
            auto ten1 = tens1[ind];
            auto ten2 = tens2[ind];
            CHECK(ten1->metadata() == ten2->metadata());
            // std::cerr << "tensor sizes: " << ten1->size() << " =?= " << ten2->size() << "\n";
            CHECK(ten1->size() == ten2->size());            
        }
    }

}


void test_tensors(IFrame::pointer signal, IFrame::pointer noise)
{
    const size_t nchan = 4;
    const size_t nunified = 1;
    const size_t nsparse = nchan;
    const size_t ntagged = 2;

    // 3x2x2

    const std::vector<bool> tf = {true, false};

    for (auto truncate : tf) {
        // signal has 4 channels, 2 unique
        test_tensor(signal, FrameTensorMode::unified, truncate, nunified, nchan, "signal");
        test_tensor(signal, FrameTensorMode::tagged,  truncate, ntagged,  nchan, "signal");
        test_tensor(signal, FrameTensorMode::sparse,  truncate, nsparse,  nchan, "signal");

        // noise has 4 unique channels
        test_tensor(noise, FrameTensorMode::unified, truncate, nunified, nchan, "noise");
        test_tensor(noise, FrameTensorMode::tagged,  truncate, ntagged,  nchan, "noise");
        test_tensor(noise, FrameTensorMode::sparse,  truncate, nsparse,  nchan, "noise");
    }
}

TEST_CASE("tensordm frame")
{
    debug("doctest_tensordm_frame");
    const double tick = 0.5 * units::us;
    const double readout = 5 * units::ms;
    const int nsamples = round(readout / tick);

    // "Signal"
    const int signal_frame_ident = 100;
    const double signal_start_time = 6 * units::ms;
    std::vector<float> signal{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 0.0};

    ITrace::vector signal_traces{std::make_shared<SimpleTrace>(1, 10, signal),
                                 std::make_shared<SimpleTrace>(1, 20, signal),
                                 std::make_shared<SimpleTrace>(2, 30, signal),
                                 std::make_shared<SimpleTrace>(2, 33, signal)};
    {
        auto sbr = tbin_range(signal_traces);
        debug("doctest_tensordm_frame: signal tbin range=[{},{}]", sbr.first, sbr.second);
        CHECK(sbr.first == 10);
        CHECK(sbr.second == 44);
    }

    auto f1 = std::make_shared<SimpleFrame>(signal_frame_ident, signal_start_time, signal_traces, tick);
    f1->tag_frame("signal");
    f1->tag_traces("oddch", {0,2});
    f1->tag_traces("evench", {1,3});

    // "Noise"
    const int noise_frame_ident = 0;
    const double noise_start_time = 5 * units::ms;
    std::vector<float> noise(nsamples, -0.5);
    ITrace::vector noise_traces{std::make_shared<SimpleTrace>(0, 0, noise),
                                std::make_shared<SimpleTrace>(1, 0, noise),
                                std::make_shared<SimpleTrace>(2, 0, noise),
                                std::make_shared<SimpleTrace>(3, 0, noise)};
    {
        auto sbr = tbin_range(noise_traces);
        debug("doctest_tensordm_frame: noise tbin range=[{},{}]", sbr.first, sbr.second);
        CHECK(sbr.first == 0);
        CHECK(sbr.second == nsamples);
    }
    auto f2 = std::make_shared<SimpleFrame>(noise_frame_ident, noise_start_time, noise_traces, tick);
    f2->tag_frame("noise");
    f2->tag_traces("oddch", {0,2});
    f2->tag_traces("evench", {1,3});

    const int summed_ident = 1;
    auto f3 = Aux::sum({f1, f2}, summed_ident);
    CHECK(f3);
    CHECK(f3->ident() == summed_ident);
    CHECK(f3->time() == noise_start_time);
    CHECK(f3->tick() == tick);
    auto traces = f3->traces();
    CHECK(traces->size() == 4);
    for (auto trace : *traces) {
        const auto& charge = trace->charge();
        for (size_t ind = 0; ind < charge.size(); ++ind) {
            const float q = charge[ind];
            if (q < 0) {
                continue;
            }
        }
    }

    test_tensors(f1, f2);
}
