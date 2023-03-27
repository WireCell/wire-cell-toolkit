#include "WireCellAux/TaggedFrameTensorSet.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/Util.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellUtil/Logging.h"

#include <unordered_set>

WIRECELL_FACTORY(TaggedFrameTensorSet, WireCell::Aux::TaggedFrameTensorSet, WireCell::IFrameTensorSet,
                 WireCell::IConfigurable)

using namespace WireCell;

WireCell::Configuration Aux::TaggedFrameTensorSet::default_configuration() const
{
    // Each tag will produce a tensor.
    Configuration cfg;
    cfg["tensors"] = Json::arrayValue;
    return cfg;
}

void Aux::TaggedFrameTensorSet::configure(const WireCell::Configuration& config) { m_cfg = config; }

bool Aux::TaggedFrameTensorSet::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {  // pass on EOS
        return true;
    }

    auto log = Log::logger("tens");

    ITensor::vector* itv = new ITensor::vector;

    Configuration set_md;
    set_md["time"] = in->time();
    set_md["tick"] = in->tick();

    std::unordered_set<std::string> tags_seen;

    for (auto jten : m_cfg["tensors"]) {
        auto tag = get<std::string>(jten, "tag", "");
        if (tags_seen.find(tag) == tags_seen.end()) {
            tags_seen.insert(tag);
        }
        else {
            log->warn("Frame->Tensor: skipping duplicate tag '{}'", tag);
            continue;
        }

        auto traces = Aux::tagged_traces(in, tag);
        if (traces.empty()) {
            log->warn("Frame->Tensor: no traces for tag '{}', skipping", tag);
            continue;
        }
        size_t ntraces = traces.size();
        // log->debug("Frame->Tensor: {} traces for tag '{}'", ntraces, tag);

        float pad = get<float>(jten, "pad", 0.0);

        auto mm_tbin = Aux::tbin_range(traces);
        int tbin0 = mm_tbin.first;

        // traces may be degenerate/overlapping in channels and ticks.
        // Do a little dance to map chid to an index sorted by chid.
        std::vector<int> uchid;
        for (const auto& trace : traces) {
            uchid.push_back(trace->channel());
        }
        // log->debug("Frame->Tensor: tag '{}': {} channels total, {} traces", tag, uchid.size(), ntraces);

        std::sort(uchid.begin(), uchid.end());
        std::vector<int>::iterator it = std::unique(uchid.begin(), uchid.end());
        const size_t nchans = std::distance(uchid.begin(), it);
        // log->debug("Frame->Tensor: tag '{}': {} channels unique", tag, nchans);
        uchid.resize(nchans);

        // map channel ID to consecutive, sorted index
        std::map<int, size_t> chid2ind;
        for (auto chid : uchid) {
            size_t ind = chid2ind.size();
            chid2ind[chid] = ind;
        }

        // "channels" tensor
        SimpleTensor* cht = new SimpleTensor({nchans}, (int*)nullptr);
        Eigen::Map<Eigen::ArrayXi> charr((int*) cht->data(), nchans);
        assert((size_t) charr.size() == nchans);
        charr.setZero();

        // "summary" tensor, create but only fill and save if we have summary
        SimpleTensor* sumt = new SimpleTensor({nchans}, (double*)nullptr);
        Eigen::Map<Eigen::ArrayXd> sumarr((double*) sumt->data(), nchans);
        assert((size_t) sumarr.size() == nchans);
        sumarr.setZero();

        const auto& ts = in->trace_summary(tag);
        bool have_summary = ts.size() > 0;

        size_t nticks = mm_tbin.second - mm_tbin.first;
        const std::vector<size_t> shape = {nchans, nticks};

        // "waveform" tensor
        SimpleTensor* st = new SimpleTensor(shape, (float*)nullptr);
        Eigen::Map<Eigen::ArrayXXf> arr((float*) st->data(), nchans, nticks);
        arr.setConstant(pad);

        for (size_t itr = 0; itr < ntraces; ++itr) {
            auto trace = traces[itr];
            const auto& qarr = trace->charge();
            size_t nq = qarr.size();

            int chid = trace->channel();
            auto it = chid2ind.find(chid);
            assert(it != chid2ind.end());
            size_t chind = it->second;
            assert(chind < (size_t) charr.size());

            // log->debug("itr: {}/{}, chid: {}, chind: {}", itr, ntraces, chid, chind);
            charr[chind] = chid;
            if (have_summary) {
                sumarr[chind] += ts[itr];
            }

            size_t tbini = trace->tbin() - tbin0;  // where in output to lay-down charge.
            size_t tbinf = std::min(nticks, tbini + nq);
            size_t tbinn = tbinf - tbini;

            Eigen::Map<const Eigen::ArrayXf> wave(qarr.data(), nq);
            arr.row(chind).segment(tbini, tbinn) += wave;
        }

        auto& wf_md = st->metadata();
        wf_md["tag"] = tag;
        wf_md["pad"] = pad;
        wf_md["tbin"] = tbin0;
        wf_md["type"] = "waveform";
        itv->push_back(ITensor::pointer(st));

        auto& ch_md = cht->metadata();
        ch_md["tag"] = tag;
        ch_md["type"] = "channels";
        itv->push_back(ITensor::pointer(cht));

        if (have_summary) {
            auto& sum_md = sumt->metadata();
            sum_md["tag"] = tag;
            sum_md["type"] = "summary";
            itv->push_back(ITensor::pointer(sumt));
        }

        // per tag cmm
        for (auto label_cmm : in->masks()) {
            auto label = label_cmm.first;
            log->trace("masks: {}", label);
            auto cmm = label_cmm.second;
            size_t nranges = 0;
            for (auto ch_bl : cmm) {
                nranges += ch_bl.second.size();
            }
            SimpleTensor* ranges = new SimpleTensor({nranges, 2}, (double*)nullptr);
            Eigen::Map<Eigen::ArrayXXd> ranges_arr((double*) ranges->data(), nranges, 2);
            SimpleTensor* channels = new SimpleTensor({nranges}, (double*)nullptr);
            Eigen::Map<Eigen::ArrayXd> channels_arr((double*) channels->data(), nranges);

            size_t ind = 0;
            for (auto ch_bl : cmm) {
                auto ch = ch_bl.first;
                auto bl = ch_bl.second;
                for (auto b : bl) {
                    ranges_arr(ind, 0) = b.first;
                    ranges_arr(ind, 1) = b.second;
                    channels_arr(ind) = ch;
                    ++ind;
                }
            }

            auto& ranges_md = ranges->metadata();
            ranges_md["tag"] = tag;
            ranges_md["type"] = label + ":cmm_range";
            itv->push_back(ITensor::pointer(ranges));

            auto& channels_md = channels->metadata();
            channels_md["tag"] = tag;
            channels_md["type"] = label + ":cmm_channel";
            itv->push_back(ITensor::pointer(channels));
        }

        set_md["tags"].append(tag);
    }

    out = std::make_shared<SimpleTensorSet>(in->ident(), set_md, ITensor::shared_vector(itv));
    log->trace(Aux::dump(out));
    return true;
}

Aux::TaggedFrameTensorSet::TaggedFrameTensorSet() {}
Aux::TaggedFrameTensorSet::~TaggedFrameTensorSet() {}
