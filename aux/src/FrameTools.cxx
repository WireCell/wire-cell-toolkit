#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <numeric>
#include <sstream>

using namespace WireCell;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

std::string Aux::name(const WireCell::IFrame::pointer& frame)
{
    if (!frame) {
        return "(null frame)";
    }

    std::stringstream ss;
    ss << "frame_" << frame->ident();
    return ss.str();
    
}

int Aux::frmtcmp(IFrame::pointer frame, double time)
{
    auto traces = frame->traces();
    auto tbin_mm = Aux::tbin_range(*traces.get());

    const double tref = frame->time();
    const double tick = frame->tick();
    const double tmin = tref + tbin_mm.first * tick;   // low edge
    const double tmax = tref + tbin_mm.second * tick;  // high edge

    if (tmax <= time) {  // high edge may be exactly equal to target time and frame does not span.
        return -1;
    }
    if (tmin >= time) {  // low edge may be exactly equal to target time and frame does not span.
        return +1;
    }
    return 0;
}
std::pair<IFrame::pointer, IFrame::pointer> Aux::split(IFrame::pointer frame, double time)
{
    if (!frame) {
        return std::pair<IFrame::pointer, IFrame::pointer>(nullptr, nullptr);
    }
    int cmp = frmtcmp(frame, time);
    if (cmp < 0) {
        return std::pair<IFrame::pointer, IFrame::pointer>(frame, nullptr);
    }
    if (cmp > 0) {
        return std::pair<IFrame::pointer, IFrame::pointer>(nullptr, frame);
    }
    // now gotta do the hard work.

    const double tref = frame->time();
    const double tick = frame->tick();
    const int ident = frame->ident();

    // Every tick equal or larger than this is in the second frame.
    const double fnticks = (time - tref) / tick;
    const int tbin_split = 0.5 + fnticks;

    ITrace::vector mtraces, ptraces;
    for (auto trace : (*frame->traces())) {
        const int tbin = trace->tbin();

        if (tbin >= tbin_split) {
            ptraces.push_back(trace);
            continue;
        }

        const ITrace::ChargeSequence& wave = trace->charge();
        const int lbin = tbin + wave.size();

        if (lbin <= tbin_split) {
            mtraces.push_back(trace);
            continue;
        }
        const int ind_split = tbin_split - tbin;
        ITrace::ChargeSequence mcharge(wave.begin(), wave.begin() + ind_split);
        ITrace::ChargeSequence pcharge(wave.begin() + ind_split, wave.end());

        const int chid = trace->channel();

        auto mtrace = std::make_shared<SimpleTrace>(chid, tbin, mcharge);
        mtraces.push_back(mtrace);
        auto ptrace = std::make_shared<SimpleTrace>(chid, tbin_split, pcharge);
        ptraces.push_back(ptrace);
    }
    // fixme: what about ident, cmm, tags....
    IFrame::pointer mframe = std::make_shared<SimpleFrame>(ident, tref, mtraces, tick);
    IFrame::pointer pframe = std::make_shared<SimpleFrame>(ident, tref, ptraces, tick);
    return std::pair<IFrame::pointer, IFrame::pointer>(mframe, pframe);
}

ITrace::vector Aux::untagged_traces(IFrame::pointer frame)
{
    ITrace::vector ret;
    if (!frame) {
        return ret;
    }
    auto traces = frame->traces();
    size_t ntraces = traces->size();

    std::vector<size_t> tagged;
    for (auto tag : frame->trace_tags()) {
        const auto& taglist = frame->tagged_traces(tag);
        tagged.insert(tagged.end(), taglist.begin(), taglist.end());
    }
    std::sort(tagged.begin(), tagged.end());
    std::vector<size_t> all(ntraces), untagged;
    std::iota(all.begin(), all.end(), 0);
    std::set_difference(all.begin(), all.end(), tagged.begin(), tagged.end(),
                        std::inserter(untagged, untagged.begin()));
    for (size_t ind : untagged) {
        ret.push_back(traces->at(ind));
    }
    return ret;
}

ITrace::vector Aux::tagged_traces(IFrame::pointer frame, IFrame::tag_t tag)
{
    ITrace::vector ret;
    if (!frame) {
        return ret;
    }
    if (tag == "") {
        return untagged_traces(frame);
    }
    auto const& all_traces = frame->traces();
    for (size_t index : frame->tagged_traces(tag)) {
        ret.push_back(all_traces->at(index));
    }
    if (!ret.empty()) {
        return ret;
    }
    auto ftags = frame->frame_tags();
    if (std::find(ftags.begin(), ftags.end(), tag) == ftags.end()) {
        return ret;
    }
    return *all_traces;  // must make copy of shared pointers
}

Aux::channel_list Aux::channels(const ITrace::vector& traces, bool uniq)
{
    const auto nchans = traces.size();
    Aux::channel_list ret(nchans, 0);
    for (size_t ind = 0; ind != nchans; ++ind) {
        ret[ind] = traces[ind]->channel();
    }
    if (uniq) {
        std::sort(ret.begin(), ret.end());
        auto end = std::unique(ret.begin(), ret.end());
        ret.resize(std::distance(ret.begin(), end));
    }

    return ret;
}

std::pair<int, int> Aux::tbin_range(const ITrace::vector& traces)
{
    const auto siz = traces.size();
    std::vector<int> tbins(siz), tlens(siz);
    for (size_t ind = 0; ind != siz; ++ind) {
        const auto& trace = traces[ind];
        const int tbin = trace->tbin();
        tbins[ind] = tbin;
        tlens[ind] = tbin + trace->charge().size();
    }
    return std::pair<int, int>(*std::min_element(tbins.begin(), tbins.end()),
                               *std::max_element(tlens.begin(), tlens.end()));
}

void Aux::fill(Array::array_xxf& array, const ITrace::vector& traces, channel_list::iterator chit,
               channel_list::iterator chend, int tbin)
{
    std::unordered_map<int, int> index;
    // one col is one tick
    // one row is one channel
    // array is indexed in order: (irow, icol)
    const int ncols = array.cols();
    const int nrows = std::min((int) array.rows(), (int) std::distance(chit, chend));
    for (int ind = 0; ind != nrows and chit != chend; ++ind, ++chit) {
        index[*chit] = ind;
    }
    for (const auto& trace : traces) {
        // resolve which row a the channel is at
        const int ch = trace->channel();
        auto it = index.find(ch);
        if (it == index.end()) {
            continue;
        }
        const int irow = it->second;

        const auto& charge = trace->charge();
        const int nticks = charge.size();
        const int dtbin = trace->tbin() - tbin;
        int icol0 = 0, itick0 = 0;
        if (dtbin < 0) {
            itick0 = -dtbin;
        }
        if (dtbin > 0) {
            icol0 = dtbin;
        }

        const int ncols_left = ncols - icol0;
        const int nticks_left = nticks - itick0;
        if (ncols_left <= 0 or nticks_left <= 0) {
            continue;
        }
        const int nleft = std::min(ncols_left, nticks_left);
        for (int ind = 0; ind != nleft; ++ind) {
            array(irow, icol0 + ind) += charge.at(itick0 + ind);
        }
    }
}

Aux::channel_list Aux::fill(Array::array_xxf& array,
                            const ITrace::vector& traces)
{
    auto ch = Aux::channels(traces, true);
    auto tbinmm = Aux::tbin_range(traces);
    const size_t ncols = tbinmm.second - tbinmm.first;
    const size_t nrows = std::distance(ch.begin(), ch.end());
    Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols);
    Aux::fill(arr, traces, ch.begin(), ch.end(), tbinmm.first);
    return ch;
}


std::string Aux::taginfo(const WireCell::IFrame::pointer& frame)
{
    if (! frame) {
        return "(null frame pointer)";
    }

    std::stringstream info;
    info << "frame: ident=" << frame->ident()
         << " time=" << frame->time() << " tick=" << frame->tick()
         << " with ";
    const auto& traces = frame->traces();
    if(!traces) {
        info << " null trace ptr!";
    } else {
        info << traces->size() << " traces.  frame tags:[ ";
    }
    for (const auto& tag : frame->frame_tags()) {
        info << "\"" << tag << "\" ";
    }

    auto ttags = frame->trace_tags();
    info << "] " << ttags.size() << " tagged trace sets:[ ";
    
    for (const auto& tag : ttags) {
        const auto& taglist = frame->tagged_traces(tag);
        info << "\"" << tag << "\":" << taglist.size() << " ";
        const auto& summary = frame->trace_summary(tag);
        info << "[" << summary.size() << "] ";
    }
    info << "]";

    info << " cmm:[ ";
    for (const auto& [name,cm] : frame->masks()) {
        info << name << ":" << cm.size() << " ";
    }
    info << "]";
    return info.str();
}


void Aux::dump_frame(WireCell::IFrame::pointer frame, Log::logptr_t log)
{
    std::stringstream info;

    auto traces = frame->traces();
    const size_t ntraces = traces->size();
    std::vector<double> means, rmses, lengths, tbins;
    for (auto trace : *traces) {
        auto const& charge = trace->charge();

        auto mr = Waveform::mean_rms(charge);
        means.push_back(mr.first);
        rmses.push_back(mr.second * mr.second);
        const int nsamps = charge.size();
        lengths.push_back(nsamps);
        tbins.push_back(trace->tbin());

        if (std::isnan(mr.second)) {
            info << "Frame: channel " << trace->channel() << " rms is NaN\n";
        }

        for (int ind = 0; ind < nsamps; ++ind) {
            float val = charge[ind];
            if (std::isnan(val)) {
                info << "Frame: channel " << trace->channel() << " sample " << ind << " is NaN\n";
            }
            if (std::isinf(val)) {
                info << "Frame: channel " << trace->channel() << " sample " << ind << " is INF\n";
            }
        }
    }
    double meanmean = Waveform::sum(means) / ntraces;
    double totrms = sqrt(Waveform::sum(rmses));
    double meanlen = Waveform::sum(lengths) / ntraces;
    double meantbin = Waveform::sum(tbins) / ntraces;
    info << "Frame: " << ntraces << " traces,"
         << " <mean>=" << meanmean << " TotRMS=" << totrms << " <tbin>=" << meantbin << " <len>=" << meanlen
         << "\n";
    for (auto it : frame->masks()) {
        info << "\t" << it.first << " : " << it.second.size() << "\n";
    }
    if (log) {
        log->debug(info.str());
    }
    else {
        std::cerr << info.str();
    }
}

void Aux::raster(WireCell::Array::array_xxf& block, WireCell::ITrace::vector traces,
                 const std::vector<int>& channels)
{
    const size_t nchannels = channels.size();
    std::unordered_map<int, size_t> ch2col;
    for (size_t ind = 0; ind < nchannels; ++ind) {
        ch2col[channels[ind]] = ind;
    }

    const size_t ncols = block.cols();
    for (auto trace : traces) {
        const auto& samples = trace->charge();
        const size_t nsamples = samples.size();
        const size_t tbin = trace->tbin();

        if (tbin >= ncols) {  // underflow impossible as they are unsigned.
            continue;         // trace is off screen
        }

        // find the row to fill
        const int ch = trace->channel();
        const auto chit = ch2col.find(ch);
        if (chit == ch2col.end()) {
            continue;
        }
        const size_t irow = chit->second;
        WireCell::Array::array_xf tofill(ncols - tbin);
        const size_t maxbin = std::min(nsamples, ncols - tbin);
        for (size_t bin = 0; bin < maxbin; ++bin) {
            tofill[bin] = samples[bin];
        }
        block.block(irow, tbin, 1, ncols - tbin).row(0) += tofill;
    }
}

typedef std::pair<int, const std::vector<float>*> tbin_charge_t;
typedef std::map<int, std::vector<tbin_charge_t> > channel_index_t;

IFrame::pointer Aux::sum(std::vector<IFrame::pointer> frames, int ident)
{
    // Extract starting times and ticks of all frames
    const int nframes = frames.size();
    std::vector<double> times(nframes), ticks(nframes);
    for (int ind = 0; ind < nframes; ++ind) {
        times[ind] = frames[ind]->time();
        ticks[ind] = frames[ind]->tick();
    }
    auto tick_mm = std::minmax_element(ticks.begin(), ticks.end());
    if (*tick_mm.first != *tick_mm.second) {
        return nullptr;
    }
    const double tick = ticks[0];
    auto times_mm = std::minmax_element(times.begin(), times.end());
    const double start_time = *times_mm.first;

    // Make an temporary index mapping channel to all the trace data,
    // applying time offsets as we fill.
    channel_index_t index;
    for (int ind = 0; ind < nframes; ++ind) {
        IFrame::pointer frame = frames[ind];
        const double dt = frame->time() - start_time;
        const int tbinoff = dt / tick;
        for (auto trace : *frame->traces()) {
            const int ch = trace->channel();
            const int new_tbin = tbinoff + trace->tbin();
            index[ch].push_back(make_pair(new_tbin, &(trace->charge())));
        }
    }

    // Process each channel to make flattened trace.
    ITrace::vector out_traces;
    for (auto ctc : index) {
        auto ch = ctc.first;
        auto tcv = ctc.second;  // vector of tbin_charge_t
        std::vector<int> tbins;
        for (auto tc : tcv) {  // go through once to find tbin bounds
            const int tbin = tc.first;
            tbins.push_back(tbin);
            tbins.push_back(tbin + tc.second->size());
        }
        auto tbin_mm = std::minmax_element(tbins.begin(), tbins.end());
        const int trace_tbin = *tbin_mm.first;
        const int trace_nbins = *tbin_mm.second - trace_tbin;
        std::vector<float> charge(trace_nbins, 0.0);
        for (auto tc : tcv) {
            int ind = tc.first - trace_tbin;
            for (auto q : *tc.second) {
                charge[ind] += q;
                ++ind;
            }
        }
        ITrace::pointer trace = std::make_shared<SimpleTrace>(ch, trace_tbin, charge);
        out_traces.push_back(trace);
    }

    return std::make_shared<SimpleFrame>(ident, start_time, out_traces, tick);
}

