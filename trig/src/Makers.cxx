#include "WireCellTrig/Makers.h"

using namespace WireCell::Trig;

namespace {
    // Prescale and strip inputs of the objects in out[first:].
    template <typename Vec>
    void postprocess_range(Vec& out, size_t first, const MakerConfig& cfg, uint64_t& count)
    {
        if (cfg.prescale > 1) {
            auto iter = out.begin() + first;
            while (iter != out.end()) {
                ++count;
                if (count % cfg.prescale != 0) {
                    iter = out.erase(iter);
                    continue;
                }
                ++iter;
            }
        }
        if (!cfg.store_inputs) {
            for (size_t ind = first; ind < out.size(); ++ind) {
                out[ind].inputs.clear();
            }
        }
    }
}  // namespace

TAMaker::~TAMaker() {}

bool TAMaker::preprocess(const input_type& tp)
{
    return tp.samples_over_threshold <= m_common.max_samples_over_threshold;
}

void TAMaker::do_flush(Daq::timestamp_t, output_vector&) {}

void TAMaker::operator()(const input_type& tp, output_vector& out)
{
    if (!preprocess(tp)) return;
    const size_t first = out.size();
    process(tp, out);
    postprocess(out, first);
}

void TAMaker::flush(Daq::timestamp_t until, output_vector& out)
{
    const size_t first = out.size();
    do_flush(until, out);
    postprocess(out, first);
}

void TAMaker::postprocess(output_vector& out, size_t first) { postprocess_range(out, first, m_common, m_count); }

TCMaker::~TCMaker() {}

bool TCMaker::preprocess(const input_type&) { return true; }

void TCMaker::do_flush(Daq::timestamp_t, output_vector&) {}

void TCMaker::operator()(const input_type& ta, output_vector& out)
{
    if (!preprocess(ta)) return;
    const size_t first = out.size();
    process(ta, out);
    postprocess(out, first);
}

void TCMaker::flush(Daq::timestamp_t until, output_vector& out)
{
    const size_t first = out.size();
    do_flush(until, out);
    postprocess(out, first);
}

void TCMaker::postprocess(output_vector& out, size_t first) { postprocess_range(out, first, m_common, m_count); }
