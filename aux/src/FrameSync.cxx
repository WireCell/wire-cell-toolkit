#include "WireCellAux/FrameSync.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(FrameSync,
                 WireCell::Aux::FrameSync,
                 WireCell::IFrameMerge)


using namespace WireCell;

Aux::FrameSync::FrameSync(const size_t multiplicity)
    : Aux::Logger("FrameSync", "glue")
    , m_multiplicity(multiplicity)
    , m_iqs(m_multiplicity) { }
Aux::FrameSync::~FrameSync() { }

std::vector<std::string> Aux::FrameSync::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

void Aux::FrameSync::flush(input_queues& iqs, output_queues& oqs)
{
    size_t neos=0, nmin=0, nempty=0;
    int min_ident=std::numeric_limits<int>::max();
    size_t min_index=0;

    const size_t nin = iqs.size();
    for (size_t ind=0; ind<nin; ++ind) {
        auto& iq = iqs[ind];
        if (iq.empty()) {
            ++nempty;
            log->debug("port {} empty", ind);
            continue;
        }
        auto frame = iq.front();
        if (!frame) {
            ++neos;
            log->debug("port {} eos", ind);
            continue;
        }
        log->debug("port {} frame {}", ind, frame->ident());
        /// FIXME: how about equal?
        if (min_ident > frame->ident()) {
            min_ident = frame->ident();
            min_index = ind;
            ++nmin;
        }
    }
    // log->debug("neos: {} nmin: {} nempty: {}", neos, nmin, nempty);

    // found min (good)
    if (nmin > 0) {
        auto frame = iqs[min_index].front();
        iqs[min_index].pop_front();
        std::get<0>(oqs).push_back(frame);

        // May have more
        return flush(iqs, oqs);
    }

    // eos sync
    if (neos == iqs.size()) {   // we have EOS sync
        for (auto& iq : iqs) {
            iq.pop_front();     // pop all EOS
        }
        std::get<0>(oqs).push_back(nullptr); // forward EOS

        // May have more behind the EOS
        return flush(iqs, oqs);
    }
}

bool Aux::FrameSync::operator()(input_queues& iqs, output_queues& oqs)
{
    // concatenate input queues to buffer
    for (size_t ind=0; ind<m_iqs.size(); ++ind) {
        m_iqs[ind].insert(m_iqs[ind].end(), iqs[ind].begin(), iqs[ind].end());
    }
    
    // try to flush
    flush(m_iqs, oqs);
    log->debug("oqs size: {}", std::get<0>(oqs).size());

    return true;
}
