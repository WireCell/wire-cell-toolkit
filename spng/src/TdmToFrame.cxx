#include "WireCellSpng/TdmToFrame.h"
#include "WireCellSpng/TdmFrame.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellSpng/Util.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellUtil/Waveform.h"

#include <regex>

namespace WireCell::SPNG {

    TdmToFrame::TdmToFrame()
        : Aux::Logger("TdmToFrame", "spng")
    {
    }

    TdmToFrame::~TdmToFrame()
    {
    }

    WireCell::Configuration TdmToFrame::default_configuration() const
    {
        Configuration cfg = this->ContextBase::default_configuration();
        cfg["frame"] = 0;
        return cfg;
    }

        
    void TdmToFrame::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);

        auto jff = cfg["frame"];
        if (jff.isInt()) {
            m_find_frame = jff.asInt(); // index
        }
        else if (jff.isString()) {
            m_find_frame = jff.asString(); // regex
        }
    }
    

    bool TdmToFrame::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            log->debug("EOS at call={}", m_count);
            ++m_count;
            return true;
        }

        TensorIndex ti(in);
        TdmFrame frame = std::holds_alternative<int>(m_find_frame)
            ? frame_at_index(ti, get<int>(m_find_frame))
            : frame_at_match(ti, get<std::string>(m_find_frame));
            
        if (!frame.frame) {
            log->critical("failed to get frame from tensor set id {} at call={}",
                          in->ident(), m_count);
            raise<ValueError>("failed to retrieve frame tensor from set.  Fix you configuration?");
        }


        Waveform::ChannelMaskMap cmm;
        for (const auto& [label, chmasks_itensor] : frame.chmasks) {
            if (!chmasks_itensor) {
                log->warn("no channel mask for label \"{}\"", label);
                continue;
            }

            auto cmten = chmasks_itensor->tensor();
            size_t nrows = cmten.size(0);
            Waveform::ChannelMasks cms;
            for (size_t irow=0; irow<nrows; ++irow) {
                auto cmvec = Torch::to_vector<int>(cmten[irow]);
                int batch = cmvec[0];
                if (batch) {
                    raise<ValueError>("TdmToFrame got batched channel masks tensor.  Precede this node with an unbatcher node");
                }
                int chid = cmvec[1];
                int beg = cmvec[2];
                int end = cmvec[3];
                cms[chid].emplace_back(beg,end);
            }
            cmm[label] = cms;
        }

        auto frame_md = frame.frame->metadata();
        // We'll load traces as we go
        auto all_traces = std::make_shared<ITrace::vector>();
        auto sf =  std::make_shared<Aux::SimpleFrame>(
            frame_md["ident"].asInt(),
            frame_md["time"].asFloat(),
            all_traces,
            frame_md["period"].asFloat(),
            cmm);

        /// Now actually fill in the traces and tagged traces.
        for (const auto& [tag, traces_itensor] : frame.traces) {
            if (!traces_itensor) {
                log->warn("no traces for tag \"{}\"", tag);
                continue;
            }
                
            auto traces_tensor = traces_itensor->tensor();
            const size_t ntraces = traces_tensor.size(0);
            if (!ntraces) {
                log->warn("empty traces tensor for tag \"{}\"", tag);
                continue;
            }

            auto chids_itensor = frame.chids[tag];
            if (!chids_itensor) {
                log->critical("frame corrupt: no chids tensor for tag \"{}\"", tag);
                raise<ValueError>("frame corrupt: no chids tensor");
            }
            auto chids_tensor = chids_itensor->tensor();
            const size_t nchannels = chids_tensor.size(0);
            if (nchannels != ntraces) {
                log->critical("frame corrupt: {} channels and {} traces", nchannels, ntraces);
                raise<ValueError>("frame corrupt: channel/trace count mismatch");
            }
            auto chids_vector = Torch::to_vector<int>(chids_tensor);
                
            auto traces_metadata = traces_itensor->metadata();
            const int tbin = traces_metadata["tbin"].asInt();
            
            const size_t beg = all_traces->size();
            for (size_t row=0; row < ntraces; ++row) {
                
                auto charge = Torch::to_vector<float>(traces_tensor);

                const int chid = chids_vector[row];
                auto itrace = std::make_shared<Aux::SimpleTrace>(chid, tbin, charge);
                all_traces->push_back(itrace);
            }
            const size_t end = all_traces->size();            
            
            IFrame::trace_list_t indices(end-beg);
            std::iota(indices.begin(), indices.end(), beg);

            auto summaries_itensor = frame.summaries[tag];
            if (summaries_itensor) {
                auto summaries_tensor = summaries_itensor->tensor();
                auto summaries = Torch::to_vector<double>(summaries_tensor);
                sf->tag_traces(tag, indices, summaries);
            }
            else {
                sf->tag_traces(tag, indices);
            }
        }

        out = sf;
        ++m_count;
        return true;
    }

}


