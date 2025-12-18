#include "WireCellSpng/TdmToFrame.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellSpng/Util.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h" // for taginfo()

#include "WireCellUtil/Waveform.h"

#include <regex>

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTdmToFrame,
                 WireCell::SPNG::TdmToFrame,
                 WireCell::INamed,
                 WireCell::IConfigurable,
                 WireCell::SPNG::ITorchSetToFrame)



namespace WireCell::SPNG {

    using HanaJsonCPP::to_json;
    using HanaJsonCPP::from_json;

    TdmToFrame::TdmToFrame()
        : Logger("TdmToFrame", "spng")
    {
    }
    TdmToFrame::TdmToFrame(const TdmToFrameConfig& cfg)
        : Logger("TdmToFrame", "spng")
        , m_cfg(cfg)
    {
        configme();
    }

    TdmToFrame::~TdmToFrame()
    {
    }

    WireCell::Configuration TdmToFrame::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = this->Logger::default_configuration();
        update(cfg, cfg2);
        cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

        
    void TdmToFrame::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);
        this->Logger::configure(cfg);
        from_json(m_cfg, cfg);
        configme();
    }

    static
    void assure_attribute(Configuration& cfg, const std::string& key, const std::string& value)
    {
        if (cfg.empty()) {
            return;
        }

        if (cfg.isArray()) {
            for (Configuration& one : cfg) {
                assure_attribute(one, key, value);
            }
            return;
        }
        if (!cfg.isObject()) {
            raise<ValueError>("assure_attribute requires object or array of object");
        }
        if (! cfg.isMember(key.c_str())) {
            cfg[key] = value;
        }
    }

    void TdmToFrame::configme()
    {
        // Give some user friendliness to not have to explicitly match on
        // datatype since user will already give this info as the attribute key.

        if (! m_cfg.frame.isObject()) {
            raise<ValueError>("single match object required for frame");
        }
        assure_attribute(m_cfg.frame, "datatype", "frame");

        for (auto& tt : m_cfg.tagged_traces) {
            if (tt.traces.empty()) {
                raise<ValueError>("no traces match object");
            }
            if (tt.chids.empty()) {
                raise<ValueError>("no chids match object");
            }
            assure_attribute(tt.traces, "datatype", "traces");
            assure_attribute(tt.chids, "datatype", "chids");
            assure_attribute(tt.summaries, "datatype", "summaries");
        }

        assure_attribute(m_cfg.chmasks, "datatype", "chmasks");
    }

    Waveform::ChannelMaskMap TdmToFrame::get_cmm(const ITorchTensor::vector& chmasks_itensors) const
    {
        Waveform::ChannelMaskMap cmm;
        for (const auto& chmasks_itensor : chmasks_itensors) {
            auto md = chmasks_itensor->metadata();
            if (md["label"].empty()) {
                log->warn("no label for channel mask, corrupt TDM, skipping");
                continue;
            }
            std::string label = md["label"].asString();

            auto cmten = chmasks_itensor->tensor().to(torch::kCPU); // do I need the semaphore here?
            size_t nrows = cmten.size(0);
            Waveform::ChannelMasks cms;
            for (size_t irow=0; irow<nrows; ++irow) {
                auto cmvec = to_vector<int>(cmten[irow]);
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
        return cmm;
    }

    std::vector<int> TdmToFrame::get_chids(const ITorchTensor::vector& chids_itensors) const
    {
        std::vector<int> all_chids;
        for (const auto& chids_itensor : chids_itensors) {
            auto chv = to_vector<int>(chids_itensor->tensor());
            all_chids.insert(all_chids.end(), chv.begin(), chv.end());
        }
        log->debug("have {} channels from {} groups", all_chids.size(), chids_itensors.size());
        return all_chids;
    }

    std::vector<double> TdmToFrame::get_summaries(const ITorchTensor::vector& summaries_itensors) const
    {
        std::vector<double> all_summaries;
        for (const auto& summaries_itensor : summaries_itensors) {
            auto summaries = to_vector<double>(summaries_itensor->tensor());
            all_summaries.insert(all_summaries.end(), summaries.begin(), summaries.end());
        }
        return all_summaries;
    }

    ITrace::vector TdmToFrame::get_traces(const ITorchTensor::vector& traces_itensors,
                                          const std::vector<int>& chids) const
    {
        log->debug("have {} traces itensors and {} chids",
                   traces_itensors.size(), chids.size());

        ITrace::vector all_traces;

        size_t chid_index = 0;
        for (const auto& traces_itensor : traces_itensors) {
            auto traces_tensor = traces_itensor->tensor().to(torch::kCPU);
            const size_t ntraces = traces_tensor.size(0);

            auto traces_metadata = traces_itensor->metadata();
            const int tbin = traces_metadata["tbin"].asInt();
                
            for (size_t row=0; row < ntraces; ++row) {
                auto charge = to_vector<float>(traces_tensor[row]);
                const int chid = chids[chid_index++];
                auto itrace = std::make_shared<Aux::SimpleTrace>(chid, tbin, charge);
                all_traces.push_back(itrace);
            }
        }
        if (all_traces.size() != chids.size()) {
            log->critical("trace / channel count mismatch: {} != {}, inconsistent selection criteria?", all_traces.size(), chids.size());
            raise<ValueError>("trace / channel count mismatch: %d != %d", all_traces.size(), chids.size());
        }
        return all_traces;
    }

    bool TdmToFrame::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "input");

        // partition to help speed up matching a bit.
        auto tensors_by_datatype = TDM::by_datatype(*in->tensors());

        // Get THE frame tensor, complain if we do not have exactly 1. 
        auto frame_itensors = tensors_by_datatype["frame"];
        log->debug("input {} frame tensors", frame_itensors.size());
        frame_itensors = TDM::select_tensors(frame_itensors, m_cfg.frame);
        // log->debug("selected {} frame tensors with {}", frame_itensors.size(), m_cfg.frame);
        if (frame_itensors.empty()) {
            log->warn("no frame tensor, returning empty frame for set ident {} at call={}",
                      in->ident(), get_count());
            out =  std::make_shared<Aux::SimpleFrame>(in->ident());
            next_count();
            return true;
        }
        if (frame_itensors.size() > 1) {        
            log->warn("multiple frame tensors, using first in ident {} at call={}",
                      in->ident(), get_count());
            // fixme: if multi-frame becomes a thing, TdmToFrame can be changed
            // to a queued out node.  Or, we make a plural TdmToFrames and an
            // IFrameSet.
        };
        auto frame_itensor = frame_itensors[0];


        // okay, down to business.


        Waveform::ChannelMaskMap cmm = get_cmm(
            TDM::select_tensors(tensors_by_datatype["chmasks"], m_cfg.frame));

        auto frame_md = frame_itensor->metadata();
        // We'll load traces as we go
        auto all_traces = std::make_shared<ITrace::vector>();
        auto sf =  std::make_shared<Aux::SimpleFrame>(
            frame_md["ident"].asInt(),
            frame_md["time"].asFloat(),
            all_traces,
            frame_md["period"].asFloat(),
            cmm);

        /// Loop over each requested tagged traces set.
        for (const auto& tt_cfg : m_cfg.tagged_traces) {

            auto chids_vector = get_chids(
                TDM::select_tensors(tensors_by_datatype["chids"], tt_cfg.chids));
            
            auto traces_vector = get_traces(
                TDM::select_tensors(tensors_by_datatype["traces"], tt_cfg.traces),
                chids_vector);

            // Form trace indices for this block of tagged traces.
            const size_t traces_beg = all_traces->size();
            all_traces->insert(all_traces->end(), traces_vector.begin(), traces_vector.end());
            IFrame::trace_list_t indices(traces_vector.size());
            std::iota(indices.begin(), indices.end(), traces_beg);

            auto summaries_vector = get_summaries(
                TDM::select_tensors(tensors_by_datatype["summaries"], tt_cfg.summaries));

            // summaries can be empty.
            sf->tag_traces(tt_cfg.tag, indices, summaries_vector);

        }

        out = sf;
        log->debug("call={} output frame: {}", get_count(), Aux::taginfo(out));

        next_count();
        return true;
    }

}


