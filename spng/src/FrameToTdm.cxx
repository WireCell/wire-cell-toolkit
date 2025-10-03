#include "WireCellSpng/FrameToTdm.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Fmt.h"

WIRECELL_FACTORY(SPNGFrameToTdm,
                 WireCell::SPNG::FrameToTdm,
                 WireCell::INamed,
                 WireCell::IConfigurable,
                 WireCell::SPNG::IFrameToTorchSet)


namespace WireCell::SPNG {

    FrameToTdm::FrameToTdm()
        : Logger("FrameToTdm", "spng")
    {
    }
 
   FrameToTdm::~FrameToTdm()
   {
   }

    WireCell::Configuration FrameToTdm::default_configuration() const
    {
        Configuration cfg = this->Logger::default_configuration();
        
        cfg["basepath"] = m_basepath.string();
        cfg["frame_relpath"] = m_frame_relpath.string();

        // And many more.  See comments in header and spng/docs/frametotdm.org.
        return cfg;
    }
    
    void FrameToTdm::configure(const WireCell::Configuration& cfg)
    {
        this->Logger::configure(cfg);

        std::string anode = cfg["anode"].asString();
        m_anode = Factory::find_tn<IAnodePlane>(anode);


        // build global channel order.  FIXME: we may not actually need to keep
        // m_anode after this....
        m_channel_order.clear();
        for (const int chid : m_anode->channels()) {
            auto ichan = m_anode->channel(chid);
            const int wpid = ichan->planeid().ident();
            const size_t wan = ichan->index();
            m_channel_order[chid] = ((size_t)wpid << 32) & wan;
            m_wpid_channels[wpid].insert(chid);
        }

        m_basepath = get<std::string>(cfg, "basepath", m_basepath.string());
        m_frame_relpath = get<std::string>(cfg, "frame_relpath", m_frame_relpath.string());

        m_chmasks.clear();
        m_chmasks[""] = "chmasks/{label}"; // default relpath
        auto jchmasks = cfg["chmasks"];
        if (jchmasks.isObject()) {
            for (const auto& label : jchmasks.getMemberNames()) {
                m_chmasks[label] = jchmasks[label].asString();
            }
        }

        m_rules.clear();
        auto jrules = cfg["rules"];
        if (! jrules.isArray()) {
            return;
        }

        // WCT config wrangling is exhausting.
        for (const auto& jrule : jrules) {
            auto jgroups = jrule["groups"];
            if (! jgroups.isArray()) {
                log->warn("skipping rule config with no groups: {}", jrule);
                continue;
            }
            
            std::string tag = get<std::string>(jrule, "tag");
            if (tag == "null" or tag == "*") {
                tag = "";
            }
            Rule rule{tag};
            for (const auto& jgroup : jgroups) {
                Group group;
                auto relpath = get<std::string>(jgroup, "relpath");
                if (relpath.size()) {
                    group.relpath = relpath;
                }
                auto wpids = get<std::vector<int>>(jgroup, "wpids", {});
                for (int wpid : wpids) {
                    group.wpids.insert(wpid);
                }
                for (const auto& chid : get<std::vector<int>>(jgroup, "channels", {})) {
                    group.channels.insert(chid);
                }

                if (group.wpids.empty() && group.channels.empty()) {
                    log->warn("neither wpids nor channels for tag \"{}\", group at relpath {}",
                              tag, relpath);
                    continue;
                }

                rule.groups.push_back(std::move(group));
            }
            if (rule.groups.empty()) {
                log->warn("skipping rule with no groups: {}", jrule);
            }
            m_rules.push_back(std::move(rule));
        }
        if (m_rules.empty()) {
            // Maybe this is something user actually wants, but unlikely.
            log->warn("no frame conversion rules");
        }



    }
    
    bool FrameToTdm::operator()(const input_pointer& inframe, output_pointer& outtens) 
    {
        outtens = nullptr;
        if (!inframe) {
            logit("EOS");
            ++m_count;
            return true;
        }

        auto frame = frame_itensor(inframe);
        const std::string parent = frame->metadata()["datapath"].asString();

        auto tensors = std::make_shared<ITorchTensor::vector>();
        tensors->push_back(frame);

        ITorchTensor::vector rules = rules_tensors(inframe, parent);
        if (rules.size()) {
            tensors->insert(tensors->end(), rules.begin(), rules.end());
        }

        ITorchTensor::vector chmasks = chmask_tensors(inframe, parent);
        if (chmasks.size()) {
            tensors->insert(tensors->end(), chmasks.begin(), chmasks.end());
        }

        Configuration empty;
        outtens = std::make_shared<SimpleTorchTensorSet>(inframe->ident(), empty, tensors);

        logit(outtens, "output");

        ++m_count;
        return true;
    }



    ITorchTensor::pointer FrameToTdm::frame_itensor(const IFrame::pointer& iframe) const
    {
        Configuration md;
        int ident = iframe->ident();
        md["ident"] = ident;
        md["time"] = iframe->time();
        md["period"] = iframe->tick();
        md["datatype"] = "frame";
        auto fullpath = m_basepath / m_frame_relpath;
        auto fullpathstr = fmt::format(fullpath.string(), fmt::arg("ident", ident));
        log->debug("frame tensor datapath: {}", fullpathstr);
        md["datapath"] = fullpathstr;
        // no parent, no batches
        return std::make_shared<SimpleTorchTensor>(md);
    }


    ITorchTensor::vector FrameToTdm::chmask_tensors(const IFrame::pointer& iframe, const std::string& parent) const
    {
        const auto chmasks = iframe->masks();
        const int ident = iframe->ident();
        ITorchTensor::vector tens;

        for (const auto& [label, cms] : chmasks) {

            auto ten = chmask_tensor(cms);

            std::string relpath = "";
            for (const auto& maybe : std::vector<std::string>({label, ""})) {
                auto it = m_chmasks.find(maybe);
                if (it == m_chmasks.end()) {
                    continue;
                }
                relpath = it->second;
                break;
            }

            Configuration md;
            md["ident"] = ident;
            md["label"] = label;
            md["parent"] = parent;
            tens.push_back(make_datatype("chmasks", relpath, ten, md));
        }
        return tens;
    }

    torch::Tensor FrameToTdm::chmask_tensor(const Waveform::ChannelMasks& cms) const
    {


        /* We must convert this fiddly structure
           typedef std::pair<int, int> BinRange;
           typedef std::vector<BinRange> BinRangeList;
           typedef std::map<int, BinRangeList> ChannelMasks;

           into 2D tensor shape (Nbinranges, 4) holding columns
           (batch,chid,beg,end).  The "batch" is always 0 as we emit unbatched
           tensors.  A later batcher node may do batching.
        */

        std::vector<torch::Tensor> tens;

        const int64_t batch = 0;
        for (const auto& [chid_int, brl] : cms) {
            const int64_t chid = chid_int;
            for (const auto& br : brl) {
                const auto& [beg,end] = br;
                auto ten = torch::tensor({batch, chid, (int64_t)beg, (int64_t)end},
                                         torch::kInt32);
                tens.emplace_back(ten);
            }
        }
        return torch::vstack(tens);
        // Caller handles metadata.
    }

        
    // Return trace indices for tag.  
    std::vector<size_t> FrameToTdm::tag_indices(const IFrame::pointer& iframe, const std::string& tag) const
    {
        // FIXME: this could go into aux's FrameTools.h
        std::vector<size_t> inds;
        if (tag.empty()) { // all
            const size_t ntraces = iframe->traces()->size();
            inds.resize(ntraces);
            std::iota(inds.begin(), inds.end(), 0);
        }
        else {
            inds = iframe->tagged_traces(tag);
        }
        return inds;
    }

    /// Return subset of channel IDs from have that are consistent with group.
    /// Sort according to global channel order.
    std::vector<int> FrameToTdm::group_channels(const std::vector<int>& have,
                                                const Group& group) const
    {
        std::set<int> want(group.channels.begin(), group.channels.end());
        for (const auto& wpid : group.wpids) {
            auto wit = m_wpid_channels.find(wpid);
            if (wit == m_wpid_channels.end()) {
                continue;
            }
            for (int chid : wit->second) {
                want.insert(chid);
            }
        }

        std::vector<int> got;
        std::set_intersection(have.begin(), have.end(),
                              want.begin(), want.end(),
                              std::inserter(got, got.begin()));

        std::sort(got.begin(), got.end(), [&](int a, int b) {
            auto ait = m_channel_order.find(a);
            auto bit = m_channel_order.find(b);
            auto eit = m_channel_order.end();
            if (ait != eit && bit != eit) {
                return ait->second < bit->second;
            }
            if (ait != eit) {
                return true;
            }
            return false;
        });
        return got;
    }


    
    ITorchTensor::pointer FrameToTdm::make_datatype(
        const std::string& datatype, const boost::filesystem::path& relpath,
        torch::Tensor ten, Configuration md) const
    {
        md["datatype"] = datatype;
        auto fullpath = m_basepath / relpath;
        // log->debug("making datatype {} with fullpath {}", datatype, fullpath.string());
        try {
            md["datapath"] = Fmt::format(fullpath.string(), md);
        }
        catch (const fmt::format_error& err) {
            log->critical("failed to format for datatype \"{}\" using {} with {}",
                          datatype, fullpath.string(), md);
            raise<ValueError>("failed to format tensor path, bad config?");
        }

        return std::make_shared<SimpleTorchTensor>(ten, md);
    }


    /// Run through the rules
    ITorchTensor::vector FrameToTdm::rules_tensors(
        const IFrame::pointer& iframe,
        const std::string& parent) const
    {
        ITorchTensor::vector tensors;

        const ITrace::vector& all_traces = *(iframe->traces());
        const int frame_ident = iframe->ident();


        // std::unordered_map<int, std::vector<ITrace::pointer>> channel_traces;
        // std::unordered_map<int, size_t> chid_index;

        /// For each rules we need:
        /// 1. set of trace indices specific to tag and wpids/channels config.
        /// 2. map from channel ID to tensor row by the global sort.
        /// 3. iterate over indices, use map to find row, store.

        int rule_index = -1;
        for (const auto& rule : m_rules) {
            ++rule_index;
            std::string tag = rule.tag;
            
            // The trace indices relevant to the tag.
            auto tag_inds = tag_indices(iframe, tag);
            std::vector<int> tag_chids;
            for (size_t ind : tag_inds) {
                tag_chids.push_back(all_traces[ind]->channel());
            }

            int group_index = -1;
            for (const auto& group : rule.groups) {
                ++group_index;

                // Vector of output ordered channel IDs
                auto ordered_chids = group_channels(tag_chids, group);

                // Map chid to its eventual 
                std::unordered_map<int, long> chid_row;
                for (auto chid : ordered_chids) {
                    chid_row[chid] = chid_row.size();
                }

                struct TIR {
                    ITrace::pointer trace;
                    size_t index; // into a summary or channels vector
                    int row;   // tensor row index
                };
                std::vector<TIR> tirs;

                // Find subset of indices with channels in grp_chids.
                // And, find traces bounds.
                std::vector<size_t> grp_inds;
                int tbeg=0, tend=0, nrows=0;
                size_t sindex=0; // the index into a summary
                for (size_t ind : tag_inds) {
                    const auto& trace = all_traces[ind];
                    int chid = trace->channel();
                    
                    auto it = chid_row.find(chid);
                    if (it == chid_row.end()) {
                        // a tagged trace but its channel is not in our group.
                        ++sindex;
                        continue;
                    }
                    
                    int row = it->second;
                    tirs.emplace_back(TIR{trace, sindex, row});
                    ++sindex;
                    
                    int tbin = trace->tbin();
                    int nbins = trace->charge().size();
                    
                    ++nrows;
                    if (nrows == 1) { // first time
                        tbeg = tbin;
                        tend = tbin+nbins;
                        continue;
                    }
                    if (tbeg < tbin) { tbeg = tbin; }
                    if (tend < tbin+nbins) { tend = tbin+nbins;}
                }
                int ncols = tend-tbeg;
  
                // Common md for all parts: traces, summaries, chids.
                Configuration common_md;
                common_md["ident"] = frame_ident;
                common_md["rule"] = rule_index;
                common_md["group"] = group_index;
                common_md["parent"] = parent;
                common_md["tag"] = tag.empty() ? "null" : tag;

                // Build traces tensor
                {
                    torch::Tensor ten = torch::zeros({nrows, ncols}, torch::kFloat32);
                    for (const auto& tir : tirs) {
                        const auto& charge = tir.trace->charge();
                        int tbin = tir.trace->tbin();
                        int col_beg = tbin - tbeg;
                        int col_end = col_beg + charge.size();
                        auto tmp = torch::tensor(charge, torch::kFloat32);
                        ten.index({tir.row, torch::indexing::Slice(col_beg, col_end)}) += tmp;
                    }

                    // traces have extra md
                    Configuration md = common_md;
                    md["tbin"] = (int)tbeg;
                    md["period"] = iframe->tick();
                    tensors.push_back(make_datatype("traces", group.relpath, ten, md));
                }

                // Build summaries tensor
                auto summary = iframe->trace_summary(tag);
                if (summary.size()) {
                    torch::Tensor ten = torch::zeros({nrows}, torch::kFloat64);
                    for (const auto& tir : tirs) {
                        const double value = summary[tir.index];
                        ten.index({tir.row}) += value;
                    }
                    tensors.push_back(make_datatype("summaries", group.relpath, ten, common_md));
                }

                // Build chids tensor
                {
                    torch::Tensor ten = torch::zeros({nrows}, torch::kInt);
                    for (const auto& tir : tirs) {
                        int chid = tir.trace->channel();
                        ten.index_put_({tir.row}, chid); // no accumulate
                    }
                    tensors.push_back(make_datatype("chids", group.relpath, ten, common_md));
                }                
            } // groups
        } // rules
        return tensors;
    } // rules_tensor


}
        
        

