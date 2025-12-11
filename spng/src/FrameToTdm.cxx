#include "WireCellSpng/FrameToTdm.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/TensorTools.h"
#include "WireCellSpng/Util.h"
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
        Configuration cfg2 = this->ContextBase::default_configuration();
        update(cfg, cfg2);

        // And many more.  See comments in header and spng/docs/frametotdm.org.
        cfg["basepath"] = m_basepath.string();
        cfg["frame_relpath"] = m_frame_relpath.string();

        return cfg;
    }
    
    /// Return IChannels spanning wpid in anode in order of increasing
    /// IChannel::index().  If wpid_num is negative, reverse the order.
    static IChannel::vector get_ordered_channels(IAnodePlane::pointer anode, int wpid_num)
    {
        WirePlaneId wpid(std::abs(wpid_num));
        auto face = anode->face(wpid.face());
        auto plane = face->planes()[wpid.index()];
        IChannel::vector chans = plane->channels(); // already ordered by IChannel::index().
        
        if (wpid_num < 0) {
            std::reverse(chans.begin(), chans.end());
        }
        return chans;
    }

    void FrameToTdm::configure(const WireCell::Configuration& cfg)
    {
        this->Logger::configure(cfg);
        this->ContextBase::configure(cfg);

        std::string anode_tn = cfg["anode"].asString();
        auto anode = Factory::find_tn<IAnodePlane>(anode_tn);

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
        int rule_num = -1;
        for (const auto& jrule : jrules) {
            ++rule_num;

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
            int group_num = -1;
            for (const auto& jgroup : jgroups) {
                ++group_num;

                Group group;
                auto relpath = get<std::string>(jgroup, "relpath");
                if (relpath.size()) {
                    group.relpath = relpath;
                }
                auto wpids = get<std::vector<int>>(jgroup, "wpids", {});
                for (int wpid : wpids) {
                    auto ichans = get_ordered_channels(anode, wpid);

                    log->debug("rule {}, group {}, wpid {}, chids:{}->{}, wire0head:{}->{}",
                               rule_num, group_num, WirePlaneId(std::abs(wpid)),
                               ichans.front()->ident(), ichans.back()->ident(),
                               ichans.front()->wires()[0]->ray().second,
                               ichans.back()->wires()[0]->ray().second);

                    for (auto ich : ichans) {
                        group.chid2row[ich->ident()] = group.chid2row.size();
                    }
                }
                if (group.chid2row.empty()) {
                    log->warn("rule {} group {} is empty for tag \"{}\", relpath {}",
                              rule_num, group_num, tag, relpath);
                    continue;
                }
                log->warn("rule {} group {} has {} channels for tag \"{}\", relpath {}",
                          rule_num, group_num, group.chid2row.size(), tag, relpath);

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
            next_count();
            return true;
        }

        TorchSemaphore sem(context());

        auto frame = frame_itensor(inframe);
        const std::string parent = frame->metadata()["datapath"].asString();

        
        ITorchTensor::vector itensors;
        itensors.push_back(frame);

        ITorchTensor::vector rules = rules_tensors(inframe, parent);
        if (rules.size()) {
            itensors.insert(itensors.end(), rules.begin(), rules.end());
        }

        ITorchTensor::vector chmasks = chmask_tensors(inframe, parent);
        if (chmasks.size()) {
            itensors.insert(itensors.end(), chmasks.begin(), chmasks.end());
        }

        Configuration empty;
        outtens = std::make_shared<SimpleTorchTensorSet>(inframe->ident(), empty,
                                                         to_device(itensors, device()));

        logit(outtens, "output");
        log->debug("my device is: {}", to_string(device()));

        next_count();
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

        const int frame_ident = iframe->ident();


        // std::unordered_map<int, std::vector<ITrace::pointer>> channel_traces;
        // std::unordered_map<int, size_t> chid_index;

        /// For each rules we need:
        /// 1. set of trace indices specific to tag and wpids/channels config.
        /// 2. map from channel ID to tensor row by the global sort.
        /// 3. iterate over indices, use map to find row, store.

        /// Keep track of each trace temporarily below
        struct TIR {
            ITrace::pointer trace;
            int64_t index; // into a summary or channels vector
            int64_t row;   // tensor row index
        };

        int rule_index = -1;
        for (const auto& rule : m_rules) {
            ++rule_index;
            std::string tag = rule.tag;
            
            auto tagged_traces = Aux::tagged_traces(iframe, tag);

            int group_index = -1;
            for (const auto& group : rule.groups) {
                ++group_index;
                
                std::vector<TIR> group_tirs;

                int64_t summary_index = -1;
                int tbeg=0, tend=0;
                // Go through traces, find which are compatible with group,
                // record their tensor row and summary index and sus out tensor
                // bounds
                for (auto trace : tagged_traces) {
                    ++summary_index;

                    auto rit = group.chid2row.find(trace->channel());
                    if (rit == group.chid2row.end()) {
                        continue;
                    }
                    const int64_t row = rit->second;

                    
                    int tbin = trace->tbin();
                    int nbins = trace->charge().size();
                    
                    group_tirs.emplace_back(TIR{trace, summary_index, row});

                    if (group_tirs.empty()) { // first time
                        tbeg = tbin;
                        tend = tbin+nbins;
                        continue;
                    }
                    if (tbeg < tbin) { tbeg = tbin; }
                    if (tend < tbin+nbins) { tend = tbin+nbins;}

                }
                int nrows = group.chid2row.size(); // span all channels as traces may be sparse
                int ncols = tend-tbeg;

                log->debug("ntirs={} ntraces={} nrows={} ncols={} nchid2row={}",
                           group_tirs.size(), tagged_traces.size(), nrows, ncols, group.chid2row.size());

  
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
                    for (const auto& tir : group_tirs) {
                        const auto& charge = tir.trace->charge();
                        int tbin = tir.trace->tbin();
                        int col_beg = tbin - tbeg;
                        int col_end = col_beg + charge.size();
                        // log->debug("rule={} grp={} row={} tbin={} beg={} end={}",
                        //            rule_index, group_index, tir.row, tbin, col_beg, col_end);
                        auto tmp = torch::tensor(charge, torch::kFloat32);
                        ten.index({tir.row, torch::indexing::Slice(col_beg, col_end)}) += tmp;
                    }

                    // traces have extra md
                    Configuration md = common_md;
                    md["tbin"] = (int)tbeg;
                    md["period"] = iframe->tick();
                    md["time"] = iframe->time();
                    tensors.push_back(make_datatype("traces", group.relpath, ten, md));
                }

                // Build summaries tensor
                auto summary = iframe->trace_summary(tag);
                if (summary.size()) {
                    torch::Tensor ten = torch::zeros({nrows}, torch::kFloat64);
                    for (const auto& tir : group_tirs) {
                        const double value = summary[tir.index];
                        ten.index({tir.row}) += value;
                    }
                    tensors.push_back(make_datatype("summaries", group.relpath, ten, common_md));
                }

                // Build chids tensor
                {
                    torch::Tensor ten = torch::zeros({nrows}, torch::kInt);
                    for (const auto& tir : group_tirs) {
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
        
        

