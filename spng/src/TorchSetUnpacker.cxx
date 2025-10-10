#include "WireCellSpng/TorchSetUnpacker.h"
#include "WireCellUtil/NamedFactory.h"
#include <regex>

using namespace WireCell::HanaJsonCPP;                         

WIRECELL_FACTORY(SPNGTorchSetUnpacker,
                 WireCell::SPNG::TorchSetUnpacker,
                 WireCell::INamed,
                 WireCell::ITorchSetUnpacker,
                 WireCell::IConfigurable)


namespace WireCell::SPNG {

    std::string TorchSetUnpackerSelection::str() const {
        if (index >= 0) { return std::to_string(index); }
        return datapath;
    }
        

    TorchSetUnpacker::TorchSetUnpacker()
        : Logger("TorchSetUnpacker", "spng") {
    }

    TorchSetUnpacker::TorchSetUnpacker(const TorchSetUnpackerConfig& cfg)
        : Logger("TorchSetUnpacker", "spng")
        , m_cfg(cfg) {
        configme();
    }

    TorchSetUnpacker::~TorchSetUnpacker() {}

    WireCell::Configuration TorchSetUnpacker::default_configuration() const
    {
        auto cfg = this->Logger::default_configuration();
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void TorchSetUnpacker::configure(const WireCell::Configuration& cfg)
    {
        this->Logger::configure(cfg);
        from_json(m_cfg, cfg);
        configme();
    }

    void TorchSetUnpacker::configme()
    {
        if (m_cfg.selections.empty()) {
            raise<ValueError>("TorchSetUnpacker 'selections' parameter is empty, check config?");
        }
        
        // Build vector of callable selector functions that return nullptr on error.
        for (const auto& sel : m_cfg.selections) {
            if (sel.index >= 0) {
                size_t index = sel.index;
                if (index < 0) {
                    raise<ValueError>("TorchSetUnpacker, illegal negative tensor index, check config?");
                }
                m_selectors.emplace_back([index](ITorchTensorSet::pointer ts) -> ITorchTensor::pointer {
                    auto tvec = ts->tensors();
                    if (!tvec) return nullptr;
                    if (index < tvec->size()) {
                        return tvec->at(index);
                    }
                    return nullptr;
                });
                continue;
            }
            if (sel.datapath.size()) {
                std::regex match = std::regex(sel.datapath);
                m_selectors.emplace_back([match](ITorchTensorSet::pointer ts) -> ITorchTensor::pointer {
                    auto tvec = ts->tensors();
                    for (auto iten : *(tvec)) {
                        auto datapath = iten->metadata()["datapath"];
                        if (! datapath.isString()) { continue; }
                        if (std::regex_match(datapath.asString(), match)) {
                            return iten;
                        }
                    }
                    return nullptr;                        
                });
                continue;
            }
            raise<ValueError>("TorchSetUnpacker, invalid 'selections' entry, check config?");
        }
    }

    std::vector<std::string> TorchSetUnpacker::output_types()
    {
        const std::string tname = std::string(typeid(output_type).name());
        return std::vector<std::string>(multiplicity(), tname);
    }

    bool TorchSetUnpacker::operator()(const input_pointer& in, output_vector& outv)
    {
        const size_t nout = multiplicity();

        outv.clear();
        if (! in) {
            logit("EOS");
            outv.resize(nout, nullptr);
            ++m_count;
            return true;
        }

        for (size_t ind=0; ind<nout; ++ind) {
            const auto& sel = m_selectors[ind];
            auto iten = sel(in);
            if (!iten) {
                log->warn("null tensor selected for output port {} using selection: {}",
                          ind, m_cfg.selections[ind].str());
            }
            outv.push_back(iten);
        }
        return true;
    }
}
