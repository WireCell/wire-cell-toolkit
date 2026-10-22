#include "WireCellTruth/TruthTableSource.h"
#include "WireCellTruth/TruthTable.h"

#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include <boost/filesystem.hpp>

WIRECELL_FACTORY(TruthTableSource, WireCell::Truth::TruthTableSource,
                 WireCell::INamed,
                 WireCell::ITensorSetSource, WireCell::IConfigurable)

using namespace WireCell;

Truth::TruthTableSource::TruthTableSource()
    : Aux::Logger("TruthTableSource", "truth")
{
}
Truth::TruthTableSource::~TruthTableSource()
{
}

Configuration Truth::TruthTableSource::default_configuration() const
{
    Configuration cfg;
    cfg["inname"] = m_inname;
    cfg["prefix"] = m_prefix;
    cfg["idents"] = Json::arrayValue;
    for (int i : m_idents) cfg["idents"].append(i);
    cfg["stop_on_missing"] = m_stop_on_missing;
    return cfg;
}

void Truth::TruthTableSource::configure(const Configuration& cfg)
{
    m_inname = get(cfg, "inname", m_inname);
    m_prefix = get(cfg, "prefix", m_prefix);
    if (cfg.isMember("idents")) {
        m_idents.clear();
        for (const auto& j : cfg["idents"]) m_idents.push_back(j.asInt());
    }
    m_stop_on_missing = get(cfg, "stop_on_missing", m_stop_on_missing);
    if (m_inname.empty()) {
        THROW(ValueError() << errmsg{"TruthTableSource: inname is required"});
    }
}

bool Truth::TruthTableSource::operator()(ITensorSet::pointer& out)
{
    out = nullptr;
    if (m_eos) {
        return false;
    }
    if (m_count >= m_idents.size()) {
        m_eos = true;
        log->debug("EOS after {} tensor sets", m_count);
        return true;
    }
    const int ident = m_idents[m_count];
    const std::string fname = m_inname.find('%') == std::string::npos ? m_inname : String::format(m_inname, ident);
    if (!boost::filesystem::exists(fname)) {
        if (m_stop_on_missing) {
            m_eos = true;
            log->warn("no truth file {} for ident {}: EOS after {} tensor sets", fname, ident, m_count);
            return true;
        }
        THROW(IOError() << errmsg{"TruthTableSource: no truth file " + fname});
    }
    auto tt = TruthTable::load(fname, log);
    const std::string prefix = m_prefix.find('%') == std::string::npos ? m_prefix : String::format(m_prefix, ident);
    auto tv = std::make_shared<ITensor::vector>(tt.as_tensors(ident, prefix));
    Configuration setmd;
    setmd["ident"] = ident;
    setmd["truth_file"] = fname;
    setmd["truth_source"] = tt.source;
    if (!tt.meta.isNull()) setmd["truth_meta"] = tt.meta;
    out = std::make_shared<Aux::SimpleTensorSet>(ident, setmd, tv);
    log->debug("call={} ident {} from {}: {} nu rows, {} pf rows", m_count, ident, fname, tt.nu.size(), tt.pf.size());
    ++m_count;
    return true;
}
