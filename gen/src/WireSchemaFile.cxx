#include "WireCellGen/WireSchemaFile.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(WireSchemaFile, WireCell::Gen::WireSchemaFile, WireCell::IWireSchema, WireCell::IConfigurable)

using namespace WireCell;

Gen::WireSchemaFile::WireSchemaFile(const char* fname)
  : m_fname(fname)
{
}

Gen::WireSchemaFile::~WireSchemaFile() {}

void Gen::WireSchemaFile::configure(const WireCell::Configuration& cfg)
{
    m_fname = get(cfg, "filename", m_fname);
    if (m_fname.empty()) {
        THROW(ValueError() << errmsg{"must give a wire schema filename"});
    }
    int cor = (int)WireSchema::Correction::pitch;
    cor = get<int>(cfg, "correction", cor);
    m_store = WireSchema::load(m_fname.c_str(), (WireSchema::Correction) cor);
}

WireCell::Configuration Gen::WireSchemaFile::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_fname;
    cfg["correction"] = (int)WireSchema::Correction::pitch;
    return cfg;
}

const WireSchema::Store& Gen::WireSchemaFile::wire_schema_store() const { return m_store; }
