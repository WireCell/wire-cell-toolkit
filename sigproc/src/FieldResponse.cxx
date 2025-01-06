#include "WireCellSigProc/FieldResponse.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(FieldResponse, WireCell::SigProc::FieldResponse, WireCell::IFieldResponse, WireCell::IConfigurable)

using namespace WireCell;

SigProc::FieldResponse::FieldResponse(const char* frfilename)
  : m_fname(frfilename)
{
}

SigProc::FieldResponse::~FieldResponse() {}

WireCell::Configuration SigProc::FieldResponse::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_fname;
    return cfg;
}

void SigProc::FieldResponse::configure(const WireCell::Configuration& cfg)
{
    m_fname = get<std::string>(cfg, "filename", m_fname);
    if (m_fname.empty()) {
        THROW(ValueError() << errmsg{"must give field response filename"});
    }
    m_fr = Response::Schema::load(m_fname.c_str());
}

const Response::Schema::FieldResponse& SigProc::FieldResponse::field_response() const { return m_fr; }
