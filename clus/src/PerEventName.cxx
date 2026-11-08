#include "WireCellClus/PerEventName.h"
#include "WireCellClus/Facade_Ensemble.h"
#include "WireCellUtil/String.h"

using namespace WireCell;

namespace {
    void subst(std::string& s, const std::string& tok, int val)
    {
        const std::string v = std::to_string(val);
        for (auto p = s.find(tok); p != std::string::npos; p = s.find(tok, p + v.size())) {
            s.replace(p, tok.size(), v);
        }
    }
}

bool Clus::per_event_templated(const std::string& tmpl)
{
    // Only the known tokens: a literal '{' in a name is not a template (PR 544 review, item 10).
    if (tmpl.find('%') != std::string::npos) return true;
    for (const char* tok : {"{run}", "{subrun}", "{event}", "{ident}"}) {
        if (tmpl.find(tok) != std::string::npos) return true;
    }
    return false;
}

std::string Clus::per_event_name(const std::string& tmpl, int ident, int run, int subrun, int event)
{
    std::string s = tmpl;
    subst(s, "{run}", run);
    subst(s, "{subrun}", subrun);
    subst(s, "{event}", event);
    subst(s, "{ident}", ident);
    if (s.find('%') == std::string::npos) return s;
    return String::format(s, ident);
}

std::string Clus::per_event_name(const std::string& tmpl, const Facade::Ensemble& ensemble,
                                 int run, int subrun, int event)
{
    if (ensemble.rse_valid()) {
        run = ensemble.runNo();
        subrun = ensemble.subRunNo();
        event = ensemble.eventNo();
    }
    return per_event_name(tmpl, ensemble.ident(), run, subrun, event);
}
