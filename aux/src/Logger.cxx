#include "WireCellAux/Logger.h"

#include <sstream>
#include <iomanip>

using namespace WireCell;

Aux::Logger::Logger(const std::string& type_name,
                    const std::string& group_name)
    : m_iname("")
    , m_tname(type_name)
    , m_gname(group_name)
    , log(Log::logger(group_name)) // this will likely get a logger with shared sinks
{
}

Aux::Logger::~Logger() {}

static
std::string centered(std::string s, size_t target, char space = ' ')
{
    for (size_t now = s.size(); now < target; ++now) {
        if (now%2) {
            s.insert(s.end(), space);
        }
        else {
            s.insert(s.begin(), space);
        }
    }
    return s;
}

void Aux::Logger::set_name(const std::string& name)
{
    m_iname = name;

    std::string lname = m_gname + "/" + m_tname + "/" + m_iname;

    log = Log::logger(lname, false); // uniqe sinks so we can set unique pattern

    std::stringstream ss;
    ss << "[%H:%M:%S.%03e] %L [" << centered(m_gname, 8) << "] <" << m_tname << ":" << m_iname << "> %v %@";

    log->set_pattern(ss.str());

    log->debug("log name: \"{}\"", lname);
}

std::string Aux::Logger::get_name() const
{
    return m_iname;
}
