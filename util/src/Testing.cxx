#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Persist.h"

#include <sstream>

using namespace WireCell;

void boost::assertion_failed_msg(char const* expr, char const* msg, char const* function, char const* file, long line)
{
    // Intentionally a bare spdlog::critical (the process default logger): an
    // assertion can fail before any WCT sink is configured, and spdlog's
    // built-in default (stdout) guarantees the message is seen regardless.
    spdlog::critical("{}:{}:{} {} {}", file, function, line, expr, msg or "");

    std::stringstream ss;
    ss << "WireCell::AssertionError: \"" << expr << "\" in " << function << " " << file << ":" << line;
    if (msg and msg[0]) {
        ss << "\n" << msg;
    }
    THROW(AssertionError() << errmsg{ss.str()});
}

void boost::assertion_failed(char const* expr, char const* function, char const* file, long line)
{
    boost::assertion_failed_msg(expr, "", function, file, line);
}

void Testing::log(const char* argv0, const char* level)
{
    std::string name = argv0;
    name += ".log";
    Log::add_stderr(true, level);
    Log::add_file(name, level);
    Log::set_level(level);
}

void Testing::load_plugins(std::vector<std::string> list)
{
    PluginManager& pm = PluginManager::instance();

    if (list.empty()) {
        list = {"WireCellAux", "WireCellGen", "WireCellSigProc", "WireCellPgraph", "WireCellImg", "WireCellSio", "WireCellApps"};
    }
    for (const auto& one : list) {
        pm.add(one);
    }
}

Configuration Testing::detectors(const std::string& filename)
{
    return Persist::load(filename);
}

