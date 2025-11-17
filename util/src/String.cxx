#include "WireCellUtil/String.h"


std::vector<std::string> WireCell::String::split(const std::string& in, const std::string& delim)
{
    std::vector<std::string> chunks;
    if (in.empty()) {
        return chunks;
    }
    boost::split(chunks, in, boost::is_any_of(delim), boost::token_compress_on);
    return chunks;
}

std::string WireCell::String::replace(std::string s,
                                      const std::string& have, const std::string& want)
{
    size_t pos = 0;
    while ((pos = s.find(have, pos)) != std::string::npos) {
        s.replace(pos, have.length(), want);
        pos += want.length(); // Advance past the newly inserted text
    }
    return s;
}


std::pair<std::string, std::string> WireCell::String::parse_pair(const std::string& in, const std::string& delim)
{
    std::vector<std::string> chunks = split(in, delim);

    std::string first = chunks[0];
    std::string second = "";
    if (chunks.size() > 1) {
        second = chunks[1];
    }
    return make_pair(first, second);
}

bool WireCell::String::endswith(const std::string& whole, const std::string& part) {
    return boost::algorithm::ends_with(whole, part);
}

bool WireCell::String::startswith(const std::string& whole, const std::string& part) {
    return boost::algorithm::starts_with(whole, part);
}

// more:
// http://www.boost.org/doc/libs/1_60_0/doc/html/string_algo/quickref.html#idm45555128584624
