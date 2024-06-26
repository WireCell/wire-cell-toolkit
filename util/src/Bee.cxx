// TODO: add support for cluster_id array.

#include "WireCellUtil/Bee.h"
#include "WireCellUtil/Spdlog.h" // for fmt
#include "WireCellUtil/Persist.h" // for dumps(json)
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"
#include <boost/container_hash/hash.hpp>

using namespace WireCell;
        
Bee::Object::Object(const std::string& geom,
               const std::string& type,
               int run, int sub, int evt)
{
    detector(geom);
    algorithm(type);
    rse(run,sub,evt);
    m_data["x"] = Json::arrayValue;
    m_data["y"] = Json::arrayValue;
    m_data["z"] = Json::arrayValue;
    m_data["q"] = Json::arrayValue;
}

void Bee::Object::detector(const std::string& geom)
{
    m_data["geom"] = geom;
}

void Bee::Object::algorithm(const std::string& type)
{
    m_data["type"] = type;
}

std::string Bee::Object::algorithm() const
{
    return m_data["type"].asString();
}

void Bee::Object::rse(int run, int sub, int evt)
{
    m_data["runNo"] = run;
    m_data["subRunNo"] = sub;
    m_data["eventNo"] = evt;
}

void Bee::Object::reset(int evt)
{
    m_data["eventNo"] = evt;
    m_data["x"] = Json::arrayValue;
    m_data["y"] = Json::arrayValue;
    m_data["z"] = Json::arrayValue;
    m_data["q"] = Json::arrayValue;
}

std::vector<int> Bee::Object::rse() const
{
    return {
        m_data["runNo"].asInt(),
        m_data["subRunNo"].asInt(),
        m_data["eventNo"].asInt()
    };
}

std::string Bee::Object::json() const
{
    return Persist::dumps(m_data, 0, 6);
}

void Bee::Object::append(const Point& p, double q)
{
    // Here we take the unusual pattern to store a value in explicit units.
    // Normally, WCT code should NOT do this.  But, we consider m_data to belong
    // to Bee's domain.
    m_data["x"].append(p.x()/units::cm);
    m_data["y"].append(p.y()/units::cm);
    m_data["z"].append(p.z()/units::cm);
    m_data["q"].append(q);
}

void Bee::Object::append(const Bee::Object& obj)
{
    const int num = obj.size();
    const std::vector<std::string> xyzq = {"x","y","z","q"};
    for (const auto& key : xyzq) {
        for (int ind=0; ind<num; ++ind) {
            m_data[key].append(obj.m_data[key][ind]);
        }
    }
}

size_t Bee::Object::size() const
{
    return m_data["q"].size();
}
bool Bee::Object::empty() const
{
    return m_data["q"].empty();
}



///// Sink

Bee::Sink::Sink()
{
}

Bee::Sink::Sink(const std::string& store)
{
    reset(store);
}

Bee::Sink::~Sink()
{
    close();
}

void Bee::Sink::reset(const std::string& store)
{
    close();
    Stream::output_filters(m_out, store);
    if (m_out.empty()) {
        raise<IOError>("no output for %s", store);
    }
}

void Bee::Sink::index(const Object& obj)
{
    if (m_out.empty()) {
        raise<IOError>("Bee::Sink has no output stream");
    }

    // Get current object "fingerprints"
    size_t rse = 0;
    for (int num : obj.rse()) {
        boost::hash_combine(rse, num);
    }
    const auto alg = obj.algorithm();

    bool need_change = rse != m_rse || m_seen.find(alg) != m_seen.end();

    if (need_change) {
        if (m_rse) {
            // when m_rse is zero, it means first time in so use initial m_index==0.
            // when non-zero then it is time to advance.
            ++m_index;
        }
        m_rse = rse;
        m_seen.clear();
    }
    // Insert if we need change or not.
    m_seen.insert(alg);
}

size_t Bee::Sink::write(const Object& obj)
{
    index(obj);

    // Before, Bee relied on explicit Zip entries for any intermediate
    // directories to enumerate the Zip file tree structure.  Such directories
    // must be created explicitly.  We can tickle WCT iostream protocol and
    // miniz to do that by sending a "file" name ending in "/" and with an empty
    // body.  But, we turn that off after Bee was improved in order that this
    // hack does not lead to duplicate directories being created when this Bee
    // interface is used to produce files with multiple events and/or
    // algorithms.
    // {
    //     const std::string d1 = fmt::format("name data/\n", m_index);
    //     const std::string b1 = fmt::format("body 0\n");
    //     const std::string d2 = fmt::format("name data/{0}/\n", m_index);
    //     const std::string b2 = fmt::format("body 0\n");
    //     m_out << d1 << b1 << d2 << b2;
    // }

    // WCT stream protocol for actual file.
    {
        const std::string fname = fmt::format("name data/{0}/{0}-{1}.json\n", m_index, obj.algorithm());
        const std::string json = obj.json();
        const std::string body = fmt::format("body {}\n", json.size());

        m_out << fname << body << json.data();
        m_out.flush();
    }
    return m_index;
}

void Bee::Sink::close()
{
    if (m_out.empty()) { return; }
    m_out.flush();
    m_out.pop();
    m_out.clear();
}
