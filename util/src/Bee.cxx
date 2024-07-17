// TODO: add support for cluster_id array.

#include "WireCellUtil/Bee.h"
#include "WireCellUtil/Spdlog.h" // for fmt
#include "WireCellUtil/Persist.h" // for dumps(json)
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"
#include <boost/container_hash/hash.hpp>

#include <numeric>

using namespace WireCell;
        
std::string Bee::Object::json() const
{
    return Persist::dumps(m_data, 0, 6);
}

size_t Bee::Object::hash() const
{
    std::hash<Configuration> chash;
    return chash(m_data);
}

Bee::Points::Points()
{
    detector("");
    algorithm("");
    rse(0,0,0);
    m_data["x"] = Json::arrayValue;
    m_data["y"] = Json::arrayValue;
    m_data["z"] = Json::arrayValue;
    m_data["q"] = Json::arrayValue;
    m_data["cluster_id"] = Json::arrayValue;
}


Bee::Points::Points(const std::string& geom,
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
    m_data["cluster_id"] = Json::arrayValue;
}

void Bee::Points::detector(const std::string& geom)
{
    m_data["geom"] = geom;
}

void Bee::Points::algorithm(const std::string& type)
{
    m_name = type;
    m_data["type"] = type;
}

std::string Bee::Points::algorithm() const
{
    return m_data["type"].asString();
}

void Bee::Points::rse(int run, int sub, int evt)
{
    if (run >= 0) 
        m_data["runNo"] = run;
    if (sub >= 0)
        m_data["subRunNo"] = sub;
    if (evt >= 0)
        m_data["eventNo"] = evt;
}

void Bee::Points::reset(int evt, int sub, int run)
{
    rse(run,sub,evt);
    m_data["x"] = Json::arrayValue;
    m_data["y"] = Json::arrayValue;
    m_data["z"] = Json::arrayValue;
    m_data["q"] = Json::arrayValue;
    m_data["cluster_id"] = Json::arrayValue;
}

std::vector<int> Bee::Points::rse() const
{
    return {
        m_data["runNo"].asInt(),
        m_data["subRunNo"].asInt(),
        m_data["eventNo"].asInt()
    };
}

void Bee::Points::append(const Point& p, double q, int clid)
{
    // Here we take the unusual pattern to store a value in explicit units.
    // Normally, WCT code should NOT do this.  But, we consider m_data to belong
    // to Bee's domain.
    m_data["x"].append(p.x()/units::cm);
    m_data["y"].append(p.y()/units::cm);
    m_data["z"].append(p.z()/units::cm);
    m_data["q"].append(q);
    m_data["cluster_id"].append(clid);
}

void Bee::Points::append(const Bee::Points& obj)
{
    const int num = obj.size();
    const std::vector<std::string> xyzq = {"x","y","z","q","cluser_id"};
    for (const auto& key : xyzq) {
        const auto& arr = obj.m_data[key];
        for (int ind=0; ind<num; ++ind) {
            m_data[key].append(arr[ind]);
        }
    }
}

size_t Bee::Points::size() const
{
    return m_data["q"].size();
}
bool Bee::Points::empty() const
{
    return m_data["q"].empty();
}

int Bee::Points::back_cluster_id() const
{
    int siz = size();
    if (!siz) { return -1; }
    return m_data["cluser_id"][siz-1].asInt();
}


////// Patches

// Underlying m_data is a triple-nested array.
// array of patches
// each patch is array of pairs
// each pair gives (x,y).
Bee::Patches::Patches(const std::string& name, double tolerance, size_t minpts)
    : Object(name, Json::arrayValue)
    , m_tolerance(tolerance)
    , m_minpts(minpts)
{

}
        
void Bee::Patches::append(double y, double z)
{
    m_y.push_back(y);
    m_z.push_back(z);
}

void Bee::Patches::clear()
{
    m_y.clear();
    m_z.clear();
    m_data = Json::arrayValue;
}

void Bee::Patches::flush()
{
    if (m_z.empty()) {
        return;
    }

    const size_t npts = m_z.size();

    std::vector<size_t> inds(npts);
    std::iota(inds.begin(), inds.end(), 0);

    if (m_tolerance > 0) {
        auto far_less = [&](size_t a, size_t b) {
            if (std::abs(m_y[a] - m_y[b]) > m_tolerance) return m_y[a] < m_y[b];
            if (std::abs(m_z[a] - m_z[b]) > m_tolerance) return m_z[a] < m_z[b];
            return false;
        };
        std::set<size_t, decltype(far_less)> uindset(inds.begin(), inds.end(), far_less);
        inds.clear();
        inds.insert(inds.end(), uindset.begin(), uindset.end());
    }

    if (inds.size() < m_minpts) {
        m_y.clear();
        m_z.clear();
        return;
    }

    const double cy = std::accumulate(inds.begin(), inds.end(), 0.0,
                                      [&](double val, size_t ind) {
                                          return val + m_y[ind]; }) / inds.size();
    const double cz = std::accumulate(inds.begin(), inds.end(), 0.0,
                                      [&](double val, size_t ind) {
                                          return val + m_z[ind]; }) / inds.size();

    std::vector<double> angs(npts);
    for (auto ind : inds) {
        const double dy = m_y[ind] - cy;
        const double dz = m_z[ind] - cz;
        constexpr double pi = 3.14159265358979323846; // fixme: use std::numbers::pi in C++20
        angs[ind] = std::atan2(dz, dy) + pi; // from [-pi,pi] to [0,2pi]
    }

    std::sort(inds.begin(), inds.end(),
              [&](size_t a, size_t b) -> bool {
                  return angs[a] < angs[b];
              });

    // export to Bee json
    Json::Value jpatch = Json::arrayValue;
    for (auto ind : inds) {
        Json::Value jpt = Json::arrayValue;
        jpt.append(m_y[ind] / units::cm);
        jpt.append(m_z[ind] / units::cm);
        jpatch.append(jpt);
    }
    m_data.append(jpatch);

    m_y.clear();
    m_z.clear();
}

size_t Bee::Patches::size() const { return m_data.size(); }
bool Bee::Patches::empty() const { return m_data.empty(); }

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

    const auto name = obj.name();
    if (m_seen.find(name) != m_seen.end()) {
        flush();
    }
    m_seen.insert(name);
}

void Bee::Sink::flush()
{
    ++m_index;
    m_seen.clear();
}

std::string Bee::Sink::store_path(const Object& obj) const
{
    return fmt::format("data/{0}/{0}-{1}.json\n", m_index, obj.name());
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
        const std::string fname = "name " + store_path(obj);
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
    m_index=0;
    m_seen.clear();
    m_out.flush();
    m_out.pop();
    m_out.clear();
}
