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
    return Persist::dumps(asJson(), 0, 6);
}

size_t Bee::Object::hash() const
{
    std::hash<Configuration> chash;
    return chash(asJson());
}

Bee::Points::Points()
{
    detector("");
    algorithm("");
    rse(0,0,0);
}


Bee::Points::Points(const std::string& geom,
               const std::string& type,
               int run, int sub, int evt)
{
    detector(geom);
    algorithm(type);
    rse(run,sub,evt);
}

Configuration Bee::Points::asJson() const
{
    Configuration data = m_data;
    auto& jx = data["x"] = Json::arrayValue;
    auto& jy = data["y"] = Json::arrayValue;
    auto& jz = data["z"] = Json::arrayValue;
    auto& jq = data["q"] = Json::arrayValue;
    auto& jc = data["cluster_id"] = Json::arrayValue;
    auto& jr = data["real_cluster_id"] = Json::arrayValue;
    const size_t num = m_q.size();
    for (size_t ind = 0; ind < num; ++ind) {
        jx.append(m_x[ind]);
        jy.append(m_y[ind]);
        jz.append(m_z[ind]);
        jq.append(m_q[ind]);
        jc.append(m_clid[ind]);
        jr.append(m_real_clid[ind]);
    }
    return data;
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
    m_x.clear();
    m_y.clear();
    m_z.clear();
    m_q.clear();
    m_clid.clear();
    m_real_clid.clear();
}

std::vector<int> Bee::Points::rse() const
{
    return {
        m_data["runNo"].asInt(),
        m_data["subRunNo"].asInt(),
        m_data["eventNo"].asInt()
    };
}

void Bee::Points::append(const Point& p, double q, int clid, int real_clid)
{
    // Here we take the unusual pattern to store a value in explicit units.
    // Normally, WCT code should NOT do this.  But, we consider these columns
    // to belong to Bee's domain.
    m_x.push_back(p.x()/units::cm);
    m_y.push_back(p.y()/units::cm);
    m_z.push_back(p.z()/units::cm);
    m_q.push_back(q);
    m_clid.push_back(clid);
    m_real_clid.push_back(real_clid);
}

void Bee::Points::append(const Bee::Points& obj)
{
    m_x.insert(m_x.end(), obj.m_x.begin(), obj.m_x.end());
    m_y.insert(m_y.end(), obj.m_y.begin(), obj.m_y.end());
    m_z.insert(m_z.end(), obj.m_z.begin(), obj.m_z.end());
    m_q.insert(m_q.end(), obj.m_q.begin(), obj.m_q.end());
    m_clid.insert(m_clid.end(), obj.m_clid.begin(), obj.m_clid.end());
    m_real_clid.insert(m_real_clid.end(), obj.m_real_clid.begin(), obj.m_real_clid.end());
}

size_t Bee::Points::size() const
{
    return m_q.size();
}
bool Bee::Points::empty() const
{
    return m_q.empty();
}

int Bee::Points::back_cluster_id() const
{
    if (m_clid.empty()) { return -1; }
    return m_clid.back();
}


////// Patches

// Underlying m_data is a triple-nested array.
// array of patches
// each patch is array of pairs
// each pair gives (x,y).
//
// When m_tpc >= 0, m_data is instead a wrapper object
// {"version":2, "tpc":m_tpc, "polygons":<triple-nested array>}
// matching the wire-cell-bee3 v2 dead-area schema.
namespace {
    Json::Value init_patches_data(int tpc) {
        if (tpc < 0) return Json::Value(Json::arrayValue);
        Json::Value v(Json::objectValue);
        v["version"] = 2;
        v["tpc"] = tpc;
        v["polygons"] = Json::arrayValue;
        return v;
    }
}

Bee::Patches::Patches(const std::string& name, double tolerance, size_t minpts, int tpc)
    : Object(name, init_patches_data(tpc))
    , m_tolerance(tolerance)
    , m_minpts(minpts)
    , m_tpc(tpc)
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
    if (m_tpc < 0) {
        m_data = Json::arrayValue;
    } else {
        m_data["polygons"] = Json::arrayValue;
    }
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
    if (m_tpc < 0) {
        m_data.append(jpatch);
    } else {
        m_data["polygons"].append(jpatch);
    }

    m_y.clear();
    m_z.clear();
}

size_t Bee::Patches::size() const {
    return m_tpc < 0 ? m_data.size() : m_data["polygons"].size();
}
bool Bee::Patches::empty() const {
    return m_tpc < 0 ? m_data.empty() : m_data["polygons"].empty();
}

///// Flashes

namespace {
    void init_flashes_arrays(Configuration& d)
    {
        d["op_t"] = Json::arrayValue;
        d["op_peTotal"] = Json::arrayValue;
        d["op_pes"] = Json::arrayValue;
        d["op_pes_pred"] = Json::arrayValue;
        d["op_cluster_ids"] = Json::arrayValue;
        d["apa"] = Json::arrayValue;
    }
}

Bee::Flashes::Flashes()
    : Object("")
{
    detector("");
    rse(0, 0, 0);
    init_flashes_arrays(m_data);
}

Bee::Flashes::Flashes(const std::string& geom, const std::string& name,
                      int run, int sub, int evt)
    : Object(name)
{
    detector(geom);
    rse(run, sub, evt);
    init_flashes_arrays(m_data);
}

void Bee::Flashes::detector(const std::string& geom)
{
    m_data["geom"] = geom;
}

void Bee::Flashes::rse(int run, int sub, int evt)
{
    if (run >= 0) m_data["runNo"] = run;
    if (sub >= 0) m_data["subRunNo"] = sub;
    if (evt >= 0) m_data["eventNo"] = evt;
}

void Bee::Flashes::reset(int evt, int sub, int run)
{
    rse(run, sub, evt);
    init_flashes_arrays(m_data);
}

void Bee::Flashes::append(double t, const std::vector<double>& pes, double peTotal,
                          const std::vector<int>& cluster_ids,
                          const std::vector<double>& pes_pred,
                          int apa)
{
    m_data["op_t"].append(t);
    m_data["op_peTotal"].append(peTotal);

    Json::Value jpes(Json::arrayValue);
    for (auto v : pes) jpes.append(v);
    m_data["op_pes"].append(jpes);

    Json::Value jpred(Json::arrayValue);
    for (auto v : pes_pred) jpred.append(v);
    m_data["op_pes_pred"].append(jpred);

    Json::Value jcid(Json::arrayValue);
    for (auto c : cluster_ids) jcid.append(c);
    m_data["op_cluster_ids"].append(jcid);

    // Legacy dump_light stored the per-flash APA as a string ("0"/"1").
    m_data["apa"].append(std::to_string(apa));
}

void Bee::Flashes::set_groups(const std::vector<int>& groups)
{
    Json::Value jg(Json::arrayValue);
    for (int g : groups) jg.append(g);
    m_data["op_flash_group"] = jg;
}

void Bee::Flashes::set_t1(const std::vector<double>& t1)
{
    Json::Value jt(Json::arrayValue);
    for (double t : t1) jt.append(t);
    m_data["op_t1"] = jt;
}

size_t Bee::Flashes::size() const
{
    return m_data["op_t"].size();
}
bool Bee::Flashes::empty() const
{
    return m_data["op_t"].empty();
}


///// Sink

Bee::Sink::Sink()
{
}

Bee::Sink::Sink(const std::string& store)
{
    reset(store);
}

Bee::Sink::Sink(const std::string& store, size_t initial_index)
{
    reset(store, initial_index);
}

void Bee::Sink::reset(const std::string& store)
{
    reset(store, 0);  // Use default index 0
}

Bee::Sink::~Sink()
{
    close();
}

void Bee::Sink::reset(const std::string& store, size_t initial_index)
{
    close();
    Stream::output_filters(m_out, store);
    if (m_out.empty()) {
        raise<IOError>("no output for %s", store);
    }
    m_index = initial_index;
}

void Bee::Sink::set_index(size_t index)
{
    m_index = index;
    m_seen.clear();  // Clear seen items when changing index
}

size_t Bee::Sink::get_index() const
{
    return m_index;
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

        // Materialize the object's JSON DOM (for Points this builds the
        // point arrays from compact columns).
        Json::Value json_data = obj.asJson();

        std::string json;

        if (json_data.isObject()) {
            // Add RSE values
            json_data["runNo"] = std::to_string(m_runNo);
            json_data["subRunNo"] = std::to_string(m_subRunNo);
            json_data["eventNo"] = std::to_string(m_eventNo);
            // Convert to string
            //Json::FastWriter writer;
            json = Persist::dumps(json_data, 0, 6);
        } else {
            // dead area data is an array not an object
            // and does not have RSE values
            // ref: https://www.phy.bnl.gov/twister/bee/set/uboone/imaging/dead-region-2/event/0/channel-deadarea/
            json = Persist::dumps(json_data, 0, 6);
        }

        //const std::string json = obj.json();
        const std::string body = fmt::format("body {}\n", json.size());
        m_out << fname << body << json.data();
        m_out.flush();
    }

    return m_index;
}

void Bee::Sink::close()
{
    if (m_out.empty()) { return; }
    //m_index=0;
    m_seen.clear();
    m_out.flush();
    m_out.pop();
    m_out.clear();
}

// ---- Bee::ParticleTree ----
//
// The serialized form is a bare JSON array, exactly matching the prototype "mc"
// output read by the Bee viewer.  No object wrapper, no RSE metadata.

Bee::ParticleTree::ParticleTree(const std::string& name)
{
    m_name = name;
    m_data = Json::arrayValue;
}

void Bee::ParticleTree::reset()
{
    m_data = Json::arrayValue;
}

void Bee::ParticleTree::set_particles(const Configuration& particles_array)
{
    m_data = particles_array;
}

bool Bee::ParticleTree::empty() const
{
    return m_data.empty();
}

size_t Bee::ParticleTree::size() const
{
    return m_data.size();
}
