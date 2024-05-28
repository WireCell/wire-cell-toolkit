#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Logging.h"
#include <boost/container_hash/hash.hpp>

#include <algorithm>
#include <vector>
#include <string>

using spdlog::debug;
using namespace WireCell::PointCloud;

Tree::Points::~Points()
{
}


//
//  Scope
//

std::size_t Tree::Scope::hash() const
{
    std::size_t h = 0;
    boost::hash_combine(h, pcname);
    boost::hash_combine(h, depth);
    for (const auto& c : coords) {
        boost::hash_combine(h, c);
    }
    return h;
}

bool Tree::Scope::operator==(const Tree::Scope& other) const
{
    if (depth != other.depth) return false;
    if (pcname != other.pcname) return false;
    if (coords.size() != other.coords.size()) return false;
    for (size_t ind = 0; ind<coords.size(); ++ind) {
        if (coords[ind] != other.coords[ind]) return false;
    }
    return true;
}
bool Tree::Scope::operator!=(const Tree::Scope& other) const
{
    if (*this == other) return false;
    return true;
}

std::ostream& WireCell::PointCloud::Tree::operator<<(std::ostream& o, WireCell::PointCloud::Tree::Scope const& s)
{
    o << "<Scope \"" << s.pcname << "\" L" << s.depth;
    std::string comma = " ";
    for (const auto& cn : s.coords) {
        o << comma << cn;
        comma = ",";
    }
    o << ">";
    return o;
}

//
//  Scoped
//

Tree::ScopedBase::~ScopedBase ()
{
}

static void assure_arrays(const std::vector<std::string>& have, // ds keys
                          const Tree::Scope& scope)
{
    // check that it has the coordinate arrays
    std::vector<std::string> want(scope.coords), both, missing;
    std::sort(want.begin(), want.end());
    std::set_intersection(have.begin(), have.end(), want.begin(), want.end(),
                          std::back_inserter(both));
    if (both.size() == want.size()) {
        return;                 // arrays exist
    }

    // collect missing for exception message
    std::set_difference(have.begin(), have.end(), want.begin(), want.end(),
                        std::back_inserter(missing));
    std::string s;
    for (const auto& m : missing) {
        s += " " + m;
    }
    WireCell::raise<WireCell::IndexError>("Tree::Points data set missing arrays \"%s\" from scope %s",
                                          s, scope);    
}


void Tree::ScopedBase::append(ScopedBase::node_t* node)
{
    m_nodes.push_back(node);
}

void Tree::ScopedBase::fill_cache() const
{
    const_cast<ScopedBase*>(this)->fill_cache();
}

void Tree::ScopedBase::fill_cache()
{
    if (m_node_count == m_nodes.size()) return;

    m_node_count = 0;
    m_pcs.clear();
    m_npoints=0;
    m_selections.clear();

    const Scope& s = scope();
    for (auto* node : m_nodes) {
        Dataset& pc = node->value.local_pcs()[s.pcname];
        assure_arrays(pc.keys(), s); // sanity check
        m_pcs.push_back(std::ref(pc));
        m_npoints += pc.size_major();
        ++m_node_count;
        m_selections.emplace_back(std::make_unique<selection_t>(pc.selection(s.coords)));
    }
}

const Tree::ScopedBase::pointclouds_t& Tree::ScopedBase::pcs() const
{
    fill_cache();
    return m_pcs;
}


size_t Tree::ScopedBase::npoints() const
{
    fill_cache();
    return m_npoints;
}

const Tree::ScopedBase::selections_t& Tree::ScopedBase::selections() const
{
    fill_cache();
    return m_selections;
}

//
//  Points
//

// Tree::Scoped& Tree::Points::scoped(const Scope& scope)
// {
//     auto it = m_scoped.find(scope);
//     if (it == m_scoped.end()) {
//         raise<KeyError>("no scope %s", scope);
//     }
//     Scoped* sptr = it->second.get();
//     if (!sptr) {
//         raise<KeyError>("null scope %s", scope); // should not happen!
//     }
//     return *sptr;
// }

const Tree::ScopedBase* Tree::Points::get_scoped(const Scope& scope) const
{
    auto it = m_scoped.find(scope);
    if (it == m_scoped.end()) {
        return nullptr;
    }
    return it->second.get();
}

Tree::ScopedBase* Tree::Points::get_scoped(const Scope& scope) 
{
    return const_cast<Tree::ScopedBase*>(
        const_cast<const self_t*>(this)->get_scoped(scope));

    // auto it = m_scoped.find(scope);
    // if (it == m_scoped.end()) {
    //     return nullptr;
    // }
    // return it->second.get();
}



// Called new scoped view is created.
void WireCell::PointCloud::Tree::Points::init(const WireCell::PointCloud::Tree::Scope& scope) const
{
    auto& sv = m_scoped[scope];
    // Walk the tree in scope, adding in-scope nodes.
    for (auto& node : m_node->depth(scope.depth)) { // depth part of sceop.
        auto& value = node.value;
        auto it = value.m_lpcs.find(scope.pcname); // PC name part of scope.
        if (it == value.m_lpcs.end()) {
            continue;           // it is okay if node lacks PC
        }

        // Check for coordintate arrays on first construction. 
        Dataset& pc = it->second;
        assure_arrays(pc.keys(), scope); // throws if user logic error detected
        // Tell scoped view about its new node.
        sv->append(&node);
    }
}



// Return true one node at path_length relative to us matches the scope.
static
bool in_scope(const Tree::Scope& scope, const Tree::Points::node_t* node, size_t path_length)
{
    // path must be in our the scope depth
    // path length is 0 for me, 1 for children, 2 for ...
    // scope depth is 1 for me, 2 for children, 3 for ...

    if (scope.depth > 0 and path_length >= scope.depth) {
        // debug("not in scope: path length: {} scope:{}", path_length, scope);
        return false;
    }

    // The scope must name a node-local PC 
    auto& lpcs = node->value.local_pcs();

    // debug("in scope find pcname {}", scope.pcname);
    auto pcit = lpcs.find(scope.pcname);
    if (pcit == lpcs.end()) {
        // debug("not in scope: node has no named lpc in scope:{}", scope);
        // for (auto lit : lpcs) {
        //     debug("\tname: {}", lit.first);
        // }
        return false;
    }

    // The node-local PC must have all the coords.  coords may be
    // empty such as when the tested scope is a key in cached point
    // clouds (m_scoped_pcs).
    const auto& ds = pcit->second;
    for (const auto& name : scope.coords) {
        if (!ds.has(name)) {
            // debug("not in scope: lacks coord {} scope:{}", name, scope);
            return false;
        }
    }

    return true;    
}

bool Tree::Points::on_insert(const std::vector<node_type*>& path)
{
    auto* node = path.back();

    // Give node to any views for which the node is in scope.
    for (auto& [scope,sv] : m_scoped) {
        if (! in_scope(scope, node, path.size())) {
            continue;
        }
        sv->append(node);
    }
    return true;
}


bool Tree::Points::on_remove(const std::vector<node_type*>& path)
{
    auto* leaf = path.front();
    size_t psize = path.size();
    std::vector<Tree::Scope> dead;
    for (auto const& [scope, _] : m_scoped) {
        if (in_scope(scope, leaf, psize)) {
            dead.push_back(scope);
        }
    }
    for (auto scope : dead) {
        m_scoped.erase(scope);
    }
        
    return true;                // continue ascent
}


