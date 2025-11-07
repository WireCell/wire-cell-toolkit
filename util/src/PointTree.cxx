#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Logging.h"
#include <boost/container_hash/hash.hpp>

#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

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
    boost::hash_combine(h, name);
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
    // This leaves open a possible bug where the user adds and removes the same
    // number of nodes.  In that case, the old cache will be kept.
    if (m_node_count == m_nodes.size()) return;

    invalidate();               // propagate to subclass

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


Tree::pointcloud_t Tree::ScopedBase::flat_coords() const
{
    Tree::pointcloud_t flat;

    const Scope& s = scope();
    for (auto* node : m_nodes) {
        const auto& lpc = node->value.local_pc(s.pcname);
        flat.append(lpc.subset(s.coords));
    }
    return flat;
}

Tree::pointcloud_t Tree::ScopedBase::flat_pc(const std::string& pcname,
                                             const Dataset::name_list_t& arrnames) const
{
    Tree::pointcloud_t flat;
    
    // arbitrary PC name
    for (auto* node : m_nodes) {
        const auto& lpc = node->value.local_pc(pcname);
        if (arrnames.empty()) { // user gets all arrays
            flat.append(lpc);
        }
        else {                  // user wants select arrays
            flat.append(lpc.subset(arrnames));
        }
    }
    return flat;
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


bool Tree::Points::has_pc(const std::string& name) const
{
    auto it = m_lpcs.find(name);
    return it != m_lpcs.end();
}

const Tree::pointcloud_t& Tree::Points::local_pc(const std::string& name,
                                                 const pointcloud_t& defpc) const
{
    auto it = m_lpcs.find(name);
    if (it == m_lpcs.end()) return defpc;
    return it->second;
}

Tree::pointcloud_t& Tree::Points::local_pc(const std::string& name)
{
    return m_lpcs[name];
}

const Tree::ScopedBase* Tree::Points::get_scoped(const Scope& scope) const
{
    auto it = m_scoped.find(scope);
    if (it == m_scoped.end()) {
        return nullptr;
    }
    auto& sci = m_scoped[scope];
    if (!sci.indices_valid) {
        rebuild_indices(scope);
        sci.indices_valid = true;
    }
    return it->second.scoped.get();
}

Tree::ScopedBase* Tree::Points::get_scoped(const Scope& scope) 
{
    return const_cast<Tree::ScopedBase*>(
        const_cast<const self_t*>(this)->get_scoped(scope));
}

void WireCell::PointCloud::Tree::Points::rebuild_indices(const WireCell::PointCloud::Tree::Scope& scope) const
{
    auto& sci = m_scoped[scope];
    auto* svptr = dynamic_cast<ScopedView<double>*>(sci.scoped.get());
    if (!svptr) {
        return;
    }

    // Clear the index mappings
    svptr->clear_index_mappings();

    // For each node in the scope, add its points to the mapping
    size_t global_index = 0;
    for (auto& node : m_node->depth(scope.depth)) {
        auto& value = node.value;
        auto it = value.m_lpcs.find(scope.pcname);
        if (it == value.m_lpcs.end()) {
            continue;
        }

        bool want_if_in_scope = sci.nodes_in_scope.find(&node) != sci.nodes_in_scope.end();
        if (want_if_in_scope) {
            // For each point in this node, add its global index to the mapping
            for (size_t i = 0; i < it->second.size_major(); ++i) {
                svptr->append(global_index + i);
            }
        }
        global_index += it->second.size_major();
    }
    sci.indices_valid = true;
}

// Called new scoped view is created.
void WireCell::PointCloud::Tree::Points::init(const WireCell::PointCloud::Tree::Scope& scope) const
{
    auto& sci = m_scoped[scope]; // may create
    auto* svptr = dynamic_cast<ScopedView<double>*>(sci.scoped.get());

    size_t global_index = 0;

    // Walk the tree in scope, adding in-scope nodes.
    for (auto& node : m_node->depth(scope.depth)) { // depth part of scope.
        auto& value = node.value;
        auto it = value.m_lpcs.find(scope.pcname); // PC name part of scope.
        if (it == value.m_lpcs.end()) {
            continue;           // it is okay if node lacks PC, but such a node can not be in a scope
        }
        // Check for coordintate arrays on first construction. 
        Dataset& pc = it->second;
        assure_arrays(pc.keys(), scope); // throws if user logic error detected
        
        bool want_if_in_scope = sci.selector(node);
        // Tell scoped view about its new node if selector wants it
        if (want_if_in_scope) {
            sci.scoped->append(&node);
            sci.nodes_in_scope.insert(&node);  // Store node for quick lookup

            // If we have the ScopedView with index mapping, add point indices
            if (svptr) {
                // For each point in this node, add its global index to the mapping
                for (size_t i = 0; i < pc.size_major(); ++i) {
                    svptr->append(global_index + i);
                }
            }
        }
        global_index += pc.size_major();
    }

    sci.indices_valid = true;
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
        //     debug("\tinserted pc: {}", lit.first);
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
    // {                           // debug
    //     debug("inserting with path length {} in node with pcs:", path.size());
    //     // The scope must name a node-local PC 

    //     for (auto lit : local_pcs()) {
    //         debug("\tparent pc: {}", lit.first);
    //     }
    // }

    // auto* node = path.back();
    auto* node = path.front();

    // Give node to any views for which the node is in scope.
    for (auto& [scope,sci] : m_scoped) {
        // invalid indices ...
        sci.indices_valid = false;
        if (! in_scope(scope, node, path.size())) {
            continue;
        }
        if (sci.selector(*node)) {
            sci.scoped->append(node);
            sci.nodes_in_scope.insert(node);  // Store node for quick lookup
        }
    }
    return true;
}


bool Tree::Points::on_remove(const std::vector<node_type*>& path)
{
    auto* leaf = path.front();
    size_t psize = path.size();
    std::vector<Tree::Scope> dead;
    for (auto& [scope, sci] : m_scoped) {
        // invalid indices ...
        sci.indices_valid = false;
        if (in_scope(scope, leaf, psize)) {
            dead.push_back(scope);
        }
    }
    for (auto scope : dead) {
        m_scoped.erase(scope);
    }
        
    return true;                // continue ascent
}

std::string Tree::Points::as_string(bool recur, int level) const
{
    std::stringstream ss;
    std::string tab(level, ' ');
    ss << level << "\t" << tab << "pcs:[";
    for (auto lit : local_pcs()) {
        ss << " " << lit.first;
    }
    ss << " ] scopes:[";
    for (auto& [scope, sci] : m_scoped) {
        ss << " " << scope;
    }
    ss << " ]";
    if (! recur) {
        return ss.str();
    }
    auto children = m_node->children();
    ss << " " << children.size() << " children:\n";
    for (const auto& child : children) {
        ss << child->value.as_string(recur, level+1);
    }
    return ss.str();
}
