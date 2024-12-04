#include "WireCellUtil/NaryTesting.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell;
using namespace WireCell::NaryTesting;
using spdlog::debug;

// static member
size_t Introspective::live_count = 0;

NaryTesting::Introspective::Introspective()
    : name("(default-constructed)")
{
    ++ndef;
    ++live_count;
}

NaryTesting::Introspective::Introspective(const Introspective& d)
{
    name = d.name;
    ndef = d.ndef;
    nctor = d.nctor;
    ncopy = d.ncopy;
    nmove = d.nmove;
    nactions = d.nactions;
    ++ncopy;
    ++live_count;
}

NaryTesting::Introspective::Introspective(Introspective&& d)
{
    std::swap(name, d.name);
    std::swap(ndef, d.ndef);
    std::swap(nctor, d.nctor);
    std::swap(ncopy, d.ncopy);
    std::swap(nmove, d.nmove);
    std::swap(nactions, d.nactions);
    ++nmove;
    ++live_count;
}

NaryTesting::Introspective::Introspective(const std::string& s)
    : name(s)
{
    ++nctor;
    ++live_count;
}



NaryTesting::Introspective::~Introspective()
{
    --live_count;
}


void NaryTesting::Introspective::on_construct(node_type* n)
{
    node = n;
    debug("constructed {} with node@0x{}", name, (void*)node);
    ++nactions["constructed"];
}

bool NaryTesting::Introspective::on_insert(const std::vector<node_type*>& path)
{
    debug("insert {} {}", name, path.size());
    if (path.size() == 1) {
        ++nactions["inserted"];
    }
    return true;
}

bool NaryTesting::Introspective::on_remove(const std::vector<node_type*>& path)
{
    debug("remove {} {}", name, path.size());
    if (path.size() == 1) { 
        ++nactions["removing"];
    }
    return true;
}

bool NaryTesting::Introspective::on_ordered(const std::vector<node_type*>& path)
{
    debug("order {} {}", name, path.size());
    ++nactions["ordered"];
    return true;
}


std::ostream& NaryTesting::operator<<(std::ostream& o, const NaryTesting::Introspective& obj)
{
    o << "<Introspective " << obj.name << " #" << obj.nctor << ">";
    return o;
}

static
NaryTesting::Introspective::owned_ptr
make_a_node(const std::string name)
{
    return std::make_unique<Introspective::node_type>(Introspective(name));
}

//
//  Warning: tests rely on the structure of this tree.  Don't go changing it willy nilly. 
//
NaryTesting::Introspective::owned_ptr NaryTesting::make_simple_tree(const std::string& base_name)
{
    auto root = std::make_unique<Introspective::node_type>(Introspective(base_name));

    {
        Introspective d(base_name + ".0");
        root->insert(d);
    }

    {
        Introspective d(base_name + ".1");
        root->insert(std::move(d));
    }

    {
        root->insert(make_a_node(base_name + ".2"));
        // auto nup = std::make_unique<Introspective::node_type>(Introspective(base_name + ".2"));
        // root->insert(std::move(nup));
    }

    return root;
}

void NaryTesting::add_children(Introspective::node_type* node, std::list<size_t> layer_sizes)
{
    if (layer_sizes.empty()) return;

    std::string namedot = node->value.name + ".";
    const size_t num = layer_sizes.front();
    layer_sizes.pop_front();
    for (size_t ind = 0; ind < num; ++ind) {
        auto* child = node->insert(Introspective(namedot + std::to_string(ind)));
        add_children(child, layer_sizes);
    }
}


NaryTesting::Introspective::owned_ptr NaryTesting::make_layered_tree(std::list<size_t> layer_sizes)
{
    auto root = std::make_unique<Introspective::node_type>(Introspective("0"));
    add_children(root.get(), layer_sizes);
    return root;
}
