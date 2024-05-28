#include "WireCellPgraph/Port.h"
#include "WireCellUtil/Type.h"

#include <sstream>

#include <iostream>

using namespace std;
using namespace WireCell::Pgraph;
using WireCell::demangle;

Port::Port(Node* node, Type type, std::string signature, std::string name)
  : m_node(node)
  , m_type(type)
  , m_name(name)
  , m_sig(signature)
  , m_edge(nullptr)
{
}

bool Port::isinput() const { return m_type == Port::input; }
bool Port::isoutput() const { return m_type == Port::output; }

Edge Port::edge() const { return m_edge; }

// Connect an edge, returning any previous one.
Edge Port::plug(Edge edge)
{
    Edge ret = m_edge;
    m_edge = edge;
    return ret;
}

// return edge queue size or 0 if no edge has been plugged
size_t Port::size() const
{
    if (!m_edge) {
        return 0;
    }
    return m_edge->size();
}

// Return true if queue is empty or no edge has been plugged.
bool Port::empty() const
{
    if (!m_edge or m_edge->empty()) {
        return true;
    }
    return false;
}

// Get the next data.  By default this pops the data off
// the queue.  To "peek" at the data, pas false.
Data Port::get(bool pop)
{
    if (isoutput()) {
        THROW(RuntimeError() << errmsg{"can not get from output port"});
    }
    if (!m_edge) {
        THROW(RuntimeError() << errmsg{"port has no edge"});
    }
    if (m_edge->empty()) {
        THROW(RuntimeError() << errmsg{"edge is empty"});
    }
    Data ret = m_edge->front();
    if (pop) {
        m_edge->pop_front();
    }
    return ret;
}

// Put the data onto the queue.
void Port::put(Data& data)
{
    if (isinput()) {
        THROW(RuntimeError() << errmsg{"can not put to input port"});
    }
    if (!m_edge) {
        THROW(RuntimeError() << errmsg{"port has no edge"});
    }
    m_edge->push_back(data);
}

const std::string& Port::name() const { return m_name; }
const std::string& Port::signature() const { return m_sig; }

std::string Port::ident() const
{
    std::stringstream ss;
    ss << "<Port"
       << " \"" << m_name << "\""
       << " sig:" << demangle(m_sig);
    if (isinput()) {
        ss << " (input)";
    }
    if (isoutput()) {
        ss << " (output)";
    }
    ss << ">";
    return ss.str();
}
