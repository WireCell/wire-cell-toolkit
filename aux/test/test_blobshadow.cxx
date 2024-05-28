#include "WireCellAux/BlobShadow.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellIface/ICluster.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellUtil/Graph.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::Aux;

struct FakeSlice : public ISlice {
    FakeSlice(int id) : _ident(id) {}
    virtual ~FakeSlice() {}
    virtual IFrame::pointer frame() const { THROW (RuntimeError()); }
    int _ident{-1};
    virtual int ident() const { return _ident; }
    virtual double start() const { THROW (RuntimeError()); }
    virtual double span() const { THROW (RuntimeError()); }
    virtual map_t activity() const { THROW (RuntimeError()); }
};

struct FakeBlob : public IBlob {
    FakeBlob(int id) : _ident(id) {}
    virtual ~FakeBlob() {}
    int _ident{-1};
    virtual int ident() const { return _ident; }
    virtual float value() const { THROW (RuntimeError()); }
    virtual float uncertainty() const { THROW (RuntimeError()); }
    virtual IAnodeFace::pointer face() const { THROW (RuntimeError()); }
    virtual ISlice::pointer slice() const { THROW (RuntimeError()); }
    virtual const RayGrid::Blob& shape() const { THROW (RuntimeError()); }
};

struct FakeWire : public IWire {
    FakeWire(int id) : _ident(id) {}
    virtual ~FakeWire() {}
    int _ident{-1};
    virtual int ident() const { return _ident; }
    virtual WirePlaneId planeid() const { return WirePlaneId(0); }
    virtual int index() const { return _ident; } // normally index != ident!
    virtual int channel() const { THROW (RuntimeError()); }
    virtual int segment() const { THROW (RuntimeError()); }
    virtual WireCell::Ray ray() const { THROW (RuntimeError()); }
};

struct FakeChannel : public IChannel {
    FakeChannel(int id) : _ident(id) {}
    virtual ~FakeChannel() {}
    int _ident{-1};
    virtual int ident() const { return _ident; }
    virtual int index() const { THROW (RuntimeError()); }
    virtual const IWire::vector& wires() const { THROW (RuntimeError()); }
};


template<typename IType, typename FType>
cluster_vertex_t add_node(int ident, cluster_graph_t& cg)
{
    auto ptr = new FType(ident);
    auto iptr = std::shared_ptr<const IType>(ptr);
    return boost::add_vertex(iptr, cg);
}

cluster_vertex_t add_node(char code, int ident, cluster_graph_t& cg)
{
    switch (code) {
        case 's':
            return add_node<ISlice, FakeSlice>(ident, cg);
        case 'b':
            return add_node<IBlob, FakeBlob>(ident, cg);
        case 'w':
            return add_node<IWire, FakeWire>(ident, cg);
        case 'c':
            return add_node<IChannel, FakeChannel>(ident, cg);
    }
    THROW (RuntimeError());
}

std::string dotify(const BlobShadow::graph_t& gr, const cluster_graph_t& cg)
{
    std::stringstream ss;
    // using vertex_t = boost::graph_traits<BlobShadow::graph_t>::vertex_descriptor;

    boost::write_graphviz(ss, gr,
                          [&](std::ostream& out, cluster_vertex_t v) {
                              const auto& dat = gr[v];
                              const auto& nod = cg[dat.desc];
                              out << "[label=\"" << nod.code() << nod.ident() << "\"]";
                          },
                          [&](std::ostream& out, cluster_edge_t e) {
                              const auto& p = gr[e];
                              int num = p.end-p.beg;
                              out << "[label=\"" << num << "\"]";
                          });
    return ss.str() + "\n";
}

int main(int argc, const char* argv[])
{
    cluster_graph_t cg;
    auto s0 = add_node('s', 0, cg);
    auto b0 = add_node('b', 0, cg);
    auto w0 = add_node('w', 0, cg);
    auto b1 = add_node('b', 1, cg);
    auto b2 = add_node('b', 2, cg);

    auto s1 = add_node('s', 1, cg);
    auto b3 = add_node('b', 3, cg);

    auto w1 = add_node('w', 1, cg);
    auto w2 = add_node('w', 2, cg);
    auto w3 = add_node('w', 3, cg);
    
    boost::add_edge(s0, b0, cg);
    boost::add_edge(s0, b1, cg);
    boost::add_edge(s0, b2, cg);

    boost::add_edge(b0, w0, cg);
    boost::add_edge(b0, w2, cg);
    boost::add_edge(b0, w3, cg);

    boost::add_edge(b1, w1, cg);
    boost::add_edge(b1, w2, cg);
    boost::add_edge(b1, w3, cg);
    boost::add_edge(b2, w3, cg);

    boost::add_edge(s1, b3, cg);
    boost::add_edge(b3, w0, cg);
    boost::add_edge(b3, w1, cg);

    std::string name = argv[0];
    {
        std::ofstream outf(name + "_orig.dot");
        outf << dotify(cg);
    }

    {
        auto dgraph = cluster::directed::type_directed(cg);
        std::ofstream outf(name + "_cdir.dot");
        outf << dotify(dgraph);
    }

    BlobShadow::graph_t bsg = BlobShadow::shadow(cg, 'w');

    {
        std::ofstream outf(name + "_blob.dot");
        outf << dotify(bsg, cg);
    }

    std::cerr
        << "for name in orig cdir blob ; do \n"
        << "  dot -Tpdf -o " << name << "_${name}.pdf " << name << "_${name}.dot\n"
        << "done\n";

    return 0;
}
