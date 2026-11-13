#ifndef WIRECELLAUX_SIMPLEPRIMARYVERTEXSET
#define WIRECELLAUX_SIMPLEPRIMARYVERTEXSET

#include "WireCellIface/IPrimaryVertexSet.h"

namespace WireCell::Aux {

    // A primary vertex set that simply holds all the data it presents.
    class SimplePrimaryVertexSet : public WireCell::IPrimaryVertexSet {
        int m_ident;
        WireCell::IPrimaryVertex::shared_vector m_vertices;
        WireCell::Configuration m_metadata;

       public:
        SimplePrimaryVertexSet(int ident, const WireCell::IPrimaryVertex::vector& vertices,
                               const WireCell::Configuration& metadata = WireCell::Configuration())
          : m_ident(ident)
          , m_vertices(
                std::make_shared<WireCell::IPrimaryVertex::vector>(vertices.begin(), vertices.end()))
          , m_metadata(metadata)
        {
        }
        virtual ~SimplePrimaryVertexSet();
        virtual int ident() const { return m_ident; }
        virtual WireCell::IPrimaryVertex::shared_vector vertices() const { return m_vertices; }
        virtual WireCell::Configuration metadata() const { return m_metadata; }
    };

}  // namespace WireCell::Aux

#endif
