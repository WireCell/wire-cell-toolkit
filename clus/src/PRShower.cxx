#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRGraph.h"
namespace WireCell::Clus::PR {

    Shower::Shower(Graph& graph)
        : TrajectoryView(graph)
    {
    }

    Shower::~Shower()
    {
    }


    VertexPtr Shower::start_vertex()
    {
        return m_start_vertex;
    }

    SegmentPtr Shower::start_segment()
    {
        return m_start_segment;
    }


    // Chainable setters

    /// Chainable setter of start vertex
    Shower& Shower::start_vertex(VertexPtr vtx)
    {
        if (! vtx->descriptor_valid()) {
            m_start_vertex = nullptr;
            return *this;
        }
        this->add_vertex(vtx);
        m_start_vertex = vtx;
        return *this;
    }
    
    
    /// Chainable setter of start segment
    Shower& Shower::start_segment(SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) {
            m_start_segment = nullptr;
            return *this;
        }
        this->add_segment(seg);
        m_start_segment = seg;
        return *this;
    }



}
