/** Visit a PointGraph, possibly mutating it */

#ifndef WIRECELL_IPOINTGRAPHVISITOR
#define WIRECELL_IPOINTGRAPHVISITOR

#include "WireCellUtil/PointGraph.h"
#include "WireCellUtil/IComponent.h"

namespace WireCell {
    class IPointGraphVisitor : virtual public IComponent<IPointGraphVisitor> {
      public:

        virtual ~IPointGraphVisitor() = default;

        virtual void visit_pointgraph(PointGraph& pg) = 0;

    };
}

#endif
