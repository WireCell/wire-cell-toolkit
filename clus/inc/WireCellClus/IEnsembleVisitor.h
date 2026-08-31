/**

   Interface providing an ensemble visitor.

   Primarily this is the interface to a "clustering method component" but is
   given a more generic name as some non-clustering pattern-recognition related
   components will use it.

 */
#ifndef WIRECELLCLUS_IENSEMBLEVISITOR
#define WIRECELLCLUS_IENSEMBLEVISITOR


#include "WireCellUtil/IComponent.h"
#include "WireCellClus/Facade_Ensemble.h"

#include <set>

namespace WireCell::Clus {

    class IEnsembleVisitor : public IComponent<IEnsembleVisitor> {
    public:

        virtual ~IEnsembleVisitor() {};

        /// Mutate the ensemble 
        virtual void visit(Facade::Ensemble& ensemble) const = 0;
    };
}


#endif
