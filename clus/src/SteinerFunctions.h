/**
   These hold steiner-related free functions.

 */

#ifndef WIRECELLCLUS_STEINERFUNCTIONS
#define WIRECELLCLUS_STEINERFUNCTIONS

#include "SteinerGrapher.h"

namespace WireCell::Clus::Steiner {

    // Reserved stubs — ports of the prototype single-cluster Improve_PR3DCluster
    // overload.  That prototype path is unused in the NeutrinoID pipeline so
    // these are not yet implemented.  Do NOT call them.
    void improve_grapher_2(Grapher& grapher/*, ...*/);
    void improve_grapher_1(Grapher& grapher/*,...*/);
    void improve_grapher(Grapher& grapher/*,...*/);
    void improve_grapher(Grapher& grapher, Grapher& other_grapher/*,...*/);

}


#endif 
