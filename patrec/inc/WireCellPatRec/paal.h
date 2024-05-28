// This collects the few #includes of PAAL that are actually used by
// WireCellPatRec so that we can protect against warnings->errors due to their
// inclusion of boost.

#ifndef WIRECELLPATREC_PAAL
#define WIRECELLPATREC_PAAL

// These provide pragmas to tell the compiler to ignore warnings in boost.
#include "WireCellUtil/MultiArray.h"
#include "WireCellUtil/Graph.h"

// The actual few headers PatRec needs.
#include "WireCellPatRec/paal/data_structures/metric/graph_metrics.hpp"
#include "WireCellPatRec/paal/data_structures/metric/euclidean_metric.hpp"
#include "WireCellPatRec/paal/utils/irange.hpp"

#include "WireCellPatRec/paal/iterative_rounding/steiner_tree/steiner_tree.hpp"

#endif
