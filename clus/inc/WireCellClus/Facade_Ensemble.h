#ifndef WIRECELLCLUS_FACADE_ENSEMBLE
#define WIRECELLCLUS_FACADE_ENSEMBLE

#include "WireCellClus/Facade_Util.h"

#include "WireCellUtil/PointTree.h"


namespace WireCell::Clus::Facade {
    class Grouping;

    struct EnsembleCache {
        /* nothing for now */
    };

    /** Give a node "Ensemble" semantics.
     *
     * This node has an "ensemble" of "grouping" nodes.  It does not have a
     * strong meaning other than a "group of groupings".  For example, an
     * ensemble may collect "live" and "dead" groupings and later a "shadow"
     * grouping may be added (eg, by retiling).
     *
     * Each child grouping is made or added through the ensemble with an
     * associated "name" by which the grouping may later be retrieved.  O.w.,
     * users are free to query the children for a desired grouping.  A grouping
     * holds its own name via metadata "name" entry.
     *
     */
    class Ensemble : public NaryTree::FacadeParent<Grouping, points_t>,
                     public Mixin<Ensemble, EnsembleCache> {
    public:

        Ensemble() : Mixin<Ensemble, EnsembleCache>(*this, "ensemble_scalar") {}
        virtual ~Ensemble() {}

        bool has(const std::string& name) const;
        std::vector<std::string> names() const;

        /// Return all children with a given name.
        std::vector<Grouping*> with_name(const std::string& name);
        std::vector<const Grouping*> with_name(const std::string& name) const;

        /// Make a named grouping.
        Grouping& make_grouping(const std::string& name);

        // Add and take ownership of existing grouping node, return its facade.
        Grouping& add_grouping_node(const std::string& name, points_t::node_ptr&& node);
    };
}
#endif
