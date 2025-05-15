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
     * users are free to query the children for a desired grouping.
     *
     * A grouping holds its own name (via metadata "name" entry) and thus it is
     * possible for more than one child Grouping to have the same name.  It is
     * up the to user to avoid or utilize this feature.
     *
     */
    class Ensemble : public NaryTree::FacadeParent<Grouping, points_t>,
                     public Mixin<Ensemble, EnsembleCache> {
    public:

        Ensemble() : Mixin<Ensemble, EnsembleCache>(*this, "ensemble_scalar") {}
        virtual ~Ensemble() {}

        /// Return false if no child Groupings have the name, else true.
        bool has(const std::string& name) const;

        /// Return all child Grouping names in child-order.  Redundant names
        /// will appear multiple times.
        std::vector<std::string> names() const;

        std::set<std::string> unique_names() const;

        /// Return all children with a given name.  In the case of multiple
        /// Groupings of the same name, the returned vector is in as-seen,
        /// child-order.  When no grouping has the name, the vector is empty.
        std::vector<Grouping*> with_name(const std::string& name);
        std::vector<const Grouping*> with_name(const std::string& name) const;

        /// Make a named grouping child node and return its grouping facade.
        Grouping& make_grouping(const std::string& name);

        // Add and take ownership of existing grouping node, return its facade.
        Grouping& add_grouping_node(const std::string& name, points_t::node_ptr&& node);

        // Return the FIRST grouping found for each name.
        std::map<std::string, Grouping*> groupings_by_name();
    };
}
#endif
