#ifndef WIRECELLCLUS_FACADE_ENSEMBLE
#define WIRECELLCLUS_FACADE_ENSEMBLE

#include "WireCellClus/Facade_Mixins.h"
#include "WireCellClus/Facade_Util.h"

#include "WireCellUtil/PointTree.h"
#include "WireCellIface/ITensor.h"

#include <map>
#include <string>


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
    class Ensemble : public NaryTree::FacadeParent<Grouping, points_t>
                   , public Mixins::Cached<Ensemble, EnsembleCache> {
    public:

        Ensemble() : Mixins::Cached<Ensemble, EnsembleCache>(*this, "ensemble_scalar") {}
        virtual ~Ensemble() {}

        /// The ident of the ITensorSet this ensemble was loaded from -- for
        /// SBND that is the event number.  An Ensemble is built fresh for each
        /// event, so this is the one place a visitor can learn which event it
        /// is looking at.  Set by MultiAlgBlobClustering; -1 when unset.
        ///
        /// Its reason for existing: a process that streams many events needs
        /// its per-event file writers (the tracking ROOT file, the PR display
        /// dump) to name their output after the event rather than have every
        /// event overwrite one fixed name.
        int ident() const { return m_ident; }
        void set_ident(int ident) { m_ident = ident; }

        /// The (run, subrun, event) MultiAlgBlobClustering resolved for this
        /// event.  A visitor that stamps RSE into its output (the tracking ROOT
        /// trees, the PR display metadata) otherwise has only its own
        /// configure()-time numbers, which are one constant triplet for the
        /// whole process -- right for one-event-per-process, wrong for every
        /// event but the first of a group.  Valid only when set_rse() has been
        /// called; rse_valid() says whether it has.
        bool rse_valid() const { return m_rse_valid; }
        int runNo() const { return m_runNo; }
        int subRunNo() const { return m_subRunNo; }
        int eventNo() const { return m_eventNo; }
        void set_rse(int run, int subrun, int event) {
            m_runNo = run; m_subRunNo = subrun; m_eventNo = event; m_rse_valid = true;
        }

        /// Auxiliary per-event tensors that ride along with the point-cloud
        /// tree without being part of it, keyed by the tensor's metadata
        /// "datatype".  MultiAlgBlobClustering fills this from its input set
        /// (config "aux_datatypes", default the truth tables "truth_nu" /
        /// "truth_pf" of larwirecell's wclsTruthInformationAttacher) and
        /// forwards the same tensors to its output.  IEnsembleVisitor::visit()
        /// receives only the Ensemble, so this is how a visitor (e.g.
        /// SbndPrMagnifyTrackingVisitor's truth trees) sees them.
        /// aux_tensor() returns nullptr when the event carries none.
        void set_aux_tensor(const std::string& datatype, ITensor::pointer ten) { m_aux[datatype] = ten; }
        ITensor::pointer aux_tensor(const std::string& datatype) const {
            auto it = m_aux.find(datatype);
            return it == m_aux.end() ? nullptr : it->second;
        }
        const std::map<std::string, ITensor::pointer>& aux_tensors() const { return m_aux; }

        /// The metadata of the input tensor set this ensemble was made from,
        /// as MultiAlgBlobClustering received it (e.g. larwirecell's
        /// wclsTruthInformationAttacher "runNo"/"subRunNo"/"eventNo" and
        /// "input_file", ai-helper issue 39).  Read-only for visitors; null
        /// Configuration when the producer published none.
        void set_input_metadata(const Configuration& md) { m_input_metadata = md; }
        const Configuration& input_metadata() const { return m_input_metadata; }

        /// The Bee event index under which MultiAlgBlobClustering writes THIS
        /// event into its (possibly shared) Bee sink, published before the
        /// pipeline visitors run.  A visitor that writes Bee objects into the
        /// same sink (TaggerBeeVisitor) must use this index rather than keep a
        /// counter of its own, so the two writers cannot drift.  -1 = unset.
        void set_bee_index(int index) { m_bee_index = index; }
        int bee_index() const { return m_bee_index; }

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

       private:
        int m_ident{-1};
        bool m_rse_valid{false};
        int m_runNo{0};
        int m_subRunNo{0};
        int m_eventNo{0};
        std::map<std::string, ITensor::pointer> m_aux;
        Configuration m_input_metadata;
        int m_bee_index{-1};
    };
}
#endif
