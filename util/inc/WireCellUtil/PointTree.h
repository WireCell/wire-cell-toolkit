/**
   A point cloud with points organized in an n-ary tree structure.

   See the Point class.
 */

#ifndef WIRECELLUTIL_POINTTREE
#define WIRECELLUTIL_POINTTREE

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointCloudCoordinates.h"
#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/NaryTreeFacade.h"
#include "WireCellUtil/KDTree.h"

#include "WireCellUtil/Logging.h" // debug

#include <boost/range/adaptors.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <ostream>

namespace WireCell::PointCloud::Tree {

    /** A point cloud scope describes selection of point cloud data
     * formed from a subset of tree nodes reached by a depth-first
     * descent.  
     */
    struct Scope {

        // The name of the node-local point clouds.
        std::string pcname{""};

        // The list of PC attribute array names to interpret as coordinates.
        using name_list_t = std::vector<std::string>;
        name_list_t coords{};

        // The depth of the descent.
        size_t depth{0};

        std::size_t hash() const;
        bool operator==(const Scope& other) const;
        bool operator!=(const Scope& other) const;

    };
    std::ostream& operator<<(std::ostream& o, Scope const& s);


}

namespace std {
    template<>
    struct hash<WireCell::PointCloud::Tree::Scope> {
        std::size_t operator()(const WireCell::PointCloud::Tree::Scope& scope) const {
            return scope.hash();
        }
    };
}

namespace WireCell::PointCloud::Tree {

    // An atomic, contiguous point cloud.
    using pointcloud_t = Dataset;

    // A set of point clouds each identified by a name.
    using named_pointclouds_t = std::map<std::string, pointcloud_t>;

    // We refer to point clouds held elsewhere.
    using pointcloud_ref = std::reference_wrapper<Dataset>;

    // Collect scoped things, see below
    class ScopedBase;
    template<typename ElementType> class ScopedView;

    class Points;
    using PointsNode = NaryTree::Node<Points>;

    /** Points is a payload value type for a NaryTree::Node.

        A Points instance stores a set of point clouds local to the node.
        Individual point clouds in the set are accessed by a name of type
        string.

        A Points also provides access to "scoped" objects.

        Points is also a Faced and so can accept a NaryTree::Facade.
     */
    class Points : public NaryTree::Faced<Points>
    {
        
      public:

        using self_t = Points;
        using base_t = NaryTree::Notified<Points>;

        using node_t = NaryTree::Node<Points>;
        using node_ptr = std::unique_ptr<node_t>;
        using node_path_t = std::vector<node_t*>;

        // template<typename ElementType>
        // using kdtree_t = typename KDTree<ElementType>::kdtree_type;

        Points() = default;
        virtual ~Points();

        // Copy constructor disabled due to holding unique k-d tree 
        Points(const Points& other) = delete;
        /// Move constructor.
        Points(Points&& other) = default;

        /// Copy assignment is deleted.
        Points& operator=(const Points& other) = delete;
        /// Move assignment
        Points& operator=(Points&& other) = default;

        /// Construct with local point clouds by copy
        explicit Points(const named_pointclouds_t& pcs, facade_ptr fac = nullptr)
            : m_lpcs(pcs.begin(), pcs.end()) {
            if (fac) set_facade(std::move(fac));
        }

        /// Construct with local point clouds by move
        explicit Points(named_pointclouds_t&& pcs, facade_ptr fac = nullptr)
            : m_lpcs(std::move(pcs)) {
            if (fac) set_facade(std::move(fac));
        }

        /// Access the set of point clouds local to this node.
        named_pointclouds_t& local_pcs() { return m_lpcs; }
        const named_pointclouds_t& local_pcs() const { return m_lpcs; }

        /// Access a scoped view.  The returned Scoped instance WILL be mutated
        /// when any new node is inserted into the scope.  The Scoped WILL be
        /// invalidated if any existing node is removed from the scope.
        template<typename ElementType=double>
        const ScopedView<ElementType>& scoped_view(const Scope& scope) const;
        template<typename ElementType=double>
        ScopedView<ElementType>& scoped_view(const Scope& scope);

        // Receive notification from n-ary tree to update existing
        // NFKDs if node is in any existing scope.
        virtual bool on_insert(const std::vector<node_type*>& path);

        // This is a brutal response to a removed node.  Any scope
        // containing the removed node will be removed from the cached
        // scoped data sets and k-d trees.  This will invalidate any
        // references and iterators from these objects that the caller
        // may be holding.
        virtual bool on_remove(const std::vector<node_type*>& path);

        // Return scoped view without creation, nullptr returned if the scope is
        // not yet created.
        const ScopedBase* get_scoped(const Scope& scope) const;
        ScopedBase* get_scoped(const Scope& scope);

      private:

        // our node-local point clouds
        named_pointclouds_t m_lpcs;

        // mutable cache
        using unique_scoped_t = std::unique_ptr<ScopedBase>;
        mutable std::unordered_map<Scope, unique_scoped_t> m_scoped;

        void init(const Scope& scope) const;
        
    };                          // Points

    template<typename ElementType>
    const ScopedView<ElementType>& Points::scoped_view(const Scope& scope) const
    {
        return const_cast<const ScopedView<ElementType>&>(
            const_cast<self_t*>(this)->scoped_view(scope));
    }
    template<typename ElementType>
    ScopedView<ElementType>& Points::scoped_view(const Scope& scope) 
    {
        using SV = ScopedView<ElementType>;
        auto * sbptr = get_scoped(scope);
        if (sbptr) {
            auto * svptr = dynamic_cast<SV*>(sbptr);
            if (svptr) {
                return *svptr;
            }
        }
        auto uptr = std::make_unique<SV>(scope);
        auto& sv = *uptr;
        m_scoped[scope] = std::move(uptr);
        init(scope);
        return sv;
    }

    // A scoped view on a subset of nodes in a NaryTree::Node<Points>.
    //
    // See also ScopedView<ElementType>.
    class ScopedBase {
      public:

        // The disjoint point clouds in scope
        using pointclouds_t = std::vector<pointcloud_ref>;

        // The nodes providing those point clouds
        using node_t = typename Points::node_t;
        using nodes_t = std::vector<node_t*>;

        // The arrays from each local PC that are in scope
        using selection_t = Dataset::selection_t;

        // This is a bit awkward but important.  We can not store the
        // selection_t by value AND assure stable memory locations AND
        // allow the store to grow AND keep the k-d tree updated as it
        // grows.  So, we store by pointer.
        using unique_selection_t = std::unique_ptr<selection_t>;

        // The selected arrays over all our PCs
        using selections_t = std::vector<unique_selection_t>;

        explicit ScopedBase(const Scope& scope) : m_scope(scope) {}
        virtual ~ScopedBase();

        const Scope& scope() const { return m_scope; }

        // Access the scoped nodes.
        const nodes_t& nodes() const { return m_nodes; }
        nodes_t& nodes() { return m_nodes; }

        // Access the scoped point cloud
        const pointclouds_t& pcs() const;

        // Total number of points across the scoped point cloud
        size_t npoints() const;

        // Access scoped point cloud as colleciton of selections.
        const selections_t& selections() const;

        // Resolve a node holding the point at the given index (eg as may be
        // returned from a k-d tree query).
        virtual       node_t* node_with_point(size_t point_index) = 0;
        virtual const node_t* node_with_point(size_t point_index) const = 0;

      protected:
        friend class Points;
        // Add a node that has been added to the tree in our scope.
        // This will also append to this scope's selections.
        virtual void append(node_t* node);

      private:

        void fill_cache();
        void fill_cache() const;

        Scope m_scope;
        nodes_t m_nodes;        // eager
        // The rest are filled lazily.
        // Mismatch with m_nodes.size() trigger cache fill.
        size_t m_node_count{0};
        size_t m_npoints{0};
        pointclouds_t m_pcs;
        selections_t m_selections;

    };


    // A scoped view with additional information that requires a specific
    // element type.
    template<typename ElementType=double>
    class ScopedView : public ScopedBase {
      public:
        
        //using nfkd_t = NFKDVec::Tree<ElementType>; // dynamic
        using nfkd_t = NFKDVec::Tree<ElementType, NFKDVec::IndexStatic>; // static

        explicit ScopedView(const Scope& scope)
            : ScopedBase(scope)
            , m_nfkd(nullptr)   // this is lazily created.
        {
        }
        virtual ~ScopedView() {}

        // Access the kdtree.
        //
        // The underlying k-d tree interface is lazy-loaded.  If force=true, any
        // previously created k-d tree is purged and a new one is made from the
        // existing selections.
        const nfkd_t& kd(bool force=false) const {
            return const_cast<ScopedView<ElementType>*>(this)->kd();
        }
        const nfkd_t& kd(bool force=false) {
            if (force) {
                m_nfkd = nullptr;
            }
            if (m_nfkd) { return *m_nfkd; }

            const auto& s = scope();
            m_nfkd = std::make_unique<nfkd_t>(s.coords.size());

            for (auto& sel : selections()) {
                m_nfkd->append(*sel);
            }
            return *m_nfkd;
        }

        // Return the node with the given point.
        const node_t* node_with_point(size_t point_index) const {
            const size_t block_index = kd().major_index(point_index);
            return nodes().at(block_index);
        }
        node_t* node_with_point(size_t point_index) {
            // This is truly a const operation.
            return const_cast<node_t*>( const_cast<const ScopedView<ElementType>*>(this)->node_with_point(point_index));
        }

      private:

        // This is actually mutable do to lazy behavior but the mutability is
        // assured via const_cast's in the methods..
        std::unique_ptr<nfkd_t> m_nfkd{nullptr};
    };

}

#endif
