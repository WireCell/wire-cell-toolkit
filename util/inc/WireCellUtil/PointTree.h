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
#include "WireCellUtil/Limits.h" // SIZE_MAX

#include "WireCellUtil/Logging.h" // debug

#include <boost/range/adaptors.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <ostream>

namespace WireCell::PointCloud::Tree {

    /** A "scope" describes selection of point cloud data formed from a subset
        of tree nodes reached by a depth-first descent.
     
        A "scope" is typically used to created a "scoped view".  See
        Tree::scoped_view().
     */
    struct Scope {

        /// The name of the node-local point clouds that are in-scope.
        std::string pcname{""};

        /// The list of PC attribute array names to interpret as coordinates.
        using name_list_t = std::vector<std::string>;
        name_list_t coords{};

        /// The depth of the depth-first descent.  0 is unlimited.
        size_t depth{0};

        /**
           An additional name may be provided to distinguish otherwise identical
           scopes.  A unique name is typically required when a "filtered scope
           view" is constructed.  See Tree::scoped_view().
           them.
        */
        std::string name{""};

        /// Return a hash on the attributes.
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

        /// Access a mutable local PC by name.  If PC does not exist, create it.
        pointcloud_t& local_pc(const std::string& name);
        /// Access a const local PC by name.  If PC does not exist, return defpc.
        const pointcloud_t& local_pc(const std::string& name, const pointcloud_t& defpc = pointcloud_t()) const;

        bool has_pc(const std::string& name) const;

        /** Create a scoped view.

            Scoped views come in two types:

            - unfiltered :: the scope alone determines the included nodes.

            - filtered :: the scope and a selector function determines the included nodes.

            An unfiltered scoped view is equivalent to a filtered scope view
            with a selector function that returns only true.  Other selector
            functions may produce a scoped view that is a subset of the
            unfiltered scoped view produced from the same scope.

            When producing a scoped view, the user is strongly suggested to set
            the "scope.name" to a string that is unique among all the otherwise
            identical scopes.

            The ScopedViews are owned and cached by the PC tree using their
            scope as the identifying key.

            Existing scoped views are subject to tree mutations involving
            in-scope nodes being inserted or removed.

            - insertion :: if a newly inserted node is in a scope and if the
              selector function returns true for that new node then the new node
              is appended to the scope.  Note, this call of the selector
              function is isolated and not part of a DFS.

            - removal :: if an existing in-scope node is removed, the scoped
              view is INVALIDATED.  Any outstanding references to the scope will
              refer to garbage memory.

            See also get_scoped().
         */
        using SelectorFunction = std::function<bool(const Points::node_t& node)>;
        template<typename ElementType=double>
        const ScopedView<ElementType>& scoped_view(
            const Scope& scope,
            SelectorFunction selector = [](const Points::node_t&){return true;}) const;
        template<typename ElementType=double>
        ScopedView<ElementType>& scoped_view(
            const Scope& scope,
            SelectorFunction selector = [](const Points::node_t&){return true;});

        // Receive notification from n-ary tree to update existing
        // NFKDs if node is in any existing scope.
        virtual bool on_insert(const std::vector<node_type*>& path);

        // This is a brutal response to a removed node.  Any scope
        // containing the removed node will be removed from the cached
        // scoped data sets and k-d trees.  This will invalidate any
        // references and iterators from these objects that the caller
        // may be holding.
        virtual bool on_remove(const std::vector<node_type*>& path);

        // Return scoped view from cache without creation, nullptr returned if
        // the scope is not yet created or has been invalidated.
        const ScopedBase* get_scoped(const Scope& scope) const;
        ScopedBase* get_scoped(const Scope& scope);

        // Return representation of this as a string.  By default this recurs
        // into children.  The "level" is used to provide indentation.
        std::string as_string(bool recur=true, int level=0) const;

      private:

        // our node-local point clouds
        named_pointclouds_t m_lpcs;

        // mutable cache
    //  public: // debugging purpose ...
        using unique_scoped_t = std::unique_ptr<ScopedBase>;
        struct ScopedViewCacheItem {
            unique_scoped_t scoped;
            SelectorFunction selector;

            // cache of indices ...
            bool indices_valid{false};  // Flag to track if indices need rebuilding
            // nodes in scope ...
            std::unordered_set<node_t*> nodes_in_scope;
        };
        mutable std::unordered_map<Scope, ScopedViewCacheItem> m_scoped;

        void init(const Scope& scope) const;
        void rebuild_indices(const Scope& scope) const;
        
    };                          // Points

    template<typename ElementType>
    const ScopedView<ElementType>& Points::scoped_view(const Scope& scope, SelectorFunction selector) const
    {
        return const_cast<const ScopedView<ElementType>&>(
            const_cast<self_t*>(this)->scoped_view(scope, selector));
    }
    template<typename ElementType>
    ScopedView<ElementType>& Points::scoped_view(const Scope& scope, SelectorFunction selector) 
    {
        using SV = ScopedView<ElementType>;
        auto * sbptr = get_scoped(scope);
        if (sbptr) {
            auto * svptr = dynamic_cast<SV*>(sbptr);
            if (svptr) {
                auto& sci = m_scoped[scope];
                if (!sci.indices_valid) {
                    rebuild_indices(scope);
                    sci.indices_valid = true;
                }
                return *svptr;
            }
        }
        auto uptr = std::make_unique<SV>(scope);
        auto& sv = *uptr;
        auto& sci = m_scoped[scope];
        sci.scoped = std::move(uptr);
        sci.selector = selector;
        sci.indices_valid = false; // Start with invalid indices

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
        
        /// Return a contiguous, monolithic dataset of the scoped coordinate PC.
        ///
        /// This "flat" PC represents a COPY for the per-node local PC
        /// coordinate arrays.  In particular, no connection with the scoped
        /// view is retained and the returned dataset will not be mutated if the
        /// underlying PC tree is later mutated.
        pointcloud_t flat_coords() const;

        /// Return a contiguous, monolithic dataset built from the PCs of the
        /// given name local to the scoped nodes.  The arrays in the PC are
        /// limited to those named in the name list the list is empty in which
        /// case all arrays in the named local PCs are returned.
        ///
        /// This "flat" PC represents a COPY for the per-node local PC arrays.
        /// In particular, no connection with the scoped view is retained and
        /// the returned dataset will not be mutated if the underlying PC tree
        /// is later mutated.
        pointcloud_t flat_pc(const std::string& pcname,
                             const Dataset::name_list_t& arrnames = Dataset::name_list_t()) const;


        /// Return a contiguous, monolithic vector built from the arrays of the
        /// given "aname" from the PCs of given pcname on the nodes in the
        /// scope.
        ///
        /// If pcname is not given, the scope's PCs are used.
        ///
        /// If pcname matches scope's PC name and aname names a scope coord then
        /// the returned vector will be aligned with the in-scope points.  If
        /// either names are out-of-scope, there is no alignment guarantee.  If
        /// any in-scope node lacks the PC or the array, its contribution is
        /// simply omitted.
        template<typename T>
        std::vector<T> flat_vector(const std::string& aname, std::string pcname="") const {
            if (pcname.empty()) {
                pcname = m_scope.pcname;
            }

            std::vector<T> flat;
    
            for (auto* node : m_nodes) {
                if (! node->value.has_pc(pcname)) {
                    continue;
                }
                const auto& lpc = node->value.local_pc(pcname);
                auto aptr = lpc.get(aname);
                if (aptr == nullptr) {
                    continue;
                }
                auto ele = aptr->elements<T>();
                flat.insert(flat.end(), ele.begin(), ele.end());
            }
            return flat;
        }

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
        // Available for ScopedView<T> to implement to invalidate typed portion of the class.
        virtual void invalidate() {}

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

            // A little tricky order of operations here.  This invalidates the kd:
            const auto& sels = selections();

            // NOW we can make and use it without triggering invalidate().
            m_nfkd = std::make_unique<nfkd_t>(s.coords.size());
            assert(m_nfkd);
            
            for (auto& sel : sels) {
                assert(m_nfkd);
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

        // New methods for index mapping
        void clear_index_mappings() {
            m_global_to_local.clear();
            m_local_to_global.clear(); 
        }

        // Add a global index to the mapping - maintains both maps
        void append(size_t global_index) {
            size_t local_index = m_local_to_global.size();
            m_local_to_global.push_back(global_index);
            m_global_to_local[global_index] = local_index;
        }

        size_t global_to_local(size_t global_idx, size_t def=std::numeric_limits<size_t>::max()) const {
            auto it = m_global_to_local.find(global_idx);
            if (it != m_global_to_local.end()) {
                return it->second;
            }
            return def;
        }

        size_t local_to_global(size_t local_idx, size_t def=std::numeric_limits<size_t>::max()) const {
            if (local_idx < m_local_to_global.size()) {
                return m_local_to_global[local_idx];
            }
            return def;
        }

        virtual void invalidate() {
            m_nfkd = nullptr;
        }

      private:

        // This is actually mutable do to lazy behavior but the mutability is
        // assured via const_cast's in the methods..
        std::unique_ptr<nfkd_t> m_nfkd{nullptr};

        // New member variables for indices ...
        std::unordered_map<size_t, size_t> m_global_to_local;
        std::vector<size_t> m_local_to_global;
    };

}

template <> struct fmt::formatter<WireCell::PointCloud::Tree::Scope> : fmt::ostream_formatter {};

#endif
