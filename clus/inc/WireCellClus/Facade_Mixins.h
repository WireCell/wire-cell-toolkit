#ifndef WIRECELL_CLUS_FACADEMIXINS
#define WIRECELL_CLUS_FACADEMIXINS

#include "WireCellClus/Graphs.h"
#include "WireCellUtil/PointTree.h"

#include <string>
#include <map>

namespace WireCell::Clus::Facade::Mixins {

    struct DummyCache{};


    /// The Ensemble/Grouping/Cluster/Facade classes inherit from this to gain
    /// common methods.  The mixin itself needs to know its facade type and
    /// instance but specifically does not cover any operation that requires
    /// knowledge of parent and children.
    ///
    /// It provides helper functions to deal with local PCs and an optional
    /// caching mechanism.  See comments on cache() and fill_cache() and
    /// clear_cache().  Note, using the cache mechanism does not preclude facade
    /// doing DIY caching.
    template<typename SelfType, typename CacheType=DummyCache>
    class Cached {
        SelfType& self;
        std::string scalar_pc_name, ident_array_name;
        mutable std::unique_ptr<CacheType> m_cache;
    public:
        Cached(SelfType& self, const std::string& scalar_pc_name, const std::string& ident_array_name = "ident")
            : self(self)
            , scalar_pc_name(scalar_pc_name)
            , ident_array_name(ident_array_name) {
            
        }

    protected:
        /// Facade cache management has a few simple rules, some optional:
        ///
        /// Cache rule 1:
        ///
        /// The SelfType MAY call this to access the cache instance.  It is
        /// guaranteed to have fill_cache() called.  See below.
        CacheType& cache() const
        {
            if (! m_cache) {
                m_cache = std::make_unique<CacheType>();
                fill_cache(* const_cast<CacheType*>(m_cache.get()));
            }
            return *m_cache.get();
        }

        /// Cache rule 2:
        ///
        /// Optionally, the SelfType may override this method in order to "bulk
        /// fill" the cache instance.
        virtual void fill_cache(CacheType& cache) const {}
        

        /// Cache rule 3:
        ///
        /// Optionally, the SelfType MAY implement lazy, fine-grained caching.
        /// This can be done with code such as:
        ///
        /// auto& mydata = cache().mydata;
        /// if (mydata.empty()) { /* fill/set mydata */ }


        /// Cache rule 4:
        ///
        /// Optionally, but not recommended, the SelfType may provide this
        /// method in order to do something when the clear() method is called.
        /// The cache instance will be removed just after this call returns.
        ///
        /// This is not recommended because you should be putting all cached
        /// items in the cache.
        ///
        /// DO NOT EXPOSE THIS METHOD.
        virtual void clear_cache() const
        {
            m_cache = nullptr;
        }

    public:

        /// Clear my node of all children nodes and purge my local PCs.
        /// Invalidates any cache.
        void clear()
        {
            // node level:
            self.node()->remove_children();
            // value level:
            self.local_pcs().clear();
            // facade cache level:
            clear_cache();
        }

        // Get the map from name to PC for all local PCs.
        WireCell::PointCloud::Tree::named_pointclouds_t& local_pcs()
        {
            return self.value().local_pcs();
        }
        const WireCell::PointCloud::Tree::named_pointclouds_t& local_pcs() const
        {
            return self.value().local_pcs();
        }

        // Return an "identifying number" from the "scalar" PC of the node.  As
        // with all "ident" values in WCT, there is no meaning ascribed to the
        // actual value (by WCT).  It is meant to refer to some external
        // identity.  If the scalar PC or the ident array are not found, the
        // default is returned.
        //
        // This is a special case method that merely delegates to get_scalar().
        int ident(int def = -1) const
        {
            return get_scalar<int>(ident_array_name, def);
        }

        // Set an ident number, delegating to set_scalar().
        void set_ident(int id)
        {
            set_scalar<int>(ident_array_name, id);
        }

        template <typename T = int>
        T get_element(const std::string& pcname, const std::string& aname, size_t index, T def = 0) const {
            const auto& lpcs = local_pcs();
            auto it = lpcs.find(pcname);
            if (it == lpcs.end()) {
                return def;
            }
            const auto arr = it->second.get(aname);
            if (!arr) {
                return def;
            }
            return arr->template element<T>(index);
        }

        // Return a value from the scalar PC
        template <typename T = int>
        T get_scalar(const std::string& aname, T def = 0) const {
            return get_element(scalar_pc_name, aname, 0, def);
        }
        
        // Set a value on the scalar PC
        template <typename T = int>
        void set_scalar(const std::string& aname, T val = 0) {
            auto& lpcs = local_pcs();
            auto& cs = lpcs[scalar_pc_name]; // create if not existing
            auto arr = cs.get(aname);
            if (!arr) {
                cs.add(aname, PointCloud::Array({(T)val}));
                return;
            }
            arr->template element<T>(0) = (T)val;
        }

        /// A flag is a name that can be "set" on a facade.  It is simply an
        /// entry in the scalar PC.  Most imply, a flag is Boolean false if
        /// unset (not defined) or has value 0 and set if defined with non-zero
        /// value.  Non-boolean values are allowed.  The flag name has a prefix
        /// (default "flag_") to provide a namespace.
        void set_flag(const std::string& name, int value=1, const std::string& prefix="flag_") {
            set_scalar<int>(prefix + name, value);
        }

        /// Get the value of a flag.  If the flag is unset, return the
        /// default_value.  See set_flag().
        int get_flag(const std::string& name, int default_value=0, const std::string& prefix="flag_") const {
            return get_scalar<int>(prefix + name, default_value);
        }

        /// Get all set flag names with a given prefix.
        std::vector<std::string> flag_names(const std::string& prefix="flag_") const {
            std::vector<std::string> ret;
            const auto& spc = get_pc(scalar_pc_name);
            for (const auto& key : spc.keys()) {
                if (String::startswith(key, prefix)) {
                    ret.push_back(key.substr(0, prefix.size()));
                }
            }
            return ret;
        }

        // Any flag set on the other will be set on this.
        void flags_from(const SelfType& other, const std::string& prefix="flag_") {
            for (const auto& fname : other.flag_names(prefix)) {
                set_flag(fname, other.get_flag(fname, 0, prefix), prefix);
            }
        }


        bool has_pc(const std::string& pcname) const
        {
            static PointCloud::Dataset dummy;
            const auto& lpcs = local_pcs();
            auto it = lpcs.find(pcname);
            if (it == lpcs.end()) {
                return false;
            }
            return true;
        }

        // Const access to a local PC/Dataset.  If pcname is missing return
        // reference to an empty dataset.
        const PointCloud::Dataset& get_pc(std::string pcname) const
        {
            if (pcname.empty()) {
                pcname = scalar_pc_name;
            }

            static PointCloud::Dataset dummy;
            const auto& lpcs = local_pcs();
            auto it = lpcs.find(pcname);
            if (it == lpcs.end()) {
                return dummy;
            }
            return it->second;
        }
        // Mutable access to a local PC/Dataset.  If pcname is missing, a new
        // dataset of that name will be created.
        PointCloud::Dataset& get_pc(std::string pcname)
        {
            if (pcname.empty()) {
                pcname = scalar_pc_name;
            }

            static PointCloud::Dataset dummy;
            auto& lpcs = local_pcs();
            return lpcs[pcname];
        }

        // Return true if this cluster has a PC array and PC of given names and type.
        template<typename ElementType=int>
        bool has_pcarray(const std::string& aname, std::string pcname) const {
            if (pcname.empty()) {
                pcname = scalar_pc_name;
            }

            auto& lpc = local_pcs();
            auto lit = lpc.find(pcname);
            if (lit == lpc.end()) {
                return false;
            }

            auto arr = lit->second.get(aname);
            if (!arr) {
                return false;
            }
            return arr->template is_type<ElementType>();
        }

        // Return as a span an array named "aname" stored in clusters PC named
        // by pcname.  Returns default span if PC or array not found or there is
        // a type mismatch.  Note, span uses array data in place.
        template<typename ElementType=int>
        PointCloud::Array::span_t<ElementType>
        get_pcarray(const std::string& aname, std::string pcname) {
            if (pcname.empty()) {
                pcname = scalar_pc_name;
            }

            auto& lpc = local_pcs();
            auto lit = lpc.find(pcname);
            if (lit == lpc.end()) {
                return {};
            }

            auto arr = lit->second.get(aname);
            if (!arr) {
                return {};
            }
            return arr->template elements<ElementType>();
        }
        template<typename ElementType=int>
        const PointCloud::Array::span_t<ElementType>
        get_pcarray(const std::string& aname, std::string pcname) const {
            if (pcname.empty()) {
                pcname = scalar_pc_name;
            }

            auto& lpc = local_pcs();
            auto lit = lpc.find(pcname);
            if (lit == lpc.end()) {
                return {};
            }

            auto arr = lit->second.get(aname);
            if (!arr) {
                return {};
            }
            return arr->template elements<ElementType>();
        }

        // Store vector as an array named "aname" into this cluster's PC named "pcname".
        // Reminder, all arrays in a PC must have same major size.
        template<typename ElementType=int>
        void
        put_pcarray(const std::vector<ElementType>& vec,
                    const std::string& aname, std::string pcname) {
            if (pcname.empty()) {
                pcname = scalar_pc_name;
            }

            auto &lpc = local_pcs();
            auto& pc = lpc[pcname];

            PointCloud::Array::shape_t shape = {vec.size()};

            auto arr = pc.get(aname);
            if (arr) {
                arr->template assign(vec.data(), shape, false);
            }
            else {
                pc.add(aname, PointCloud::Array(vec, shape, false));
            }
        }

        std::string get_name() const {
            const PointCloud::Dataset& spc = get_pc(scalar_pc_name);
            const auto& md = spc.metadata();
            auto jname = md["name"];
            if (jname.isString()) {
                return jname.asString();
            }
            return "";
        }

        void set_name(const std::string& name) {
            PointCloud::Dataset& spc = get_pc(scalar_pc_name);
            auto& md = spc.metadata();
            md["name" ] = name;
        }

    };

    class Graphs {
    public:
        /**
           Graph support.

           The facade owns every graph produced by this support and the graph
           dies with the facade.  
         */
        using graph_type = WireCell::Clus::Graphs::Weighted::Graph;
        using graph_store_type = std::map<std::string, graph_type>;

        /** Return true if named graph exists. */
        bool has_graph(const std::string& name) const;

        /** Return known graphs.

            This is only available as const.  User may use it to test for
            existence of a graph or iterate.
         */
        const graph_store_type& graph_store() const { return m_graph_store; }

        /**
           Create a graph of the given name.

           Replaces graph if it exists.
         */
        graph_type& make_graph(const std::string& name, size_t nvertices=0);

        /** Transfer a graph to the facade.

           Replaces graph if it exists.
        */
        graph_type& give_graph(const std::string& name, graph_type&& gr);

        /** Return a graph by name.

            Creates empty graph if one does not exist
        */
        graph_type& get_graph(const std::string& name);

        /** Return a graph by name.

            Throw KeyError if named graph does not exist.
        */
        const graph_type& get_graph(const std::string& name) const;

        /** Transfer ownership of graph.

            If the graph exists it is returned by a moved value.  Else an empty
            graph is returned.
        */
        graph_type take_graph(const std::string& name);

    private:
        graph_store_type m_graph_store;

    };


}

#endif
