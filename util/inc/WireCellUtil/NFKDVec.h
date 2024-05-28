// An interface to nanoflan using std::vector stores

#ifndef WIRECELLUTIL_NFKDVEC
#define WIRECELLUTIL_NFKDVEC

#include "WireCellUtil/nanoflann.hpp"
#include "WireCellUtil/DetectionIdiom.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/PointCloudDataset.h"

#include <vector>
#include <numeric>              // iota

namespace WireCell::NFKDVec {

    /// Provide nanoflan index adaptor traits to make class templates
    /// easier to supply.
    struct IndexTraits {    };
    struct IndexStatic : public IndexTraits {
        template <typename Distance, 
                  class DatasetAdaptor, int32_t DIM = -1,
                  typename IndexType = size_t>
        struct traits {
            using index_t = nanoflann::KDTreeSingleIndexAdaptor<Distance, DatasetAdaptor, DIM, IndexType>;
        };
    };
    struct IndexDynamic : public IndexTraits {
        template <typename Distance, 
                  class DatasetAdaptor, int32_t DIM = -1,
                  typename IndexType = size_t>
        struct traits {
            using index_t = nanoflann::KDTreeSingleIndexDynamicAdaptor<Distance, DatasetAdaptor, DIM, IndexType>;
        };
    };

    /// Likewise for distance.  Here, nanoflann provides them so we
    /// simply forward their names for the ones supported here.
    using DistanceTraits = nanoflann::Metric;
    // L1 - sum of absolute linear difference on each coordinate.  The
    // radius and distances in results sets are in units of [LENGTH].
    using DistanceL1 = nanoflann::metric_L1;
    // L2 - sum of squared difference on each coordinate.  NOTE: the
    // value for this metric is in length units SQUARED.  The "radius"
    // given to a radius() query and the distances between the query
    // point and the points in the result set (including that from
    // knn() query) are all in units of [LENGTH]^2.  They are not NOT
    // in units of [LENGTH].
    using DistanceL2 = nanoflann::metric_L2;
    // L2 but optimize for low-dimension L2, 2D or 3D.  This is
    // default.
    using DistanceL2Simple = nanoflann::metric_L2_Simple;    


    // Interface to a k-d tree with points accessed by an iterator
    // range.
    template<typename ElementType
             ,typename IndexTraits = IndexDynamic
             ,typename DistanceTraits = DistanceL2Simple
             >    
    class Tree {

      public:

        using self_type = Tree<ElementType, IndexTraits, DistanceTraits>;

        // scalar point element type
        using element_type = ElementType;
        // 1D elements of a single coordinate
        using coordinates_type = std::vector<element_type>;
        // 2D collection of coordinates
        using points_type = std::vector<coordinates_type>;
        // Use int, instead of size_t, to save a little memory.
        using block_number_type = int;
        // Track from where a point came from by a major/minor index.
        using block_indices_type = std::vector<block_number_type>;

        // nanoflann types
        using metric_type = typename DistanceTraits::template traits<element_type, self_type>::distance_t;
        using nfkdindex_type = typename IndexTraits::template traits<metric_type, self_type>::index_t;
        using distance_type = typename metric_type::DistanceType;

        /// A result of a k-d query.  The first part of a pair gives
        /// an index of the point in the flat collection.  The second
        /// gives the distance from the query point to that point.
        /// NOTICE: using L2 metrics means "distance" is squared. 
        using iterator = size_t;
        using result_item = std::pair<iterator, distance_type>;
        using results_type = std::vector<result_item>;

        // Must not use this for static, may use for dynamic.
        explicit Tree(size_t dimensionality)
            : m_points(dimensionality)
        {
        }

        // Must use this for static, may use it for dynamic.
        explicit Tree(const points_type& points)
            : m_points(points.size())
        {
            append(points);
        }

        ~Tree()
        {
        }

        // Return the number of dimensions of the K-D space
        const size_t ndim() const { return m_points.size(); }

        // Access the collection of points.
        const points_type& points() const { return m_points; }
        size_t npoints() const {
            if (m_points.empty()) { return 0; }
            return m_points[0].size();
        }

        // Build and return point at index across coordinate arrays.
        std::vector<element_type> point(size_t index) const {
            const size_t ndims = ndim();
            std::vector<element_type> ret(ndims);
            for (size_t dim=0; dim<ndims; ++dim) {
                ret[dim] = m_points.at(dim).at(index);
            }
            return ret;
        }
        // Build and return a 3D point.  This is only a valid call if ndim==3.
        Point point3d(size_t index) const {
            if (ndim() != 3) {
                raise<LogicError>("NFKD::Tree: point3d requires 3 dimension, have %d", ndim());
            }
            return Point(m_points.at(0).at(index),
                         m_points.at(1).at(index),
                         m_points.at(2).at(index));
        }

        // Access the vector of block number by point index.
        const block_indices_type& major_indices() const { return m_major_indices; }
        const block_indices_type& minor_indices() const { return m_minor_indices; }

        // Number of blocks that have been appended.  Note, this is not
        // necessarily the value of the last element of major_indices.  Empty
        // blocks can be appended and they are counted but have no entries in
        // the {major,minor_indices.
        size_t nblocks() const {
            return m_nblocks;
        }

        // Return the number of the block that provided the point at the given index
        block_number_type major_index(size_t point_index) const {
            return m_major_indices.at(point_index);
        }

        // Return the index in its block that provided the point at index
        block_number_type minor_index(size_t point_index) const {
            return m_minor_indices.at(point_index);
        }

        // Append one PointCloud::Dataset selection (vector of PC arrays)
        void append(const PointCloud::Dataset::selection_t& sel) {
            const size_t block_index = m_nblocks++;

            if (sel.empty()) {
                return;
            }
            const size_t adding = sel[0]->size_major();
            if (!adding) {
                return;
            }

            // Extend the coordinate points vectors
            const size_t ndims = ndim();
            for (size_t dim=0; dim<ndims; ++dim) {
                auto& vdim = m_points[dim];
                auto sdim = sel[dim]->elements<element_type>();
                vdim.insert(vdim.end(), sdim.begin(), sdim.end());
            }
            const size_t newsize = m_points[0].size();

            // Update the major/minor map
            const size_t oldsize = m_major_indices.size();
            m_major_indices.resize(newsize, block_index);
            m_minor_indices.resize(newsize);
            std::iota(m_minor_indices.begin() + oldsize, m_minor_indices.end(), 0);
            this->addn<nfkdindex_type>(oldsize, adding);
        }

        // Append with vector of vector of element.
        void append(const points_type& pts) {
            const size_t block_index = m_nblocks++;

            if (pts.empty()) {
                return;
            }
            const size_t adding = pts[0].size();
            if (!adding) {
                return;
            }

            // Extend the coordinate vectors
            const size_t ndims = ndim();
            for (size_t dim=0; dim<ndims; ++dim) {
                auto& vdim = m_points[dim];
                const auto& sdim = pts[dim];
                vdim.insert(vdim.end(), sdim.begin(), sdim.end());
            }
            const size_t newsize = m_points[0].size();

            // Update the major/minor map
            const size_t oldsize = m_major_indices.size();
            m_major_indices.resize(newsize, block_index);
            m_minor_indices.resize(newsize);
            std::iota(m_minor_indices.begin() + oldsize, m_minor_indices.end(), 0);
            this->addn<nfkdindex_type>(oldsize, adding);
        }

        // Return the number calls made so far to resolve a point
        // coordinate.  Mostly for debugging/perfing.
        size_t point_calls() const { return m_point_calls; }        


        template<typename VectorLike>
        results_type knn(size_t kay, const VectorLike& query_point) const {
            results_type ret;
            if (kay==0 || !npoints() || query_point.size() != ndim()) {
                return ret;
            }
            this->prepquery<nfkdindex_type>();

            std::vector<size_t> indices(kay,0);
            std::vector<distance_type> distances(kay, 0);
            nanoflann::KNNResultSet<element_type> nf(kay);
            nf.init(&indices[0], &distances[0]);
            m_nfkdindex->findNeighbors(nf, query_point.data(),
                                       nanoflann::SearchParameters());
            const size_t nfound = nf.size();
            ret.resize(nfound);
            for (size_t ind=0; ind<nfound; ++ind) {
                ret[ind] = std::make_pair(indices[ind], distances[ind]);
            }
            return ret;
        }


        template<typename VectorLike>
        results_type radius(distance_type rad, const VectorLike& query_point) const {
            results_type ret;
            if (rad==0 || !npoints() || query_point.size() != ndim()) {
                return ret;
            }
            this->prepquery<nfkdindex_type>();

            std::vector<nanoflann::ResultItem<size_t, element_type>> res;
            nanoflann::RadiusResultSet<element_type, size_t> rs(rad, res);
            m_nfkdindex->findNeighbors(rs, query_point.data());

            const size_t nfound = res.size();
            ret.resize(nfound);
            for (size_t ind=0; ind<nfound; ++ind) {
                const auto& one = res[ind];
                ret[ind] = std::make_pair(one.first, one.second);
            }
            return ret;
        }

        // nanoflann API.  Total number of points.
        inline size_t kdtree_get_point_count() const {
            return npoints();
        }

        // nanoflann API.  Must provide, but return false to let nanoflann
        // calculate.
        template <class BBOX>
        bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

        // nanoflann API.  Value of a point's dimension coordinate.
        inline element_type kdtree_get_pt(size_t idx, size_t dim) const {
            ++m_point_calls;
            return m_points.at(dim).at(idx);
        }

      private:
        points_type m_points;
        size_t m_nblocks{0};     // how many times we have been appended to.
        // The index is made lazily
        mutable std::unique_ptr<nfkdindex_type> m_nfkdindex;
        block_indices_type m_major_indices, m_minor_indices;

        mutable size_t m_point_calls{0};
        

        // discovery idiom to figure out how to add points given the
        // k-d tree is dynamic or static.
        template <typename T, typename ...Ts>
        using addpoints_type = decltype(std::declval<T>().addPoints(std::declval<Ts>()...));

        template<typename T>
        using has_addpoints = is_detected<addpoints_type, T, uint32_t, uint32_t>;

        // Dynamic index case
        template <class T, std::enable_if_t<has_addpoints<T>::value>* = nullptr>
        void addn(size_t beg, size_t n) {
            if (!n) {
                return;
            }
            // No query yet so no index, so nothing to add.
            if (!m_nfkdindex) {
                return;
            }
            // nanoflann takes an inclusive range of indices and NOT the usual
            // C++ begin/end convention!
            m_nfkdindex->addPoints(beg, beg+n-1);
        }
        // Dynamic
        template <class T, std::enable_if_t<has_addpoints<T>::value>* = nullptr>
        void prepquery() const {
            if (m_nfkdindex) {
                return;
            }
            m_nfkdindex = std::make_unique<nfkdindex_type>(ndim(), *this);
        }

        // Static index case
        template <class T, std::enable_if_t<!has_addpoints<T>::value>* = nullptr>
        void addn(size_t beg, size_t n) {
            // Appending after we have made a k-d tree invalidates the k-d tree.
            if (m_nfkdindex) {
                m_nfkdindex = nullptr;
            }
        }
        template <class T, std::enable_if_t<!has_addpoints<T>::value>* = nullptr>
        void prepquery() const {
            if (m_nfkdindex) {
                return;
            }
            m_nfkdindex = std::make_unique<nfkdindex_type>(ndim(), *this);
        }
    }; // Tree
} // WireCell::NFKDVec

#endif
