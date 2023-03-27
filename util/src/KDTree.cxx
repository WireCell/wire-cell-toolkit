#include "WireCellUtil/KDTree.h"

#include "nanoflann.hpp"

using namespace WireCell::KDTree;

using WireCell::PointCloud::selection_t;
using WireCell::PointCloud::name_list_t;
using WireCell::PointCloud::Array;
using WireCell::PointCloud::Dataset;

/** This class may be used as a nanoflan "dataset adaptor" in
    order to use a "selection" of arrays from a Dataset as the
    coordinates of points.
*/
template<typename ElementType>
class NanoflannAdaptor {
    selection_t m_pts;
  public:
    using element_t = ElementType;

    /** Construct with a reference to a dataset and the ordered
        list of arrays in the dataset to use as the coordinates.
    */
    explicit NanoflannAdaptor(const selection_t& points)
        : m_pts(points)
    {
    }
    /// No default constructor
    NanoflannAdaptor() = delete;
    /// Copies and moves are okay
    NanoflannAdaptor(const NanoflannAdaptor&) = default;
    NanoflannAdaptor(NanoflannAdaptor&&) = default;
    NanoflannAdaptor& operator=(const NanoflannAdaptor&) = default;
    NanoflannAdaptor& operator=(NanoflannAdaptor&&) = default;
    ~NanoflannAdaptor() = default;

    inline size_t kdtree_get_point_count() const
    {
        if (m_pts.empty()) {
            return 0;
        }
        return m_pts[0].get().size_major();
    }

    inline element_t kdtree_get_pt(size_t idx, size_t dim) const
    {
        if (dim < m_pts.size()) {
            const Array& arr = m_pts[dim];
            if (idx < arr.size_major()) {
                return arr.element<ElementType>(idx);
            }
        }
        THROW(WireCell::IndexError() << WireCell::errmsg{"index out of bounds"});
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};                              // NanoflannAdaptor


template<typename IndexType>
struct QueryStatic : public Query<typename IndexType::ElementType>
{
    using element_t = typename IndexType::ElementType;
    using distance_t = typename IndexType::ElementType;
    using dataset_adaptor_t = NanoflannAdaptor<element_t>;
    using results_t = typename Query<element_t>::results_t;
    using point_t = typename Query<element_t>::point_t;

    dataset_adaptor_t m_dataset_adaptor;
    name_list_t m_selection;
    IndexType m_index;
    Metric m_metric;

    QueryStatic() = delete;
    virtual ~QueryStatic() = default;

    QueryStatic(Dataset& dataset, const name_list_t& selection, Metric mtype)
        : m_dataset_adaptor(dataset.selection(selection))
        , m_selection(selection)
        , m_index(selection.size(), m_dataset_adaptor)
        , m_metric(mtype)
    {
    }

    virtual const name_list_t& selection() const
    {
        return m_selection;
    }
    virtual bool dynamic() const
    {
        return false;
    }
    virtual Metric metric() const
    {
        return m_metric;
    }

    // knn and radius have semingly capricious differences in
    // interface!  Why 2 arrays for knn and vec<pair> for rad?

    virtual results_t knn(size_t kay,
                          const point_t& query_point)
    {
        results_t ret(kay);
        nanoflann::KNNResultSet<element_t> nf(kay);
        nf.init(&ret.index[0], &ret.distance[0]);
        m_index.findNeighbors(nf, query_point.data(),
                            nanoflann::SearchParams());
        const size_t nfound = nf.size();
        ret.index.resize(nfound, -1);
        ret.distance.resize(nfound, -1);
        return ret;
    }


    virtual results_t radius(element_t rad,
                             const point_t& query_point)
    {
        std::vector<std::pair<size_t, element_t>> ids;
        nanoflann::RadiusResultSet<element_t> nf(rad, ids);
        m_index.findNeighbors(nf, query_point.data(),
                            nanoflann::SearchParams());
        const size_t nfound = ids.size();
        results_t ret(nfound);
        for (size_t ifound=0; ifound<nfound; ++ifound) {
            const auto& [indx, dist] = ids[ifound];
            ret.index[ifound] = indx;
            ret.distance[ifound] = dist;
        }
        return ret;
    }

};

template<typename IndexType>
struct QueryDynamic : public QueryStatic<IndexType>
{
  public:

    QueryDynamic() = delete;
    virtual ~QueryDynamic() = default;

    QueryDynamic(Dataset& dataset, const name_list_t& selection, Metric mtype)
        : QueryStatic<IndexType>(dataset, selection, mtype)
    {
        dataset.register_append([this](size_t beg, size_t end) {
            // NOTE: we subtract one as nanoflann expects a closed
            // range not the standard C++ right-open range.
            this->m_index.addPoints(beg, end-1); });
    }

    virtual bool dynamic() const
    {
        return true;
    }
};


template<typename DistanceType>
std::unique_ptr<Query<typename DistanceType::ElementType>>
make_query_with_distance(Dataset& dataset,
                         const name_list_t& selection,
                         bool dynamic, Metric mtype)
{
    using element_t = typename DistanceType::ElementType;
    using dataset_adaptor_t = NanoflannAdaptor<element_t>;

    if (dynamic) {
        using index_t = nanoflann::KDTreeSingleIndexDynamicAdaptor<DistanceType, dataset_adaptor_t>;
        using query_t = QueryDynamic<index_t>;
        return std::make_unique<query_t>(dataset, selection, mtype);
    }
    using index_t = nanoflann::KDTreeSingleIndexAdaptor<DistanceType, dataset_adaptor_t>;
    using query_t = QueryStatic<index_t>;
    return std::make_unique<query_t>(dataset, selection, mtype);
}

template<typename ElementType>
std::unique_ptr<Query<ElementType>>
make_query(Dataset& dataset, const name_list_t& selection,
           bool dynamic, Metric mtype)
{
    using element_t = ElementType;
    using dataset_adaptor_t = NanoflannAdaptor<element_t>;

    if (mtype == Metric::l2simple) {
        using distance_t = nanoflann::L2_Simple_Adaptor<element_t, dataset_adaptor_t>;
        return make_query_with_distance<distance_t>(dataset, selection, dynamic, mtype);
    }
    if (mtype == Metric::l1) {
        using distance_t = nanoflann::L1_Adaptor<element_t, dataset_adaptor_t>;
        return make_query_with_distance<distance_t>(dataset, selection, dynamic, mtype);
    }
    if (mtype == Metric::l2) {
        using distance_t = nanoflann::L2_Adaptor<element_t, dataset_adaptor_t>;
        return make_query_with_distance<distance_t>(dataset, selection, dynamic, mtype);
    }
    // if (mtype == Metric::so2) {
    //     using distance_t = nanoflann::SO2_Adaptor<element_t, dataset_adaptor_t>;
    //     return make_query_with_distance<distance_t>(dataset, selection, dynamic);
    // }
    // if (mtype == Metric::so3) {
    //     using distance_t = nanoflann::SO3_Adaptor<element_t, dataset_adaptor_t>;
    //     return make_query_with_distance<distance_t>(dataset, selection, dynamic);
    // }

    return nullptr;
}

template<>
std::unique_ptr<Query<int>>
WireCell::KDTree::query<int>(
    Dataset& dataset, const name_list_t& selection,
    bool dynamic, Metric mtype)
{
    return make_query<int>(dataset, selection, dynamic, mtype);
}


template<>
std::unique_ptr<Query<float>>
WireCell::KDTree::query<float>(
    Dataset& dataset, const name_list_t& selection,
    bool dynamic, Metric mtype)
{
    return make_query<float>(dataset, selection, dynamic, mtype);
}

template<>
std::unique_ptr<Query<double>>
WireCell::KDTree::query<double>(
    Dataset& dataset, const name_list_t& selection,
    bool dynamic, Metric mtype)
{
    return make_query<double>(dataset, selection, dynamic, mtype);
}

