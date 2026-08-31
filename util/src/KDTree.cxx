#include "WireCellUtil/KDTree.h"
#include "WireCellUtil/Exceptions.h"

using namespace WireCell::KDTree;

using WireCell::raise;
using WireCell::PointCloud::Array;
using WireCell::PointCloud::Dataset;

/** This class may be used as a nanoflan "dataset adaptor" in
    order to use a "selection" of arrays from a Dataset as the
    coordinates of points.
*/
template<typename ElementType>
struct DatasetSelectionAdaptor {
    Dataset::selection_t points;
    
    using element_t = ElementType;

    explicit DatasetSelectionAdaptor(const Dataset::selection_t& selection)
        : points(selection)
    {
    }

    inline size_t kdtree_get_point_count() const
    {
        size_t num = 0;
        if (points.size() > 0) {
            num = points[0]->size_major();
        }
        return num;
    }

    inline element_t kdtree_get_pt(size_t idx, size_t dim) const
    {
        for (auto arr : points) {
            assert(arr);
            assert(arr->template is_type<ElementType>());
        }
        if (dim < points.size()) {
            auto arr = points[dim];
            if (!arr) {
                raise<WireCell::IndexError>("my array at dim %d/%d went away on index %d", dim, points.size(), idx);
            }
            if (!arr->bytes().data()) {
                raise<WireCell::IndexError>("my array data at dim %d/%d went away on index %d", dim, points.size(), idx);
            }
                
            if (idx < arr->size_major()) {
                element_t val = arr->template element<ElementType>(idx);
                return val;
            }
        }
        raise<WireCell::IndexError>("index %d dim %d out of bounds", idx, dim);
        return 0;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};                              // DatasetSelectionAdaptor


template<typename IndexType>
struct QueryStatic : public Query<typename IndexType::ElementType>
{
    using element_t = typename IndexType::ElementType;
    using dataset_adaptor_t = DatasetSelectionAdaptor<element_t>;
    using results_t = typename Query<element_t>::results_t;
    using point_t = typename Query<element_t>::point_t;

    dataset_adaptor_t m_dataset_adaptor;
    IndexType m_index;
    Metric m_metric;

    QueryStatic() = delete;
    virtual ~QueryStatic() = default;

    QueryStatic(const Dataset::selection_t& selection, Metric mtype)
        : m_dataset_adaptor(selection)
        , m_index(selection.size(), m_dataset_adaptor)
        , m_metric(mtype)
    {
        for (auto arr : selection) {
            assert(arr);
        }
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
                            nanoflann::SearchParameters());
        const size_t nfound = nf.size();
        ret.index.resize(nfound, -1);
        ret.distance.resize(nfound, -1);
        return ret;
    }


    virtual results_t radius(element_t rad,
                             const point_t& query_point)
    {
        // std::vector<std::pair<size_t, element_t>> ids;
        std::vector<nanoflann::ResultItem<size_t, element_t>> ids;
        nanoflann::RadiusResultSet<element_t, size_t> rs(rad, ids);
        m_index.findNeighbors(rs, query_point.data(),
                            nanoflann::SearchParameters());
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

    QueryDynamic(const Dataset::selection_t& selection, Metric mtype)
        : QueryStatic<IndexType>(selection, mtype)
    {
        for (auto arr : selection) {
            assert(arr);
        }
    }

    virtual void update(size_t beg, size_t end)
    {
        // Subtract 1 because NF uses a closed interval.
        this->m_index.addPoints(beg, end-1); 
    }

    virtual bool dynamic() const
    {
        return true;
    }
};

template<typename DistanceType>
std::unique_ptr<Query<typename DistanceType::ElementType>>
make_query_with_distance(Dataset& ds,
                         const Dataset::name_list_t& names,
                         bool dynamic, Metric mtype)
{
    using element_t = typename DistanceType::ElementType;
    using dataset_adaptor_t = DatasetSelectionAdaptor<element_t>;

    Dataset::selection_t selection = ds.selection(names);
    for (auto arr : selection) {
        assert(arr);
    }

    if (dynamic) {
        using index_t = nanoflann::KDTreeSingleIndexDynamicAdaptor<DistanceType, dataset_adaptor_t>;
        using query_t = QueryDynamic<index_t>;
        auto ret = std::make_unique<query_t>(selection, mtype);
        auto* raw = ret.get();
        ds.register_append([raw](size_t beg, size_t end) {
            raw->update(beg, end);
        });
        return std::move(ret);
    }

    using index_t = nanoflann::KDTreeSingleIndexAdaptor<DistanceType, dataset_adaptor_t>;
    using query_t = QueryStatic<index_t>;
    return std::make_unique<query_t>(selection, mtype);
}

template<typename ElementType>
std::unique_ptr<Query<ElementType>>
make_query(Dataset& ds,
           const Dataset::name_list_t& names,
           bool dynamic, Metric mtype)
{
    using element_t = ElementType;
    using dataset_adaptor_t = DatasetSelectionAdaptor<element_t>;

    if (mtype == Metric::l2simple) {
        using distance_t = nanoflann::L2_Simple_Adaptor<element_t, dataset_adaptor_t>;
        return make_query_with_distance<distance_t>(ds, names, dynamic, mtype);
    }
    if (mtype == Metric::l1) {
        using distance_t = nanoflann::L1_Adaptor<element_t, dataset_adaptor_t>;
        return make_query_with_distance<distance_t>(ds, names, dynamic, mtype);
    }
    if (mtype == Metric::l2) {
        using distance_t = nanoflann::L2_Adaptor<element_t, dataset_adaptor_t>;
        return make_query_with_distance<distance_t>(ds, names, dynamic, mtype);
    }
    // Metric::so2 ...
    // Metric::so3 ...

    return nullptr;
}

template<>
std::unique_ptr<Query<int>>
WireCell::KDTree::query<int>(
    Dataset& ds,
    const Dataset::name_list_t& names,
    bool dynamic, Metric mtype)
{
    return make_query<int>(ds, names, dynamic, mtype);
}


template<>
std::unique_ptr<Query<float>>
WireCell::KDTree::query<float>(
    Dataset& ds,
    const Dataset::name_list_t& names,
    bool dynamic, Metric mtype)
{
    return make_query<float>(ds, names, dynamic, mtype);
}

template<>
std::unique_ptr<Query<double>>
WireCell::KDTree::query<double>(
    Dataset& ds,
    const Dataset::name_list_t& names,
    bool dynamic, Metric mtype)
{
    return make_query<double>(ds, names, dynamic, mtype);
}

