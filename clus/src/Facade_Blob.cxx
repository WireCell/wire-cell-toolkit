#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

std::ostream& Facade::operator<<(std::ostream& os, const Facade::Blob& blob)
{
    os << "<Blob [" << (void*) blob.hash() << "]:" << " npts=" << blob.npoints() << " r=" << blob.center_pos()
       << " q=" << blob.charge() << " t=[" << blob.slice_index_min() << "," << blob.slice_index_max() << "]" << " u=["
       << blob.u_wire_index_min() << "," << blob.u_wire_index_max() << "]" << " v=[" << blob.v_wire_index_min() << ","
       << blob.v_wire_index_max() << "]" << " w=[" << blob.w_wire_index_min() << "," << blob.w_wire_index_max() << "]"
       << ">";
    return os;
}

Cluster* Blob::cluster() { return this->m_node->parent->value.template facade<Cluster>(); }
const Cluster* Blob::cluster() const { return this->m_node->parent->value.template facade<Cluster>(); }

size_t Blob::hash() const
{
    std::size_t h = 0;
    boost::hash_combine(h, npoints());
    // boost::hash_combine(h, charge());
    boost::hash_combine(h, center_x());
    boost::hash_combine(h, center_y());
    boost::hash_combine(h, center_z());
    boost::hash_combine(h, face());

    boost::hash_combine(h, slice_index_min());
    boost::hash_combine(h, slice_index_max());
    boost::hash_combine(h, u_wire_index_min());
    boost::hash_combine(h, u_wire_index_max());
    boost::hash_combine(h, v_wire_index_min());
    boost::hash_combine(h, v_wire_index_max());
    boost::hash_combine(h, w_wire_index_min());
    boost::hash_combine(h, w_wire_index_max());
    return h;
}

void Blob::fill_cache(BlobCache& cache) const
{
    const auto& pc_scalar = get_pc("scalar");

    if (pc_scalar.size_major() != 1) {
        raise<ValueError>("scalar PC is not scalar but size %d", pc_scalar.size_major());
    }

    cache.charge = pc_scalar.get("charge")->elements<float_t>()[0];
    cache.center_x = pc_scalar.get("center_x")->elements<float_t>()[0];
    cache.center_y = pc_scalar.get("center_y")->elements<float_t>()[0];
    cache.center_z = pc_scalar.get("center_z")->elements<float_t>()[0];
    cache.face = pc_scalar.get("face")->elements<int_t>()[0];
    cache.npoints = pc_scalar.get("npoints")->elements<int_t>()[0];
    cache.slice_index_min = pc_scalar.get("slice_index_min")->elements<int_t>()[0];
    cache.slice_index_max = pc_scalar.get("slice_index_max")->elements<int_t>()[0];
    cache.u_wire_index_min = pc_scalar.get("u_wire_index_min")->elements<int_t>()[0];
    cache.u_wire_index_max = pc_scalar.get("u_wire_index_max")->elements<int_t>()[0];
    cache.v_wire_index_min = pc_scalar.get("v_wire_index_min")->elements<int_t>()[0];
    cache.v_wire_index_max = pc_scalar.get("v_wire_index_max")->elements<int_t>()[0];
    cache.w_wire_index_min = pc_scalar.get("w_wire_index_min")->elements<int_t>()[0];
    cache.w_wire_index_max = pc_scalar.get("w_wire_index_max")->elements<int_t>()[0];
    cache.max_wire_interval = pc_scalar.get("max_wire_interval")->elements<int_t>()[0];
    cache.min_wire_interval = pc_scalar.get("min_wire_interval")->elements<int_t>()[0];
    cache.max_wire_type = pc_scalar.get("max_wire_type")->elements<int_t>()[0];
    cache.min_wire_type = pc_scalar.get("min_wire_type")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change the above.
    ///
    const auto& lpcs = m_node->value.local_pcs();
    const auto& pc_corner = lpcs.at("corner");
    const auto& x = pc_corner.get("x")->elements<float_t>();
    const auto& y = pc_corner.get("y")->elements<float_t>();
    const auto& z = pc_corner.get("z")->elements<float_t>();
    const size_t size = x.size();
    cache.corners_.resize(size);
    for (size_t ind = 0; ind < size; ++ind) {
        cache.corners_[ind] = {x[ind], y[ind], z[ind]};
        // std::cout << "corner " << cache.corners_[ind] << std::endl;
    }
}

bool Blob::overlap_fast(const Blob& b, const int offset) const
{
    // check face ...
    if (face() != b.face()) return false;
    if (u_wire_index_min() > b.u_wire_index_max()-1 + offset) return false;
    if (b.u_wire_index_min() > u_wire_index_max()-1 + offset) return false;
    if (v_wire_index_min() > b.v_wire_index_max()-1 + offset) return false;
    if (b.v_wire_index_min() > v_wire_index_max()-1 + offset) return false;
    if (w_wire_index_min() > b.w_wire_index_max()-1 + offset) return false;
    if (b.w_wire_index_min() > w_wire_index_max()-1 + offset) return false;
    return true;
}

geo_point_t Blob::center_pos() const
{
    return {cache().center_x, cache().center_y, cache().center_z};
}

size_t Blob::nbpoints() const
{
    const auto& pc = get_pc("3d");
    return pc.size_major();
}

bool Blob::sanity(Log::logptr_t log) const
{
    if (nbpoints() == (size_t) npoints()) return true;
    if (log) log->debug("blob sanity: blob points mismatch: {}", *this);
    return false;
}

std::vector<geo_point_t> Blob::points() const
{
    const auto& pc = m_node->value.local_pcs()["3d"];
    auto sel = pc.selection({"x", "y", "z"});
    const size_t npts = sel[0]->size_major();

    std::vector<geo_point_t> ret(npts);
    for (int dim = 0; dim < 3; ++dim) {
        auto coord = sel[dim]->elements<double>();
        for (size_t ind = 0; ind < npts; ++ind) {
            ret[ind][dim] = coord[ind];
        }
    }
    return ret;
}

bool Facade::blob_less(const Facade::Blob* a, const Facade::Blob* b)
{
    if (a == b) return false;
    {
        const auto naf = a->face();
        const auto nbf = b->face(); 
        if (naf < nbf) return true;
        if (nbf < naf) return false;
    }
    {
        const auto na = a->npoints();
        const auto nb = b->npoints();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->charge();
        const auto nb = b->charge();
        if (na < nb) return true;
        if (nb < na) return false;
    }

    {
        const auto na = a->slice_index_min();
        const auto nb = b->slice_index_min();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->slice_index_max();
        const auto nb = b->slice_index_max();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->u_wire_index_min();
        const auto nb = b->u_wire_index_min();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->v_wire_index_min();
        const auto nb = b->v_wire_index_min();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->w_wire_index_min();
        const auto nb = b->w_wire_index_min();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->u_wire_index_max();
        const auto nb = b->u_wire_index_max();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->v_wire_index_max();
        const auto nb = b->v_wire_index_max();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const auto na = a->w_wire_index_max();
        const auto nb = b->w_wire_index_max();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    return a < b;
}



void Facade::sort_blobs(std::vector<const Blob*>& blobs) { std::sort(blobs.rbegin(), blobs.rend(), blob_less); }
void Facade::sort_blobs(std::vector<Blob*>& blobs) { std::sort(blobs.rbegin(), blobs.rend(), blob_less); }

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
