#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type

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
    boost::hash_combine(h, wpid().ident());

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

    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.charge = pc_scalar.get("charge")->elements<float_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.center_x = pc_scalar.get("center_x")->elements<float_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.center_y = pc_scalar.get("center_y")->elements<float_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.center_z = pc_scalar.get("center_z")->elements<float_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.wpid = WirePlaneId(pc_scalar.get("wpid")->elements<int>()[0]);
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.npoints = pc_scalar.get("npoints")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.slice_index_min = pc_scalar.get("slice_index_min")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.slice_index_max = pc_scalar.get("slice_index_max")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.u_wire_index_min = pc_scalar.get("u_wire_index_min")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.u_wire_index_max = pc_scalar.get("u_wire_index_max")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.v_wire_index_min = pc_scalar.get("v_wire_index_min")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.v_wire_index_max = pc_scalar.get("v_wire_index_max")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.w_wire_index_min = pc_scalar.get("w_wire_index_min")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.w_wire_index_max = pc_scalar.get("w_wire_index_max")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.max_wire_interval = pc_scalar.get("max_wire_interval")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.min_wire_interval = pc_scalar.get("min_wire_interval")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.max_wire_type = pc_scalar.get("max_wire_type")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
    ///
    cache.min_wire_type = pc_scalar.get("min_wire_type")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change.
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
    // check apa/face
    if (wpid().apa() != b.wpid().apa()) return false;
    if (wpid().face() != b.wpid().face()) return false;
    if (u_wire_index_min() > b.u_wire_index_max()-1 + offset) return false;
    if (b.u_wire_index_min() > u_wire_index_max()-1 + offset) return false;
    if (v_wire_index_min() > b.v_wire_index_max()-1 + offset) return false;
    if (b.v_wire_index_min() > v_wire_index_max()-1 + offset) return false;
    if (w_wire_index_min() > b.w_wire_index_max()-1 + offset) return false;
    if (b.w_wire_index_min() > w_wire_index_max()-1 + offset) return false;
    return true;
}


// void Blob::check_dead_wire_consistency() const{
//     auto grouping = cluster()->grouping();

//     // Find apa and face for this blob
//     int apa = wpid().apa();
//     int face = wpid().face();

   
//     // Loop over time slices and wire indices
//     // for (int slice = slice_index_min(); slice < slice_index_max(); slice += 4) {
//     //     for (int plane = 0; plane < 3; ++plane) {
//     //         int wire_min = 0, wire_max = 0;
//     //         if (plane == 0) {
//     //             wire_min = u_wire_index_min();
//     //             wire_max = u_wire_index_max();
//     //         } else if (plane == 1) {
//     //             wire_min = v_wire_index_min();
//     //             wire_max = v_wire_index_max();
//     //         } else if (plane == 2) {
//     //             wire_min = w_wire_index_min();
//     //             wire_max = w_wire_index_max();
//     //         }
//     //         for (int wire = wire_min; wire < wire_max; ++wire) {
//     //              bool is_dead =  grouping->is_wire_dead(apa, face, plane, wire, slice);
//     //              auto charge_pair = grouping->get_wire_charge(apa, face, plane, wire, slice);

//     //             //  if (is_dead == 1) std::cout << "apa: " << apa << ", face: " << face
//     //             //         << ", slice: " << slice << ", plane: " << plane
//     //             //         << ", wire: " << wire
//     //             //         << ", is_dead: " << is_dead
//     //             //         << ", charge: (" << charge_pair.first
//     //             //         << ", " << charge_pair.second << ")" << std::endl;
//     //             //  std::cout << "APA: " << apa << ", Face: " << face
//     //             //        << ", Slice: " << slice << ", Plane: " << plane
//     //             //        << ", Wire: " << wire
//     //             //        << ", Is Dead: " << is_dead
//     //             //        << ", Charge: (" << charge_pair.first
//     //             //        << ", " << charge_pair.second << ")" << std::endl;
//     //             // You can add your logic here for each (slice, plane, wire)
//     //             // Example: LogDebug("apa=" << apa << " face=" << face << " slice=" << slice << " plane=" << plane << " wire=" << wire);
//     //         }
//     //     }
//     // }

// }

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

std::vector<geo_point_t> Blob::points(const std::string& pc_name, 
                                      const std::vector<std::string>& coords) const
{
    const auto& pc = m_node->value.local_pcs()[pc_name];
    auto sel = pc.selection(coords);
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
        const auto naf = a->wpid();
        const auto nbf = b->wpid(); 
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
    // After exhausting all "content" comparison, we are left with the question,
    // are these two blobs really different or not.  We have two choices.  We
    // may compare on pointer value which will surely "break the tie" but will
    // introduce randomness.  We may return "false" which says "these are equal"
    // in which case any unordered set/map will not hold both.  Randomness is
    // the better choice as we would have a better chance to detect that in some
    // future bug.
    return a < b;
}


double Blob::estimate_total_charge() const {
    const Cluster* cluster_ptr = this->cluster();
    if (!cluster_ptr) {
        return 0.0;
    }
    const Grouping* grouping = cluster_ptr->grouping();
    if (!grouping) {
        return 0.0;
    }
    
    double total_charge = 0.0;
    int valid_plane_count = 0;
    
    const auto wpid_val = wpid();
    const int apa = wpid_val.apa();
    const int face = wpid_val.face();
    
    // Process each plane (U=0, V=1, W=2)
    for (int plane = 0; plane < 3; plane++) {
        double plane_charge = 0.0;
        bool plane_has_data = false;
        
        // Get wire ranges for this plane
        int wire_min, wire_max;
        switch (plane) {
            case 0: // U plane
                wire_min = u_wire_index_min();
                wire_max = u_wire_index_max();
                break;
            case 1: // V plane
                wire_min = v_wire_index_min();
                wire_max = v_wire_index_max();
                break;
            case 2: // W plane
                wire_min = w_wire_index_min();
                wire_max = w_wire_index_max();
                break;
            default:
                continue;
        }
        
        // Check if this plane has valid wire ranges
        if (wire_min >= wire_max) {
            continue;
        }
        
        // Iterate through time slices for this blob
        int time_slice = slice_index_min(); 
        int num_dead_wire = 0;
        // Iterate through wires in this plane
        for (int wire_index = wire_min; wire_index < wire_max; wire_index++) {
            // Check if wire is dead
            if (grouping->is_wire_dead(apa, face, plane, wire_index, time_slice)) {
                num_dead_wire++;
            }

            // Get wire charge
            auto charge_pair = grouping->get_wire_charge(apa, face, plane, wire_index, time_slice);
            double charge = charge_pair.first;            
            if (charge > 0) { // Only count positive charges
                plane_charge += charge;
                plane_has_data = true;
            }
        }
        if (num_dead_wire > 1 || num_dead_wire == wire_max - wire_min) plane_has_data = false;
        
        if (plane_has_data) {
            total_charge += plane_charge;
            valid_plane_count++;
        }
    }
    
    // Average across valid planes (equivalent to prototype's division by count)
    if (valid_plane_count > 0) {
        total_charge /= valid_plane_count;
    }
    
    return total_charge;
}

double Blob::estimate_minimum_charge() const {
    const Cluster* cluster_ptr = this->cluster();
    if (!cluster_ptr) {
        return 1e9;
    }
    
    const Grouping* grouping = cluster_ptr->grouping();
    if (!grouping) {
        return 1e9;
    }
    
    double min_charge = 1e9;
    const auto wpid_val = wpid();
    const int apa = wpid_val.apa();
    const int face = wpid_val.face();
    
    // Process each plane (U=0, V=1, W=2)
    for (int plane = 0; plane < 3; plane++) {
        double plane_charge = 0.0;
        bool plane_has_data = false;
        
        // Get wire ranges for this plane
        int wire_min, wire_max;
        switch (plane) {
            case 0: // U plane
                wire_min = u_wire_index_min();
                wire_max = u_wire_index_max();
                break;
            case 1: // V plane
                wire_min = v_wire_index_min();
                wire_max = v_wire_index_max();
                break;
            case 2: // W plane
                wire_min = w_wire_index_min();
                wire_max = w_wire_index_max();
                break;
            default:
                continue;
        }
        
        // Check if this plane has valid wire ranges
        if (wire_min >= wire_max) {
            continue; // Skip this plane (equivalent to bad_planes check)
        }
        
        // Iterate through time slices for this blob
        int time_slice = slice_index_min(); 
        int num_dead_wire = 0;
        // Iterate through wires in this plane
        for (int wire_index = wire_min; wire_index < wire_max; wire_index++) {
            // Check if wire is dead (equivalent to bad_planes check)
            if (grouping->is_wire_dead(apa, face, plane, wire_index, time_slice)) {
                num_dead_wire++;
            }
            
            // Get wire charge
            auto charge_pair = grouping->get_wire_charge(apa, face, plane, wire_index, time_slice);
            double charge = charge_pair.first;
            
            if (charge > 0) { // Only count positive charges
                plane_charge += charge;
                plane_has_data = true;
            }
        }
        if (num_dead_wire > 1 || num_dead_wire == wire_max - wire_min) plane_has_data = false;
        
        
        // Update minimum charge if this plane has data
        if (plane_has_data && plane_charge < min_charge) {
            min_charge = plane_charge;
        }
    }
    
    return min_charge;
}

double Blob::get_wire_charge(int plane, const int_t wire_index) const {
    const Cluster* cluster_ptr = this->cluster();
    if (!cluster_ptr) {
        return 0.0;
    }
    
    const Grouping* grouping = cluster_ptr->grouping();
    if (!grouping) {
        return 0.0;
    }
    
    const auto wpid_val = wpid();
    const int apa = wpid_val.apa();
    const int face = wpid_val.face();

    // Get charge for the middle time slice of this blob as representative
    const int time_slice = slice_index_min();
    
    auto charge_pair = grouping->get_wire_charge(apa, face, plane, wire_index, time_slice);
    return charge_pair.first;
}

double Blob::get_wire_charge_error(int plane, const int_t wire_index) const {
    const Cluster* cluster_ptr = this->cluster();
    if (!cluster_ptr) {
        return 1e12; // Large error indicates no data
    }
    
    const Grouping* grouping = cluster_ptr->grouping();
    if (!grouping) {
        return 1e12; // Large error indicates no data
    }
    
    const auto wpid_val = wpid();
    const int apa = wpid_val.apa();
    const int face = wpid_val.face();
    
    // Get charge error for the middle time slice of this blob as representative
    const int time_slice = slice_index_min();
    
    auto charge_pair = grouping->get_wire_charge(apa, face, plane, wire_index, time_slice);
    return charge_pair.second;
}




void Facade::sort_blobs(std::vector<const Blob*>& blobs) { std::sort(blobs.rbegin(), blobs.rend(), blob_less); }
void Facade::sort_blobs(std::vector<Blob*>& blobs) { std::sort(blobs.rbegin(), blobs.rend(), blob_less); }

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
