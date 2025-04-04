#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/PointCloudTransform.h"

#include <vector>
// #include <iostream>             // only debug

// Implementation is totally local to this comp. unit so no need for namespacing.
class DetectorVolumes;
class T0Correction;

WIRECELL_FACTORY(DetectorVolumes, DetectorVolumes,
                 WireCell::IDetectorVolumes, WireCell::IConfigurable)


using namespace WireCell;
using WireCell::PointCloud::Dataset;
using WireCell::PointCloud::Array;

class T0Correction : public WireCell::PointCloud::Transform {
    public:
     T0Correction(IDetectorVolumes::pointer dv) : m_dv(dv) {
         for (const auto& [wfid, _] : m_dv->wpident_faces()) {
             WirePlaneId wpid(wfid);
             m_time_global_offsets[wpid.apa()][wpid.face()] = m_dv->metadata(wpid)["time_offset"].asDouble();
             m_drift_speeds[wpid.apa()][wpid.face()] = m_dv->metadata(wpid)["drift_speed"].asDouble();
         }
     }
     virtual ~T0Correction() = default;
 
     /**
      * From time2drift in Facade_Util.cxx
      * x_raw = xorig + face->dirx * (time_read_out + time_global_offset) * abs_drift_speed;
      * x_corr = xorig + face->dirx * (time_read_out - clustser_t0) * abs_drift_speed;
      *x_corr - x_raw = face->dirx * (- clustser_t0 - time_global_offset) * abs_drift_speed;
      */
     virtual Point forward(const Point &pos, double clustser_t0, int face,
                                         int apa) const override
     {
         Point result(pos);
         result[0] -= m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (clustser_t0 + m_time_global_offsets.at(apa).at(face)) *
         m_drift_speeds.at(apa).at(face);
         return result;
     }
     virtual Point backward(const Point &pos, double clustser_t0, int face,
                                          int apa) const override
     {
         Point result(pos);
         result[0] += m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (clustser_t0 + m_time_global_offsets.at(apa).at(face)) *
         m_drift_speeds.at(apa).at(face);
         return result;
     }
     virtual bool filter(const Point &pos, double clustser_t0, int face,
                         int apa) const override
     {
         return (m_dv->contained_by(pos)) ? true : false;
     }
     virtual Dataset forward(const Dataset &pc, const std::vector<std::string>& arr_names, double clustser_t0, int face,
                              int apa) const override
     {
         const auto &arr_x = pc.get(arr_names[0])->elements<double>();
         const auto &arr_y = pc.get(arr_names[1])->elements<double>();
         const auto &arr_z = pc.get(arr_names[2])->elements<double>();
         std::vector<double> arr_x_corr(arr_x.size());
         for (size_t i = 0; i < arr_x.size(); ++i) {
             arr_x_corr[i] = arr_x[i] - m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (clustser_t0 + m_time_global_offsets.at(apa).at(face)) *
             m_drift_speeds.at(apa).at(face);
         }
         Dataset ds;
         ds.add("x", Array(arr_x_corr));
         ds.add("y", Array(arr_y));
         ds.add("z", Array(arr_z));
         return ds;
     }
     virtual Dataset backward(const Dataset &pc, const std::vector<std::string>& arr_names, double clustser_t0, int face,
                               int apa) const override
     {
         const auto &arr_x = pc.get(arr_names[0])->elements<double>();
         const auto &arr_y = pc.get(arr_names[1])->elements<double>();
         const auto &arr_z = pc.get(arr_names[2])->elements<double>();
         std::vector<double> arr_x_corr(arr_x.size());
         for (size_t i = 0; i < arr_x.size(); ++i) {
             arr_x_corr[i] = arr_x[i] + m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (clustser_t0 + m_time_global_offsets.at(apa).at(face)) *
             m_drift_speeds.at(apa).at(face);
         }
         Dataset ds;
         ds.add("x", Array(arr_x_corr));
         ds.add("y", Array(arr_y));
         ds.add("z", Array(arr_z));
         return ds;
     }
     virtual Dataset filter(const Dataset &pc, const std::vector<std::string>& arr_names, double clustser_t0, int face,
                            int apa) const override
     {
         std::vector<int> arr_filter(pc.size_major());
         const auto &arr_x = pc.get(arr_names[0])->elements<double>();
         const auto &arr_y = pc.get(arr_names[1])->elements<double>();
         const auto &arr_z = pc.get(arr_names[2])->elements<double>();
         for (size_t i = 0; i < arr_x.size(); ++i) {
             arr_filter[i] = (m_dv->contained_by(Point(arr_x[i], arr_y[i], arr_z[i]))) ? 1 : 0;
            //  std::cout << "yuhw: {" << arr_x[i] << " " << arr_y[i] << " " << arr_z[i] << "} -> " << arr_filter[i] << std::endl;
         }
         Dataset ds;
         ds.add("filter", Array(arr_filter));
         return ds;
     }
 
    private:
     IDetectorVolumes::pointer m_dv; // do not own
 
     // // m_time_global_offsets.at(apa).at(face) = time_global_offset
     std::map<int, std::map<int, double>> m_time_global_offsets;
     std::map<int, std::map<int, double>> m_drift_speeds;
 };
 

class DetectorVolumes : public IDetectorVolumes,
                        public IFiducial,
                        public IConfigurable {

   public:

    DetectorVolumes()  {
    }
    virtual ~DetectorVolumes() {}

    virtual Configuration default_configuration() const {
        Configuration cfg;

        // A list of IAnodePlane "type:name" identifiers to match up to wpids. 
        cfg["anodes"] = Json::arrayValue;

        // Arbitrary, user-provided, application-specific metadata that is
        // provided site-unseen to the clients of DetectorVolumes.
        //
        // The metadata as a whole must be an object with keys formed either to
        // match WirePlaneId::name() or the literal string "default".
        //
        // The attribute values are not determined by DetectorVolumes but by the
        // consuming application (the components for which the user has
        // configured to use the DetectorVolumes).
        //
        // If an attribute key is not found that matches wpid.name() of the
        // query then the attribute value "default" is returned.  If both do not
        // exist then a "null" JSON value is returned.
        cfg["metadata"] = Json::objectValue;

        // Add grid configuration parameters with defaults
        cfg["grid_size_x"] = 0.0;  // If 0, will be determined automatically
        cfg["grid_size_y"] = 0.0;  // If 0, will be determined automatically
        cfg["grid_size_z"] = 0.0;  // If 0, will be determined automatically
        cfg["grid_fraction"] = 0.5; // Default fraction of average bounding box size

        return cfg;
    }

    // Return a new wpid with layer masked to 0 so it represents only
    // anode+face.  This is the value used to index m_faces.
    WirePlaneId wfid(WirePlaneId wpid) const {
        return WirePlaneId(WirePlaneLayer_t::kAllLayers, wpid.face(), wpid.apa());
    }

    virtual void configure(const Configuration& cfg) {
        m_faces.clear();
        
        // Reset spatial structures
        m_overall_bb = BoundingBox();
        m_grid.clear();

        // Store grid configuration parameters
        m_config_grid_size_x = cfg.get("grid_size_x", 0.0).asDouble();
        m_config_grid_size_y = cfg.get("grid_size_y", 0.0).asDouble();
        m_config_grid_size_z = cfg.get("grid_size_z", 0.0).asDouble();
        m_grid_fraction = cfg.get("grid_fraction", 0.5).asDouble();

        for (const auto& janode : cfg["anodes"]) {
            const std::string anode_tn = janode.asString();
            auto anode = Factory::find_tn<IAnodePlane>(anode_tn);

            for (auto iface : anode->faces()) {
                if (! iface) continue;
                
                auto planes = iface->planes();
                auto wpid = wfid(planes[0]->planeid());
                if (!wpid) {
                    raise<ValueError>("got bogus wpid from anode %s", anode_tn);
                }

                // std::cerr << "DetectorVolumes face for: " << wpid << "\n";
                m_faces[wpid.ident()] = iface;

                // Update overall bounding box with this face's sensitive volume
                BoundingBox face_bb = iface->sensitive();
                m_overall_bb(face_bb.bounds());
            }
        }
        initialize_spatial_queries();

        m_md = cfg["metadata"];

        m_transforms["T0Correction"] = std::make_shared<T0Correction>(IDetectorVolumes::pointer(this));

        Json::FastWriter fastWriter;
        SPDLOG_TRACE("metadata: {}", fastWriter.write(m_md));
    }


    // IFiducial
    virtual bool contained(const Point& point) const {
        return contained_by(point);
    }


    // // Rest is IDetectorVolumes
    // virtual WirePlaneId contained_by(const Point& point) const {
    //     // This initial imp is perhaps too slow.  There are two options I can
    //     // think of immediately:
    //     // 
    //     // 1) Try to divine a way to represent the BBs on a regular grid of
    //     // boxes.  Calculate the 3D grid coordinates of a point directly, eg
    //     // i=floor((x-o)/w), etc for j/y and k/z.  Use ijk to look up iface to
    //     // do BB.inside() test.
    //     //
    //     // 2) Perhaps simpler, construct a k-d tree with BB corners.  Query to
    //     // find closest corner to point.  Associate corner back to iface and to
    //     // BB.inside() test.

    //     for (const auto& [wpident, iface] : m_faces) {
    //         auto bb = iface->sensitive();
    //         if (bb.inside(point)) {
    //             return WirePlaneId(wpident);
    //         }
    //     }
    //     return WirePlaneId(WirePlaneLayer_t::kUnknownLayer, -1, -1);
    // }
    
    IAnodeFace::pointer get_face(WirePlaneId wpid) const {
        wpid = wfid(wpid);
        if (!wpid) {
            // std::cerr << "get_face false wpid: " << wpid << std::endl;
            return nullptr; 
        }
        auto it = m_faces.find(wpid.ident());
        if (it == m_faces.end()) {
            // std::cerr << "get_face no face for wpid: " << wpid << std::endl;
            return nullptr;
        }
        return it->second;
    }        

    IWirePlane::pointer get_plane(WirePlaneId wpid) const {
        if (! wpid.valid()) {
            // std::cerr << "get_plane invalid wpid: " << wpid << std::endl;
            return nullptr;
        }
        auto iface = get_face(wpid);
        if (!iface) {
            return nullptr;
        }
        return iface->planes()[wpid.index()];
    }

    virtual int face_dirx(WirePlaneId wpid) const {
        auto iface = get_face(wpid);
        if (!iface) return 0;
        return iface->dirx();
    }

    virtual Vector wire_direction(WirePlaneId wpid) const {
        auto iplane = get_plane(wpid);
        if (! iplane) { return Vector(0,0,0); }
        return iplane->pimpos()->axis(1);
    }
    
    virtual Vector pitch_vector(WirePlaneId wpid) const {
        auto iplane = get_plane(wpid);
        if (! iplane) { return Vector(0,0,0); }
        const auto* pimpos = iplane->pimpos();
        auto pdir = pimpos->axis(2);
        auto pmag = pimpos->region_binning().binsize();
        return pmag*pdir;
    }

    virtual BoundingBox inner_bounds(WirePlaneId wpid) const {
        auto iface = get_face(wpid);
        if (iface) {
            return iface->sensitive();
        }
        return BoundingBox();
    }

    /// Forward any user-provided, application specific metadata for a
    /// particular wpid.  
    /// TODO: use wpid.ident() = 0 for overall metadata?
    virtual Configuration metadata(WirePlaneId wpid) const {
        const auto key = wpid.ident() == 0? "overall" : wpid.name();
        if (m_md.isNull()) {
            return Json::nullValue;
        }
        if (! m_md.isMember(key)) {
            return m_md["default"];
        }
        return m_md[key];
    }

    // Initialize spatial data structures for efficient queries
    virtual bool initialize_spatial_queries() {

        // Calculate grid dimensions based on overall bounding box
        Point min_point = m_overall_bb.bounds().first;
        Point max_point = m_overall_bb.bounds().second;
        
        // Store origin for grid calculations
        m_grid_origin = min_point;
        
        // Calculate average bounding box size if auto grid sizing is requested
        if (m_config_grid_size_x <= 0.0 || m_config_grid_size_y <= 0.0 || m_config_grid_size_z <= 0.0) {
            Vector avg_size(0, 0, 0);
            int count = 0;
            for (const auto& [wpident, iface] : m_faces) {
                BoundingBox face_bb = iface->sensitive();
                if (!face_bb.empty()) {
                    Vector size = face_bb.dimensions();
                    avg_size = avg_size + size;
                    count++;
                }
            }
            
            if (count > 0) {
                avg_size = avg_size * (1.0 / count);
                // Set grid cell size to a fraction of average bounding box size
                if (m_config_grid_size_x <= 0.0) {
                    m_grid_size_x = avg_size.x() * m_grid_fraction;
                } else {
                    m_grid_size_x = m_config_grid_size_x;
                }
                
                if (m_config_grid_size_y <= 0.0) {
                    m_grid_size_y = avg_size.y() * m_grid_fraction;
                } else {
                    m_grid_size_y = m_config_grid_size_y;
                }
                
                if (m_config_grid_size_z <= 0.0) {
                    m_grid_size_z = avg_size.z() * m_grid_fraction;
                } else {
                    m_grid_size_z = m_config_grid_size_z;
                }
            } else {
                // Fallback if no valid bounding boxes
                m_grid_size_x = 1.0;
                m_grid_size_y = 1.0;
                m_grid_size_z = 1.0;
            }
        } else {
            // Use configured sizes
            m_grid_size_x = m_config_grid_size_x;
            m_grid_size_y = m_config_grid_size_y;
            m_grid_size_z = m_config_grid_size_z;
        }
        
        // Ensure grid cell sizes are positive
        m_grid_size_x = std::max(m_grid_size_x, 0.001*units::m);
        m_grid_size_y = std::max(m_grid_size_y, 0.001*units::m);
        m_grid_size_z = std::max(m_grid_size_z, 0.001*units::m);
        
        // Calculate number of cells in each dimension
        int nx = std::ceil((max_point.x() - min_point.x()) / m_grid_size_x);
        int ny = std::ceil((max_point.y() - min_point.y()) / m_grid_size_y);
        int nz = std::ceil((max_point.z() - min_point.z()) / m_grid_size_z);
        
        // Limit grid cell count to reasonable bounds
        const int min_grid_cells = 5;    // Minimum useful number of cells per dimension
        const int max_grid_cells = 100;  // Maximum reasonable number of cells per dimension
        
        nx = std::max(min_grid_cells, std::min(nx, max_grid_cells));
        ny = std::max(min_grid_cells, std::min(ny, max_grid_cells));
        nz = std::max(min_grid_cells, std::min(nz, max_grid_cells));
        
        // Recalculate grid cell size based on new cell counts
        m_grid_size_x = (max_point.x() - min_point.x()) / nx;
        m_grid_size_y = (max_point.y() - min_point.y()) / ny;
        m_grid_size_z = (max_point.z() - min_point.z()) / nz;
        
        // Ensure at least one cell in each dimension
        nx = std::max(nx, 1);
        ny = std::max(ny, 1);
        nz = std::max(nz, 1);
        
        // Resize the grid to hold all cells
        m_grid.resize(nx);
        for (auto& grid_x : m_grid) {
            grid_x.resize(ny);
            for (auto& grid_y : grid_x) {
                grid_y.resize(nz);
            }
        }
        
        // Populate the grid with face IDs
        for (const auto& [wpident, iface] : m_faces) {
            BoundingBox face_bb = iface->sensitive();
            if (face_bb.empty()) continue;
            
            // Calculate grid cell indices for this face's bounding box
            Point bb_min = face_bb.bounds().first;
            Point bb_max = face_bb.bounds().second;
            
            // Add a small epsilon to ensure boundary boxes are included in all relevant cells
            double epsilon = 1e-6*units::m;  // Small safety margin
            
            int min_i = std::max(0, static_cast<int>((bb_min.x() - min_point.x() - epsilon) / m_grid_size_x));
            int min_j = std::max(0, static_cast<int>((bb_min.y() - min_point.y() - epsilon) / m_grid_size_y));
            int min_k = std::max(0, static_cast<int>((bb_min.z() - min_point.z() - epsilon) / m_grid_size_z));
            
            int max_i = std::min(nx - 1, static_cast<int>((bb_max.x() - min_point.x() + epsilon) / m_grid_size_x));
            int max_j = std::min(ny - 1, static_cast<int>((bb_max.y() - min_point.y() + epsilon) / m_grid_size_y));
            int max_k = std::min(nz - 1, static_cast<int>((bb_max.z() - min_point.z() + epsilon) / m_grid_size_z));
            
            // Add the face ID to all cells that overlap with its bounding box
            for (int i = min_i; i <= max_i; ++i) {
                for (int j = min_j; j <= max_j; ++j) {
                    for (int k = min_k; k <= max_k; ++k) {
                        m_grid[i][j][k].push_back(wpident);
                    }
                }
            }
        }
        
        return true;
    }

    // Optimized implementation using spatial grid
    virtual WirePlaneId contained_by(const Point& point) const {
        // Quick check if point is in overall volume first
        if (!m_overall_bb.inside(point)) {
            return WirePlaneId(WirePlaneLayer_t::kUnknownLayer, -1, -1);
        }
        
        // Calculate grid cell indices for the point
        int i = static_cast<int>((point.x() - m_grid_origin.x()) / m_grid_size_x);
        int j = static_cast<int>((point.y() - m_grid_origin.y()) / m_grid_size_y);
        int k = static_cast<int>((point.z() - m_grid_origin.z()) / m_grid_size_z);
            
        // Check if indices are within grid bounds
        if (i >= 0 && i < static_cast<int>(m_grid.size()) &&
            j >= 0 && j < static_cast<int>(m_grid[0].size()) &&
            k >= 0 && k < static_cast<int>(m_grid[0][0].size())) {
                
            // Check faces in this grid cell
            for (int wpident : m_grid[i][j][k]) {
                auto it = m_faces.find(wpident);
                if (it != m_faces.end()) {
                    auto bb = it->second->sensitive();
                    if (bb.inside(point)) {
                        return WirePlaneId(wpident);
                    }
                }
            }
                
            // Special handling for points near cell boundaries
            // Check neighboring cells if we're close to an edge
            double eps = 1e-6*units::m;  // Small distance from boundary to check neighbors
            std::vector<std::tuple<int,int,int>> neighbors;
                
            // Check if we're close to x boundaries
            double x_frac = fmod((point.x() - m_grid_origin.x()) / m_grid_size_x, 1.0);
            if (x_frac < eps || x_frac > (1.0 - eps)) {
                // Add x-neighbors
                if (i > 0) neighbors.push_back({i-1, j, k});
                if (i < static_cast<int>(m_grid.size()) - 1) neighbors.push_back({i+1, j, k});
            }
                
            // Check if we're close to y boundaries
            double y_frac = fmod((point.y() - m_grid_origin.y()) / m_grid_size_y, 1.0);
            if (y_frac < eps || y_frac > (1.0 - eps)) {
                // Add y-neighbors
                if (j > 0) neighbors.push_back({i, j-1, k});
                if (j < static_cast<int>(m_grid[0].size()) - 1) neighbors.push_back({i, j+1, k});
            }
                
            // Check if we're close to z boundaries
            double z_frac = fmod((point.z() - m_grid_origin.z()) / m_grid_size_z, 1.0);
            if (z_frac < eps || z_frac > (1.0 - eps)) {
                // Add z-neighbors
                if (k > 0) neighbors.push_back({i, j, k-1});
                if (k < static_cast<int>(m_grid[0][0].size()) - 1) neighbors.push_back({i, j, k+1});
            }
                
            // Check neighboring cells
            for (auto [ni, nj, nk] : neighbors) {
                for (int wpident : m_grid[ni][nj][nk]) {
                    auto it = m_faces.find(wpident);
                    if (it != m_faces.end()) {
                        auto bb = it->second->sensitive();
                        if (bb.inside(point)) {
                            return WirePlaneId(wpident);
                        }
                    }
                }
            }
        }
            
        return WirePlaneId(WirePlaneLayer_t::kUnknownLayer, -1, -1);
    }


    virtual const std::map<int, IAnodeFace::pointer>& wpident_faces() const {
        return m_faces;
    }

    const WireCell::PointCloud::Transform::pointer pc_transform(const std::string& transform_name) const {
        auto it = m_transforms.find(transform_name);
        if (it == m_transforms.end()) {
            raise<KeyError>("No transform found for name: %s", transform_name.c_str());
        }
        return it->second;
    }



private:
    // A map of transform names to their respective transform objects
    std::map<std::string, WireCell::PointCloud::Transform::pointer> m_transforms;

    // Map wpid with layer=0 to its face.
    std::map<int, IAnodeFace::pointer> m_faces;

    Configuration m_md;
    // Overall bounding box containing all sensitive volumes
    BoundingBox m_overall_bb;
    
    // Grid configuration parameters
    double m_config_grid_size_x;
    double m_config_grid_size_y;
    double m_config_grid_size_z;
    double m_grid_fraction;
    
    // Actual grid cell sizes (may be auto-calculated)
    double m_grid_size_x;
    double m_grid_size_y;
    double m_grid_size_z;
    
    // Grid origin point (minimum corner of overall bounding box)
    Point m_grid_origin;
    
    // 3D grid of face IDs for spatial lookups
    // Format: [x][y][z] -> list of wpidents that overlap this cell
    std::vector<std::vector<std::vector<std::vector<int>>>> m_grid;
};


