#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellClus/DynamicPointCloud.h"
namespace WireCell::Clus::PR {

 
    // Default initialization constructor following WCPPID WCShower logic
    Shower::Shower(Graph& graph)
        : TrajectoryView(graph)
        , m_full_graph(graph)
    {
        // Initialize all ShowerData members to defaults
        data.particle_type = 0;
        data.kenergy_range = 0;
        data.kenergy_dQdx = 0;
        data.kenergy_charge = 0;
        data.kenergy_best = 0;
        
        data.start_point = Point(0, 0, 0);
        data.end_point = Point(0, 0, 0);
        data.init_dir = Vector(0, 0, 0);
        
        data.start_connection_type = 0;
        
        // Initialize start vertex and segment to nullptr
        m_start_vertex = nullptr;
        m_start_segment = nullptr;
        
        // Set shower flag (following WCPPID)
        set_flags(ShowerFlags::kShower);
        unset_flags(ShowerFlags::kKinematics);
    }

    Shower::~Shower()
    {
    }


    VertexPtr Shower::start_vertex()
    {
        return m_start_vertex;
    }

    SegmentPtr Shower::start_segment()
    {
        return m_start_segment;
    }


    // Chainable setters

    /// Chainable setter of start vertex
    Shower& Shower::set_start_vertex(VertexPtr vtx, int type)
    {
        if (! vtx->descriptor_valid()) {
            m_start_vertex = nullptr;
            return *this;
        }
        this->add_vertex(vtx);
        m_start_vertex = vtx;
        data.start_connection_type = type;


        return *this;
    }

    
    
    /// Chainable setter of start segment
    Shower& Shower::set_start_segment(SegmentPtr seg, bool flag_include_vertices, const std::string& cloud_name_fit, const std::string& cloud_name_associate)
    {
        if (! seg->descriptor_valid()) {
            m_start_segment = nullptr;
            return *this;
        }
        this->add_segment(seg);
        m_start_segment = seg;
        
        // If flag_include_vertices is true, add all vertices connected to this segment
        if (flag_include_vertices) {
            // Get the two vertices connected to this segment from the full graph
            auto vertices = find_vertices(m_full_graph, seg);
            
            // Add each vertex to the view (skip the start_vertex)
            if (vertices.first && vertices.first != m_start_vertex) {
                this->add_vertex(vertices.first);
            }
            if (vertices.second && vertices.second != m_start_vertex) {
                this->add_vertex(vertices.second);
            }
        }
        
        // Merge dynamic point clouds from segment to shower with the provided names
        if (!cloud_name_fit.empty()) {
            auto seg_dpc_fit = seg->dpcloud(cloud_name_fit);
            if (seg_dpc_fit) {
                auto shower_dpc_fit = this->dpcloud(cloud_name_fit);
                if (!shower_dpc_fit) {
                    // Create new DPC if it doesn't exist in shower, copying wpid_params from segment's DPC
                    // We need the wpid_params to construct a new DPC, but it's private
                    // For now, just share the pointer - we'll need to modify this if independent DPCs are needed
                    this->dpcloud(cloud_name_fit, seg_dpc_fit);
                } else {
                    // Add points from segment's DPC to existing shower's DPC
                    shower_dpc_fit->add_points(seg_dpc_fit->get_points());
                }
            }
        }
        
        if (!cloud_name_associate.empty()) {
            auto seg_dpc_associate = seg->dpcloud(cloud_name_associate);
            if (seg_dpc_associate) {
                auto shower_dpc_associate = this->dpcloud(cloud_name_associate);
                if (!shower_dpc_associate) {
                    // Create new DPC if it doesn't exist in shower
                    this->dpcloud(cloud_name_associate, seg_dpc_associate);
                } else {
                    // Add points from segment's DPC to existing shower's DPC
                    shower_dpc_associate->add_points(seg_dpc_associate->get_points());
                }
            }
        }
        
        return *this;
    }

    void Shower::add_segment(SegmentPtr seg, bool flag_include_vertices, const std::string& cloud_name_fit, const std::string& cloud_name_associate)
    {
        if (! seg->descriptor_valid()) {
            return;
        }
        this->add_segment(seg);
        
        // If flag_include_vertices is true, add all vertices connected to this segment
        if (flag_include_vertices) {
            // Get the two vertices connected to this segment from the full graph
            auto vertices = find_vertices(m_full_graph, seg);
            
            // Add each vertex to the view (skip the start_vertex)
            if (vertices.first ) {
                this->add_vertex(vertices.first);
            }
            if (vertices.second ) {
                this->add_vertex(vertices.second);
            }
        }

        // Merge dynamic point clouds from segment to shower with the provided names
        if (!cloud_name_fit.empty()) {
            auto seg_dpc_fit = seg->dpcloud(cloud_name_fit);
            if (seg_dpc_fit) {
                auto shower_dpc_fit = this->dpcloud(cloud_name_fit);
                if (!shower_dpc_fit) {
                    // Create new DPC if it doesn't exist in shower, copying wpid_params from segment's DPC
                    // We need the wpid_params to construct a new DPC, but it's private
                    // For now, just share the pointer - we'll need to modify this if independent DPCs are needed
                    this->dpcloud(cloud_name_fit, seg_dpc_fit);
                } else {
                    // Add points from segment's DPC to existing shower's DPC
                    shower_dpc_fit->add_points(seg_dpc_fit->get_points());
                }
            }
        }
        
        if (!cloud_name_associate.empty()) {
            auto seg_dpc_associate = seg->dpcloud(cloud_name_associate);
            if (seg_dpc_associate) {
                auto shower_dpc_associate = this->dpcloud(cloud_name_associate);
                if (!shower_dpc_associate) {
                    // Create new DPC if it doesn't exist in shower
                    this->dpcloud(cloud_name_associate, seg_dpc_associate);
                } else {
                    // Add points from segment's DPC to existing shower's DPC
                    shower_dpc_associate->add_points(seg_dpc_associate->get_points());
                }
            }
        }
    }

    void Shower::set_flag_kinematics(bool val){
        if (val) {
            set_flags(ShowerFlags::kKinematics);
        } else {
            unset_flags(ShowerFlags::kKinematics);
        }
    }
    
    bool Shower::get_flag_kinematics(){
        return flags_any(ShowerFlags::kKinematics);
    }
    
    bool Shower::get_flag_shower(){
        return flags_any(ShowerFlags::kShower);
    }

    void Shower::add_shower(Shower& shower, const std::string& cloud_name_fit, const std::string& cloud_name_associate){
        // Iterate through all vertices in the input shower's view using the nodes() accessor
        for (auto vdesc : shower.nodes()) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (vtx && vtx->descriptor_valid()) {
                this->add_vertex(vtx);
            }
        }
        
        // Iterate through all segments in the input shower's view using the edges() accessor
        for (auto edesc : shower.edges()) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (seg && seg->descriptor_valid()) {
                this->add_segment(seg);
                
                // Merge point clouds from this segment
                // Handle "fit" point cloud
                if (!cloud_name_fit.empty()) {
                    auto seg_dpc_fit = seg->dpcloud(cloud_name_fit);
                    if (seg_dpc_fit) {
                        auto shower_dpc_fit = this->dpcloud(cloud_name_fit);
                        if (!shower_dpc_fit) {
                            this->dpcloud(cloud_name_fit, seg_dpc_fit);
                        } else {
                            shower_dpc_fit->add_points(seg_dpc_fit->get_points());
                        }
                    }
                }
                
                // Handle "associate_points" point cloud
                if (!cloud_name_associate.empty()) {
                    auto seg_dpc_associate = seg->dpcloud(cloud_name_associate);
                    if (seg_dpc_associate) {
                        auto shower_dpc_associate = this->dpcloud(cloud_name_associate);
                        if (!shower_dpc_associate) {
                            this->dpcloud(cloud_name_associate, seg_dpc_associate);
                        } else {
                            shower_dpc_associate->add_points(seg_dpc_associate->get_points());
                        }
                    }
                }
            }
        }
    }

    void Shower::complete_structure_with_start_segment(std::set<SegmentPtr>& used_segments, const std::string& cloud_name_fit, const std::string& cloud_name_associate) {
        if (!m_start_segment || !m_start_segment->descriptor_valid()) return;
        
        std::vector<SegmentPtr> new_segments;
        std::vector<VertexPtr> new_vertices;
        
        // Add start_segment to the view and mark as used
        used_segments.insert(m_start_segment);
        
      
        // Find vertices connected to start_segment (excluding start_vertex)
        auto vertices = find_vertices(m_full_graph, m_start_segment);
        if (vertices.first && vertices.first != m_start_vertex) {
            this->add_vertex(vertices.first);
            new_vertices.push_back(vertices.first);
        }
        if (vertices.second && vertices.second != m_start_vertex) {
            this->add_vertex(vertices.second);
            new_vertices.push_back(vertices.second);
        }
        
        // Worklist algorithm: explore connected segments and vertices
        while (!new_vertices.empty() || !new_segments.empty()) {
            // Process new vertices - find all segments connected to them
            if (!new_vertices.empty()) {
                VertexPtr vtx = new_vertices.back();
                new_vertices.pop_back();
                
                // Find all segments connected to this vertex
                if (vtx->descriptor_valid()) {
                    auto vdesc = vtx->get_descriptor();
                    for (auto edesc : boost::make_iterator_range(boost::out_edges(vdesc, m_full_graph))) {
                        SegmentPtr seg = m_full_graph[edesc].segment;
                        if (seg && seg->descriptor_valid() && used_segments.find(seg) == used_segments.end()) {
                            this->add_segment(seg);
                            new_segments.push_back(seg);
                            used_segments.insert(seg);
                            
                            // Merge point clouds from this segment
                            if (!cloud_name_fit.empty()) {
                                auto seg_dpc_fit = seg->dpcloud(cloud_name_fit);
                                if (seg_dpc_fit) {
                                    auto shower_dpc_fit = this->dpcloud(cloud_name_fit);
                                    if (!shower_dpc_fit) {
                                        this->dpcloud(cloud_name_fit, seg_dpc_fit);
                                    } else {
                                        shower_dpc_fit->add_points(seg_dpc_fit->get_points());
                                    }
                                }
                            }
                            
                            if (!cloud_name_associate.empty()) {
                                auto seg_dpc_associate = seg->dpcloud(cloud_name_associate);
                                if (seg_dpc_associate) {
                                    auto shower_dpc_associate = this->dpcloud(cloud_name_associate);
                                    if (!shower_dpc_associate) {
                                        this->dpcloud(cloud_name_associate, seg_dpc_associate);
                                    } else {
                                        shower_dpc_associate->add_points(seg_dpc_associate->get_points());
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            // Process new segments - find all vertices connected to them
            if (!new_segments.empty()) {
                SegmentPtr seg = new_segments.back();
                new_segments.pop_back();
                
                // Find vertices connected to this segment (excluding start_vertex)
                auto vertices = find_vertices(m_full_graph, seg);
                if (vertices.first && vertices.first != m_start_vertex) {
                    if (!this->has_node(vertices.first->get_descriptor())) {
                        this->add_vertex(vertices.first);
                        new_vertices.push_back(vertices.first);
                    }
                }
                if (vertices.second && vertices.second != m_start_vertex) {
                    if (!this->has_node(vertices.second->get_descriptor())) {
                        this->add_vertex(vertices.second);
                        new_vertices.push_back(vertices.second);
                    }
                }
            }
        }
    }

    void Shower::fill_sets(std::set<VertexPtr>& used_vertices, std::set<SegmentPtr>& used_segments, bool flag_exclude_start_segment){
        // Fill used_vertices with all vertices in this shower's view
        for (auto vdesc : this->nodes()) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (vtx) {
                used_vertices.insert(vtx);
            }
        }
        
        // Fill used_segments with all segments in this shower's view
        for (auto edesc : this->edges()) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (seg) {
                // Skip start_segment if flag is set
                if (flag_exclude_start_segment && seg == m_start_segment) {
                    continue;
                }
                used_segments.insert(seg);
            }
        }
    }

    void Shower::fill_point_vector(std::vector<WireCell::Point>& points, bool flag_main){
        // Get the main cluster ID if flag_main is true
        const Facade::Cluster* main_cluster = nullptr;
        if (flag_main && m_start_segment && m_start_segment->cluster()) {
            main_cluster = m_start_segment->cluster();
        }
        
        // Fill points from all segments in the shower's view
        for (auto edesc : this->edges()) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (seg) {
                // Skip if flag_main is set and segment is not in the main cluster
                if (flag_main && main_cluster && seg->cluster() != main_cluster) {
                    continue;
                }
                
                // Get segment fit points and add all except first and last
                const auto& fits = seg->fits();
                for (size_t i = 1; i + 1 < fits.size(); i++) {
                    points.push_back(fits[i].point);
                }
            }
        }
        
        // Fill points from all vertices in the shower's view
        for (auto vdesc : this->nodes()) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (vtx) {
                // Skip if flag_main is set and vertex is not in the main cluster
                if (flag_main && main_cluster && vtx->cluster() != main_cluster) {
                    continue;
                }
                
                // Add the vertex fit point
                points.push_back(vtx->fit().point);
            }
        }
    }

    TrajectoryView& Shower::fill_maps() {
        return *this;
    }

    std::pair<std::set<VertexPtr>, std::set<SegmentPtr>> Shower::get_connected_pieces(SegmentPtr seg){
        std::set<SegmentPtr> used_segments;
        std::set<VertexPtr> used_vertices;
        
        // Check if the segment is valid and in the view
        if (!seg || !seg->descriptor_valid() || !this->has_edge(seg->get_descriptor())) {
            return std::make_pair(used_vertices, used_segments);
        }
        
        std::vector<SegmentPtr> new_segments;
        std::vector<VertexPtr> new_vertices;
        
        // Start with the given segment
        new_segments.push_back(seg);
        used_segments.insert(seg);
        
        // Worklist algorithm: explore connected segments and vertices in the view
        while (!new_vertices.empty() || !new_segments.empty()) {
            // Process new vertices - find all segments connected to them in the view
            if (!new_vertices.empty()) {
                VertexPtr vtx = new_vertices.back();
                new_vertices.pop_back();
                
                // Find all segments connected to this vertex in the full graph, then check if in view
                if (vtx->descriptor_valid()) {
                    auto vdesc = vtx->get_descriptor();
                    for (auto edesc : boost::make_iterator_range(boost::out_edges(vdesc, m_full_graph))) {
                        // Check if this edge is in the view
                        if (this->has_edge(edesc)) {
                            SegmentPtr seg1 = m_full_graph[edesc].segment;
                            if (seg1 && used_segments.find(seg1) == used_segments.end()) {
                                new_segments.push_back(seg1);
                                used_segments.insert(seg1);
                            }
                        }
                    }
                }
            }
            
            // Process new segments - find all vertices connected to them in the view
            if (!new_segments.empty()) {
                SegmentPtr seg1 = new_segments.back();
                new_segments.pop_back();
                
                // Find vertices connected to this segment in the full graph
                auto vertices = find_vertices(m_full_graph, seg1);
                
                // Check if vertices are in the view and not yet visited
                if (vertices.first && this->has_node(vertices.first->get_descriptor()) 
                    && used_vertices.find(vertices.first) == used_vertices.end()) {
                    new_vertices.push_back(vertices.first);
                    used_vertices.insert(vertices.first);
                }
                if (vertices.second && this->has_node(vertices.second->get_descriptor()) 
                    && used_vertices.find(vertices.second) == used_vertices.end()) {
                    new_vertices.push_back(vertices.second);
                    used_vertices.insert(vertices.second);
                }
            }
        }
        
        return std::make_pair(used_vertices, used_segments);
    }

    std::pair<SegmentPtr, VertexPtr> Shower::get_last_segment_vertex_long_muon(std::set<SegmentPtr>& segments_in_muons) {
        VertexPtr s_vtx = m_start_vertex;
        SegmentPtr s_seg = m_start_segment;
        
        if (!s_vtx || !s_seg) {
            return std::make_pair(s_seg, s_vtx);
        }
        
        std::set<SegmentPtr> used_segments;
        used_segments.insert(s_seg);
        
        bool flag_continue = true;
        while (flag_continue) {
            flag_continue = false;
            
            // If current vertex is start_vertex, continue
            if (s_vtx == m_start_vertex) {
                flag_continue = true;
            } else {
                // Look for a new segment connected to s_vtx that is in segments_in_muons and not used
                if (s_vtx->descriptor_valid()) {
                    auto vdesc = s_vtx->get_descriptor();
                    // Iterate over segments connected to this vertex in the full graph, then check if in view
                    for (auto edesc : boost::make_iterator_range(boost::out_edges(vdesc, m_full_graph))) {
                        // Check if this edge is in the view
                        if (this->has_edge(edesc)) {
                            SegmentPtr sg = m_full_graph[edesc].segment;
                            if (sg && segments_in_muons.find(sg) != segments_in_muons.end() 
                                && used_segments.find(sg) == used_segments.end()) {
                                s_seg = sg;
                                used_segments.insert(s_seg);
                                flag_continue = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            // If we found a new segment, find the other vertex connected to it
            if (flag_continue) {
                auto vertices = find_vertices(m_full_graph, s_seg);
                if (vertices.first && vertices.first != s_vtx) {
                    s_vtx = vertices.first;
                } else if (vertices.second && vertices.second != s_vtx) {
                    s_vtx = vertices.second;
                }
            }
        }
        
        return std::make_pair(s_seg, s_vtx);
    }

    int Shower::get_num_main_segments() {
        int num = 0;
        
        // If no start segment, return 0
        if (!m_start_segment) {
            return 0;
        }
        
        // Get the start segment's cluster
        auto start_cluster = m_start_segment->cluster();
        if (!start_cluster) {
            return 0;
        }
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower
        for (auto edesc : this->edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;
            
            // Check if this segment's cluster matches the start segment's cluster
            if (seg->cluster() == start_cluster) {
                num++;
            }
        }
        
        return num;
    }

    int Shower::get_num_segments() {
        return this->edges().size();
    }

    double Shower::get_total_length(){
        double total_length = 0;
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower
        for (auto edesc : this->edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;
            
            // Add segment length
            total_length += segment_track_length(seg);
        }
        
        return total_length;
    }
    double Shower::get_total_length(Facade::Cluster* cluster){
        double total_length = 0;
        
        if (!cluster) {
            return 0;
        }
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower
        for (auto edesc : this->edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;
            
            // Check if segment's cluster matches the input cluster
            if (seg->cluster() == cluster) {
                total_length += segment_track_length(seg);
            }
        }
        
        return total_length;
    }
    double Shower::get_total_track_length(){
        double total_length = 0;
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower
        for (auto edesc : this->edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;
            
            // Only count segments that are NOT shower segments
            // Check if segment has shower flags (kShowerTrajectory or kShowerTopology)
            if (!seg->flags_any(SegmentFlags::kShowerTrajectory) && 
                !seg->flags_any(SegmentFlags::kShowerTopology)) {
                total_length += segment_track_length(seg);
            }
        }
        
        return total_length;
    }

    void Shower::update_particle_type(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
        double track_length = 0;
        double shower_length = 0;
        
        // Only process if there's more than one segment
        if (this->edges().size() <= 1) {
            return;
        }
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower
        for (auto edesc : this->edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;
            
            double length = segment_track_length(seg);
            
            // Check if segment is a shower segment OR not a proton (PDG 2212)
            bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) || 
                           seg->flags_any(SegmentFlags::kShowerTopology);
            
            bool is_not_proton = true;
            if (seg->has_particle_info()) {
                int pdg = seg->particle_info()->pdg();
                is_not_proton = (std::abs(pdg) != 2212);
            }
            
            if (is_shower || is_not_proton) {
                shower_length += length;
            } else {
                track_length += length;
            }
        }
        
        // If shower_length dominates, update start_segment to electron
        if (shower_length > track_length && m_start_segment) {
            // Calculate 4-momentum for electron (PDG = 11)
            auto four_momentum = segment_cal_4mom(m_start_segment, 11, particle_data, recomb_model);
            
            // Create ParticleInfo for electron
            auto pinfo = std::make_shared<Aux::ParticleInfo>(
                11,                                          // electron PDG
                particle_data->get_particle_mass(11),       // electron mass
                particle_data->pdg_to_name(11),             // "electron"
                four_momentum                                // 4-momentum
            );
            
            // Store particle info in start_segment
            m_start_segment->particle_info(pinfo);
        }
    }

    std::vector<double> Shower::get_stem_dQ_dx(VertexPtr vertex, SegmentPtr segment, int limit /*=20*/){
        std::vector<double> vec_dQ_dx;
        const double MIP_dQdx = 43e3 / units::cm;
        
        if (!vertex || !segment) {
            return vec_dQ_dx;
        }
        
        // Get dQ and dx from segment's fits
        const auto& fits = segment->fits();
        if (fits.empty()) {
            return vec_dQ_dx;
        }
        
        // Determine direction based on vertex position relative to segment
        // Check if vertex is at the front of the segment
        bool vertex_at_front = false;
        if (!segment->wcpts().empty()) {
            double d1 = WireCell::ray_length(WireCell::Ray{vertex->wcpt().point, segment->wcpts().front().point});
            double d2 = WireCell::ray_length(WireCell::Ray{vertex->wcpt().point, segment->wcpts().back().point});
            vertex_at_front = (d1 < d2);
        }
        
        // Fill vec_dQ_dx based on direction
        if (vertex_at_front) {
            for (size_t i = 0; i < fits.size(); i++) {
                double dQ_dx_normalized = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQdx;
                vec_dQ_dx.push_back(dQ_dx_normalized);
                if (vec_dQ_dx.size() >= (size_t)limit) break;
            }
        } else {
            for (int i = (int)fits.size() - 1; i >= 0; i--) {
                double dQ_dx_normalized = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQdx;
                vec_dQ_dx.push_back(dQ_dx_normalized);
                if (vec_dQ_dx.size() >= (size_t)limit) break;
            }
        }
        
        // If this is the start_segment and we don't have enough points, continue to next segments
        if (segment == m_start_segment && vec_dQ_dx.size() < (size_t)limit) {
            VertexPtr curr_vertex = vertex;
            SegmentPtr curr_segment = segment;
            int count = 0;
            
            while (vec_dQ_dx.size() < (size_t)limit && count < 3) {
                // Find next vertex (the other end of current segment)
                VertexPtr next_vertex = find_other_vertex(m_full_graph, curr_segment, curr_vertex);
                if (!next_vertex) break;
                
                // Direction from current vertex to next vertex
                WireCell::Vector dir1 = curr_vertex->fit().point - next_vertex->fit().point;
                
                // Find the next segment with largest angle
                SegmentPtr next_segment = nullptr;
                WireCell::Vector dir2;
                double max_angle = 0;
                
                auto next_vdesc = next_vertex->get_descriptor();
                for (auto edesc : boost::make_iterator_range(boost::out_edges(next_vdesc, m_full_graph))) {
                    if (!has_edge(edesc)) continue;  // Only consider edges in this view
                    
                    SegmentPtr seg = m_full_graph[edesc].segment;
                    if (seg == curr_segment) continue;
                    
                    WireCell::Vector tmp_dir = segment_cal_dir_3vector(seg, next_vertex->fit().point, 10 * units::cm);
                    double angle = std::acos(std::clamp(dir1.dot(tmp_dir) / (dir1.magnitude() * tmp_dir.magnitude() + 1e-9), -1.0, 1.0)) * 180.0 / M_PI;
                    
                    if (angle > max_angle) {
                        max_angle = angle;
                        next_segment = seg;
                        dir2 = tmp_dir;
                    }
                }
                
                if (!next_segment) break;
                
                // Check if there are other segments that would make this "bad"
                bool flag_bad = false;
                for (auto edesc : boost::make_iterator_range(boost::out_edges(next_vdesc, m_full_graph))) {
                    if (!has_edge(edesc)) continue;
                    
                    SegmentPtr seg = m_full_graph[edesc].segment;
                    if (seg == curr_segment || seg == next_segment) continue;
                    
                    double seg_length = segment_track_length(seg);
                    if (seg_length > 3 * units::cm) {
                        WireCell::Vector tmp_dir = segment_cal_dir_3vector(seg, next_vertex->fit().point, 10 * units::cm);
                        double angle = std::acos(std::clamp(dir2.dot(tmp_dir) / (dir2.magnitude() * tmp_dir.magnitude() + 1e-9), -1.0, 1.0)) * 180.0 / M_PI;
                        if (angle < 25) {
                            flag_bad = true;
                            break;
                        }
                    }
                }
                
                if (flag_bad) break;
                
                // Remove last element and add points from next segment
                if (!vec_dQ_dx.empty()) {
                    vec_dQ_dx.pop_back();
                }
                
                const auto& next_fits = next_segment->fits();
                if (next_fits.empty()) break;
                
                // Determine direction for next segment
                bool next_vertex_at_front = false;
                if (!next_segment->wcpts().empty()) {
                    double d1 = WireCell::ray_length(WireCell::Ray{next_vertex->wcpt().point, next_segment->wcpts().front().point});
                    double d2 = WireCell::ray_length(WireCell::Ray{next_vertex->wcpt().point, next_segment->wcpts().back().point});
                    next_vertex_at_front = (d1 < d2);
                }
                
                // Add dQ/dx from next segment
                if (next_vertex_at_front) {
                    for (size_t i = 0; i < next_fits.size(); i++) {
                        double dQ_dx_normalized = next_fits[i].dQ / (next_fits[i].dx + 1e-9) / MIP_dQdx;
                        vec_dQ_dx.push_back(dQ_dx_normalized);
                        if (vec_dQ_dx.size() >= (size_t)limit) break;
                    }
                } else {
                    for (int i = (int)next_fits.size() - 1; i >= 0; i--) {
                        double dQ_dx_normalized = next_fits[i].dQ / (next_fits[i].dx + 1e-9) / MIP_dQdx;
                        vec_dQ_dx.push_back(dQ_dx_normalized);
                        if (vec_dQ_dx.size() >= (size_t)limit) break;
                    }
                }
                
                if (vec_dQ_dx.size() >= (size_t)limit) break;
                
                // Prepare for next iteration
                curr_vertex = next_vertex;
                curr_segment = next_segment;
                count++;
            }
        }
        
        return vec_dQ_dx;
    }

    void Shower::calculate_kinematics(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
        int nsegments = this->edges().size();
        
        if (nsegments == 1) {
            // Single segment case
            if (!m_start_segment) return;
            
            // Set particle type from start segment
            if (m_start_segment->has_particle_info()) {
                data.particle_type = m_start_segment->particle_info()->pdg();
            }
            
            // Check if shower
            bool flag_shower = m_start_segment->flags_any(SegmentFlags::kShowerTrajectory) || 
                             m_start_segment->flags_any(SegmentFlags::kShowerTopology);
            
            // Calculate energies
            double seg_length = segment_track_length(m_start_segment);
            data.kenergy_range = cal_kine_range(seg_length, data.particle_type, particle_data);
            data.kenergy_dQdx = segment_cal_kine_dQdx(m_start_segment, recomb_model);
            
            // Calculate kenergy_best
            if (data.start_connection_type == 1) {
                data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
            } else {
                if (flag_shower) {
                    data.kenergy_best = 0;
                } else {
                    data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
                }
            }
            
            // Calculate start_point and end_point
            const auto& fits = m_start_segment->fits();
            if (data.start_connection_type == 1 || !this->dpcloud("fit")) {
                if (!fits.empty()) {
                    if (m_start_segment->dirsign() == 1) {
                        data.start_point = fits.front().point;
                        data.end_point = fits.back().point;
                    } else if (m_start_segment->dirsign() == -1) {
                        data.start_point = fits.back().point;
                        data.end_point = fits.front().point;
                    }
                }
            } else {
                if (m_start_vertex) {
                    data.start_point = shower_get_closest_point(*this, m_start_vertex->fit().point, "fit").second;
                    
                    // Find farthest vertex
                    double max_dis = 0;
                    const auto& view = this->view_graph();
                    for (auto vdesc : this->nodes()) {
                        VertexPtr vtx = view[vdesc].vertex;
                        if (!vtx) continue;
                        double dis = (data.start_point - vtx->fit().point).magnitude();
                        if (dis > max_dis) {
                            max_dis = dis;
                            data.end_point = vtx->fit().point;
                        }
                    }
                }
            }
            
            // Calculate init_dir
            if (data.start_connection_type == 1) {
                data.init_dir = segment_cal_dir_3vector(m_start_segment);
            } else if (data.start_connection_type == 2 || data.start_connection_type == 3) {
                if (m_start_vertex) {
                    data.init_dir = (data.start_point - m_start_vertex->fit().point).norm();
                }
            }
            
        } else {
            // Multiple segments case
            if (!m_start_segment) return;
            
            // Get number of connected segments
            auto [segs, verts] = get_connected_pieces(m_start_segment);
            int nconnected_segs = segs.size();
            
            // Set particle type
            if (m_start_segment->has_particle_info()) {
                data.particle_type = m_start_segment->particle_info()->pdg();
            }
            
            bool flag_shower = m_start_segment->flags_any(SegmentFlags::kShowerTrajectory) || 
                             m_start_segment->flags_any(SegmentFlags::kShowerTopology);
            
            if (nsegments == nconnected_segs) {
                // Single track (all connected)
                
                // Calculate start_point
                const auto& fits = m_start_segment->fits();
                if (data.start_connection_type == 1 || !this->dpcloud("fit")) {
                    if (!fits.empty()) {
                        if (m_start_segment->dirsign() == 1) {
                            data.start_point = fits.front().point;
                        } else if (m_start_segment->dirsign() == -1) {
                            data.start_point = fits.back().point;
                        }
                    }
                } else {
                    if (m_start_vertex) {
                        data.start_point = shower_get_closest_point(*this, m_start_vertex->fit().point, "fit").second;
                    }
                }
                
                // Calculate init_dir
                double seg_length = segment_track_length(m_start_segment);
                if (data.start_connection_type == 1) {
                    if (seg_length > 8 * units::cm) {
                        data.init_dir = segment_cal_dir_3vector(m_start_segment);
                    } else if (m_start_vertex) {
                        data.init_dir = shower_cal_dir_3vector(*this, m_start_vertex->fit().point, 12 * units::cm);
                    }
                } else if (data.start_connection_type == 2 || data.start_connection_type == 3) {
                    if (m_start_vertex) {
                        data.init_dir = (data.start_point - m_start_vertex->fit().point).norm();
                    }
                }
                
                // Find farthest vertex for end_point
                double max_dis = 0;
                const auto& view = this->view_graph();
                for (auto vdesc : this->nodes()) {
                    VertexPtr vtx = view[vdesc].vertex;
                    if (!vtx) continue;
                    double dis = (data.start_point - vtx->fit().point).magnitude();
                    if (dis > max_dis) {
                        max_dis = dis;
                        data.end_point = vtx->fit().point;
                    }
                }
                
                // Collect all dQ and dx from all segments
                double total_length = 0;
                std::vector<double> vec_dQ, vec_dx;
                for (auto edesc : this->edges()) {
                    SegmentPtr seg = view[edesc].segment;
                    if (!seg) continue;
                    
                    total_length += segment_track_length(seg);
                    
                    const auto& seg_fits = seg->fits();
                    for (const auto& fit : seg_fits) {
                        vec_dQ.push_back(fit.dQ);
                        vec_dx.push_back(fit.dx);
                    }
                }
                
                // Calculate energies
                data.kenergy_range = cal_kine_range(total_length, data.particle_type, particle_data);
                data.kenergy_dQdx = cal_kine_dQdx(vec_dQ, vec_dx, recomb_model);
                
                // Calculate kenergy_best
                if (data.start_connection_type == 1) {
                    data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
                } else {
                    if (flag_shower) {
                        data.kenergy_best = 0;
                    } else {
                        data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
                    }
                }
                
            } else {
                // Multiple tracks (not all connected)
                
                // Calculate start_point
                const auto& fits = m_start_segment->fits();
                if (data.start_connection_type == 1 || !this->dpcloud("fit")) {
                    if (!fits.empty()) {
                        if (m_start_segment->dirsign() == 1) {
                            data.start_point = fits.front().point;
                        } else if (m_start_segment->dirsign() == -1) {
                            data.start_point = fits.back().point;
                        }
                    }
                } else {
                    if (m_start_vertex) {
                        data.start_point = shower_get_closest_point(*this, m_start_vertex->fit().point, "fit").second;
                    }
                }
                
                // Calculate init_dir
                double seg_length = segment_track_length(m_start_segment);
                if (data.start_connection_type == 1) {
                    if (seg_length > 8 * units::cm) {
                        data.init_dir = segment_cal_dir_3vector(m_start_segment);
                    } else if (m_start_vertex) {
                        data.init_dir = shower_cal_dir_3vector(*this, m_start_vertex->fit().point, 12 * units::cm);
                    }
                } else if (data.start_connection_type == 2 || data.start_connection_type == 3) {
                    if (m_start_vertex) {
                        data.init_dir = (data.start_point - m_start_vertex->fit().point).norm();
                    }
                }
                
                // Find farthest vertex for end_point
                double max_dis = 0;
                const auto& view = this->view_graph();
                for (auto vdesc : this->nodes()) {
                    VertexPtr vtx = view[vdesc].vertex;
                    if (!vtx) continue;
                    double dis = (data.start_point - vtx->fit().point).magnitude();
                    if (dis > max_dis) {
                        max_dis = dis;
                        data.end_point = vtx->fit().point;
                    }
                }
                
                // Collect all dQ and dx from all segments
                std::vector<double> vec_dQ, vec_dx;
                for (auto edesc : this->edges()) {
                    SegmentPtr seg = view[edesc].segment;
                    if (!seg) continue;
                    
                    const auto& seg_fits = seg->fits();
                    for (const auto& fit : seg_fits) {
                        vec_dQ.push_back(fit.dQ);
                        vec_dx.push_back(fit.dx);
                    }
                }
                
                // Calculate energies
                data.kenergy_range = 0;
                data.kenergy_dQdx = cal_kine_dQdx(vec_dQ, vec_dx, recomb_model);
                data.kenergy_best = 0;
            }
        }
    }

    void Shower::calculate_kinematics_long_muon(std::set<SegmentPtr>& segments_in_muons, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
        // Get particle type from start segment
        int particle_type = abs(m_start_segment->particle_info()->pdg());
        double particle_mass = m_start_segment->particle_info()->mass();
        
        unset_flags(ShowerFlags::kKinematics);
        
        // Calculate total length ONLY from segments in segments_in_muons
        double total_length = 0;
        for (auto edesc : this->edges()) {
            if (!has_edge(edesc)) continue;
            auto seg = view_graph()[edesc].segment;
            if (segments_in_muons.find(seg) != segments_in_muons.end()) {
                total_length += segment_track_length(seg);
            }
        }
        
        // Collect dQ and dx from ALL segments
        std::vector<double> vec_dQ;
        std::vector<double> vec_dx;
        
        for (auto edesc : this->edges()) {
            if (!has_edge(edesc)) continue;
            auto seg = view_graph()[edesc].segment;
            
            std::vector<double> seg_dQ;
            std::vector<double> seg_dx;
            
             const auto& seg_fits = seg->fits();
            for (const auto& fit : seg_fits) {
                vec_dQ.push_back(fit.dQ);
                vec_dx.push_back(fit.dx);
            }
            
            vec_dQ.insert(vec_dQ.end(), seg_dQ.begin(), seg_dQ.end());
            vec_dx.insert(vec_dx.end(), seg_dx.begin(), seg_dx.end());
        }
        
        // Calculate kinetic energies
        data.kenergy_range = cal_kine_range(total_length, particle_type, particle_data);
        data.kenergy_dQdx = cal_kine_dQdx(vec_dQ, vec_dx, recomb_model);
        
        // For long muon, use dQdx as best energy
        data.kenergy_best = data.kenergy_dQdx;
        
        // Calculate initial direction from start segment
        data.init_dir = segment_cal_dir_3vector(m_start_segment);
        
        // Set start point based on direction
        auto& fits = m_start_segment->fits();
        int dirsign_val = m_start_segment->dirsign();
        if (dirsign_val == 1) {
            data.start_point = fits.front().point;
        } else {
            data.start_point = fits.back().point;
        }
        
        // Find farthest vertex that has at least one segment in segments_in_muons
        double max_dis = 0;
        VertexPtr farthest_vertex = nullptr;
        
        for (auto vdesc : this->nodes()) {
            if (!has_node(vdesc)) continue;
            auto vtx = view_graph()[vdesc].vertex;
            
            // Check if this vertex has at least one segment in segments_in_muons
            bool flag_contain = false;
            for (auto out_edge : boost::make_iterator_range(boost::out_edges(vdesc, m_full_graph))) {
                if (!has_edge(out_edge)) continue;
                auto seg = view_graph()[out_edge].segment;
                if (segments_in_muons.find(seg) != segments_in_muons.end()) {
                    flag_contain = true;
                    break;
                }
            }
            
            if (flag_contain) {
                double dis = (vtx->fit().point - data.start_point).magnitude();
                if (dis > max_dis) {
                    max_dis = dis;
                    farthest_vertex = vtx;
                }
            }
        }
        
        // Set end point to the farthest vertex
        if (farthest_vertex) {
            data.end_point = farthest_vertex->fit().point;
        } else {
            // Fallback: use the other end of start segment
            auto& fits = m_start_segment->fits();
            if (m_start_segment->dirsign() == 1) {
                data.end_point = fits.back().point;
            } else {
                data.end_point = fits.front().point;
            }
        }
    }

}