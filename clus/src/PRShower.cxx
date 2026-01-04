#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRGraph.h"
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


}