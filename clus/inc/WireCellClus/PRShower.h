#ifndef WIRECELL_CLUS_PR_SHOWER
#define WIRECELL_CLUS_PR_SHOWER

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRTrajectoryView.h"

#include "WireCellUtil/Flagged.h"
#include "WireCellUtil/Point.h"

#include "WireCellIface/IRecombinationModel.h"
#include "WireCellClus/ParticleDataSet.h"

namespace WireCell::Clus::PR {

    /** The "flags" that may be set on a shower.
     */
    enum class ShowerFlags {
        /// The segment has no particular category.
        kUndefined = 0,
        /// The shower flag
        kShower = 1<<1,
        /// The kinematics flag
        kKinematics = 1<<2,
    };

    /** The data attributes of a PR::Shower
        
        The original WCShower has a huge number of attributes that are merely
        carried by the shower.  This struct factors out that data to facilitate
        writing stand-alone functions.

        Anything that is not part of the core shower-as-graph-view concept gets
        lumped into this ShowerData struct.

        Note, "flags" are set via Shower's Flagged base class.

    */
    struct ShowerData
    {
        int particle_type;
        double kenergy_range;
        double kenergy_dQdx;
        double kenergy_charge;
        double kenergy_best;
    
        WireCell::Point start_point;
        WireCell::Point end_point;
        WireCell::Vector init_dir;

        /// 1 for direct connection, 2 for indirection connection with a gap, 3
        /// for associations (not clear if this should be connected or not
        int start_connection_type;
    };


    /** Model a shower-like view of a trajectory.

        This is the WCT equivalent to a WCT WCShower.
     */
    class Shower
        : public TrajectoryView
        , public Flagged<ShowerFlags> // can set flags
        , public HasDPCs<Shower>      // has associated DynamicPointClouds.
    {
    public:

        Shower(Graph& graph);
        
        virtual ~Shower();

        // The bag of attributes is directly exposed to user.
        ShowerData data;


        // Getters

        /** Get the vertex that is considered the start of the shower.
         */
        VertexPtr start_vertex();

        /** Get the segment that is considered the start of the shower.
         */
        SegmentPtr start_segment();

        // Chainable setters

        /** Chainable setter of start vertex.

            The vertex must already be added to the underlying graph that this
            shower views.

            The vertex will be added to the view.

            The vertex will replace any existing start vertex and not remove the
            prior vertex from the shower's view.  Use `Shower::remove_vertex()`
            to explicitly remove from the view and `PR::remove_vertex()` to
            remove it from the underlying graph, if either are required.

            If the vertex lacks a valid descriptor, eg has yet to be added to
            the underlying graph, this function is a no-op and the stored
            start_vertex is nullified.
        */
        Shower& set_start_vertex(VertexPtr vtx, int type);
        std::pair<VertexPtr, int> get_start_vertex_and_type() {
            return std::make_pair(m_start_vertex, data.start_connection_type);
        }

        /** Chainable setter of start segment.

            This has the same semantics and caveats as the chainable setter:
            `start_vertex(VertexPtr)`.
        */
        Shower& set_start_segment(SegmentPtr seg, bool flag_include_vertices = false, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");

        Shower& set_start_point( WireCell::Point pt ) {
            data.start_point = pt;
            return *this;
        }
        WireCell::Point get_start_point() const {
            return data.start_point;
        }
        WireCell::Point get_end_point() const {
            return data.end_point;
        }

        void add_segment(SegmentPtr seg, bool flag_include_vertices = false, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");

        // particle type
        int get_particle_type(){return data.particle_type;};
        void set_particle_type(int val){data.particle_type = val;};

        // access flags
        void set_flag_kinematics(bool val);
        bool get_flag_kinematics();
        bool get_flag_shower();

        // return kinematic energy estimation
        double get_kine_range(){return data.kenergy_range;};
        double get_kine_dQdx(){return data.kenergy_dQdx;};
        void set_kine_charge(double val){data.kenergy_charge = val;};
        double get_kine_charge(){return data.kenergy_charge;};
        double get_kine_best(){return data.kenergy_best;};

        // return initial direction ...
        WireCell::Vector& get_init_dir(){return data.init_dir;};

        // Get point cloud by name (convenience wrapper around dpcloud)
        std::shared_ptr<Facade::DynamicPointCloud> get_pcloud(const std::string& cloud_name = "fit") {
            return this->dpcloud(cloud_name);
        }
        std::shared_ptr<const Facade::DynamicPointCloud> get_pcloud(const std::string& cloud_name = "associate_points") const {
            return this->dpcloud(cloud_name);
        }

        // Add all segments and vertices from another shower to this one
        void add_shower(Shower& temp_shower, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");
        void complete_structure_with_start_segment(std::set<SegmentPtr>& used_segments, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");


        // get the information from the shower
        void fill_sets(std::set<VertexPtr>& used_vertices, std::set<SegmentPtr>& used_segments, bool flag_exclude_start_segment = true);
        void fill_point_vector(std::vector<WireCell::Point>& points, bool flag_main = true);
        TrajectoryView& fill_maps();

        std::pair<std::set<VertexPtr>, std::set<SegmentPtr>> get_connected_pieces(SegmentPtr seg);

        std::pair<SegmentPtr, VertexPtr> get_last_segment_vertex_long_muon(std::set<SegmentPtr>& segments_in_muons);

        // some simple get functions
        int get_num_main_segments();
        int get_num_segments();

        double get_total_length();
        double get_total_length(Facade::Cluster* cluster);
        double get_total_track_length();

        std::vector<double> get_stem_dQ_dx(VertexPtr vertex, SegmentPtr segment, int limit = 20);

        // calculate the kinematics
        void update_particle_type(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void calculate_kinematics(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void calculate_kinematics_long_muon(std::set<SegmentPtr>& segments_in_muons, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

    private:

        Graph& m_full_graph;
        VertexPtr m_start_vertex;
        SegmentPtr m_start_segment;

    };

    using ShowerPtr = std::shared_ptr<Shower>;

}
#endif
