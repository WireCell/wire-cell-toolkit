#ifndef WIRECELLCLUS_TRACKFITTING_H
#define WIRECELLCLUS_TRACKFITTING_H

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Logging.h"
#include "WireCellClus/PRGraph.h"

#include <Eigen/IterativeLinearSolvers>


namespace WireCell::Clus {

    /**
     * Dedicated TrackFitting class that can be instantiated and used by 
     * other ensemble visitors without needing to be configured as a component.
     * 
     * This class encapsulates track fitting algorithms that can work on
     * individual clusters or collections of clusters.
     */
    class TrackFitting {
    public:
        
        enum class FittingType {
            Single,
            Multiple
        };

        /**
         * Structure to hold all track fitting parameters in one place
         */
        struct Parameters {
            // Diffusion coefficients (LArTPC standard values)
            double DL = 6.4* pow(units::cm,2)/units::second;                    // m²/s, longitudinal diffusion
            double DT = 9.8* pow(units::cm,2)/units::second;                    // m²/s, transverse diffusion

            // Software filter effects (wire dimension broadening)
            double col_sigma_w_T = 0.188060 * 3*units::mm * 0.2; // Collection plane, units: wire pitch
            double ind_sigma_u_T = 0.402993 * 3*units::mm * 0.3; // U induction plane
            double ind_sigma_v_T = 0.402993 * 3*units::mm * 0.5; // V induction plane

            // Uncertainty parameters
            double rel_uncer_ind = 0.075;          // Relative uncertainty for induction planes
            double rel_uncer_col = 0.05;           // Relative uncertainty for collection plane
            double add_uncer_ind = 0.0;            // Additional uncertainty for induction
            double add_uncer_col = 300.0;          // Additional uncertainty for collection
            
            // Longitudinal filter effects (time dimension)
            double add_sigma_L = 1.428249 *0.5505*units::mm /0.5;   
             
            // Additional useful parameters for charge err estimation ...
            double rel_charge_uncer = 0.1; // 10% 
            double add_charge_uncer = 600; // electrons

            double default_charge_th = 100;
            double default_charge_err = 1000; 

            double scaling_quality_th = 0.5;
            double scaling_ratio = 0.05;

            double area_ratio1 = 1.8*units::mm;
            double area_ratio2 = 1.7;

            double skip_default_ratio_1 = 0.25;
            double skip_ratio_cut = 0.97;
            double skip_ratio_1_cut = 0.75;

            double skip_angle_cut_1 = 160;
            double skip_angle_cut_2 = 90;
            double skip_angle_cut_3 = 45;
            double skip_dis_cut = 0.5*units::cm;

            double default_dQ_dx = 5000; 

            double end_point_factor=0.6;
            double mid_point_factor=0.9;
            int nlevel=3;
            double charge_cut=2000;
            
            double low_dis_limit = 1.2*units::cm;            // cm, lower distance limit for point organization
            double end_point_limit = 0.6*units::cm;          // cm, extension distance for end points
            double time_tick_cut = 20;   //            //  tick cut for point association

            // addition parameters
            double share_charge_err = 8000;
            double min_drift_time = 50*units::us;
            double search_range = 10; // wires, or time slices (not ticks)

            double dead_ind_weight = 0.3;
            double dead_col_weight = 0.9;
            double close_ind_weight = 0.15;
            double close_col_weight = 0.45;
            double overlap_th = 0.5;
            double dx_norm_length = 0.6*units::cm;
            double lambda= 0.0005;

            double div_sigma = 0.6*units::cm;            
        };
 
        /**
         * Constructor
         * @param fitting_type The type of fitting to perform (single or multiple tracks)
         */
        explicit TrackFitting(FittingType fitting_type = FittingType::Single);   
        virtual ~TrackFitting() = default;

        /**
         * Set the fitting type
         * @param fitting_type The new fitting type to use
         */
        void set_fitting_type(FittingType fitting_type) { m_fitting_type = fitting_type; }

        /**
         * Get the current fitting type
         * @return The current fitting type
         */
        FittingType get_fitting_type() const { return m_fitting_type; }

        // Parameter management methods
        
        /**
         * Get read-only access to current parameters
         */
        const Parameters& get_parameters() const { return m_params; }
        
        /**
         * Set new parameters (replaces all current parameters)
         */
        void set_parameters(const Parameters& params) { m_params = params; }
        
        /**
         * Set specific parameter by name
         */
        void set_parameter(const std::string& name, double value);
        
        /**
         * Get specific parameter by name
         */
        double get_parameter(const std::string& name) const;

        // single track fitting utilizes the segments ... 
        void add_segment(std::shared_ptr<PR::Segment> segment);
        /**
         * Get the set of segments currently stored in this TrackFitting instance.
         * @return Set of shared pointers to PR::Segment
         */
        std::set<std::shared_ptr<PR::Segment>> get_segments() const { return m_segments; }
        void clear_segments();
 
        // multi-track fitting utilized the Graph ... 
        void add_graph(std::shared_ptr<PR::Graph> graph);
        std::shared_ptr<PR::Graph> get_graph() const { return m_graph; }
        void clear_graph();



        // collect charge
        void prepare_data();

        // Fill the global readout map
        void fill_global_rb_map();

        /**
         * Organize original path from segment points with distance limits
         * @param segment Pointer to PR::Segment containing the path points
         * @param low_dis_limit Lower distance limit for point organization
         * @param end_point_limit Extension distance for end points
         * @return Vector of organized 3D points
         */
        std::vector<WireCell::Point> organize_orig_path(std::shared_ptr<PR::Segment> segment, double low_dis_limit=1.2*units::cm, double end_point_limit=0.6*units::cm);

        std::vector<WireCell::Point> examine_end_ps_vec(std::shared_ptr<PR::Segment> segment, const std::vector<WireCell::Point>& pts, bool flag_start, bool flag_end);

        void organize_ps_path(std::shared_ptr<PR::Segment> segment, std::vector<WireCell::Point>& pts, double low_dis_limit, double end_point_limit);


    

                /// Internal coordinate (can be more complex)
        struct Coord2D {
            int apa, face, time, wire, channel;
            WirePlaneLayer_t plane;  // Additional internal information

            Coord2D(int a, int f, int t, int w, int c, WirePlaneLayer_t p)
                : apa(a), face(f), time(t), wire(w), channel(c), plane(p) {}

            bool operator<(const Coord2D& other) const {
                if (apa != other.apa) return apa < other.apa;
                if (face != other.face) return face < other.face;
                if (time != other.time) return time < other.time;
                if (wire != other.wire) return wire < other.wire;
                if (channel != other.channel) return channel < other.channel;
                return plane < other.plane;
            }
        };

        /// Per-plane data for 3D points (exactly matches prototype)
        struct PlaneData {
            std::set<Coord2D> associated_2d_points;
            double quantity;
            
            PlaneData() : quantity(0.0) {}
        };

        /// 3D point with per-plane associations (corrected structure)
        struct Point3DInfo {
            std::map<WirePlaneLayer_t, PlaneData> plane_data;
            
            const PlaneData& get_plane_data(WirePlaneLayer_t plane) const {
                static PlaneData empty;
                auto it = plane_data.find(plane);
                return (it != plane_data.end()) ? it->second : empty;
            }
            
            void set_plane_data(WirePlaneLayer_t plane, const PlaneData& data) {
                plane_data[plane] = data;
            }
        };

        struct CoordReadout {
            int apa, time, channel;

            CoordReadout(int a, int t, int c)
            : apa(a), time(t), channel(c) {}

            bool operator<(const CoordReadout& other) const {
            if (apa != other.apa) return apa < other.apa;
            if (time != other.time) return time < other.time;
            return channel < other.channel;
            }
        };


        /// Simple charge measurement (in ternal interface)
        struct ChargeMeasurement {
            double charge, charge_err;
            int flag;
            
            ChargeMeasurement(double q = 0.0, double qe = 0.0, int f = 0) 
                : charge(q), charge_err(qe), flag(f) {}
        };



        // point associations
        void form_point_association(std::shared_ptr<PR::Segment> segment, WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, double dis_cut, int nlevel, double time_tick_cut );

        void examine_point_association(std::shared_ptr<PR::Segment> segment, WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, bool flag_end_point = false, double charge_cut = 2000);

        void form_map(std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& ptss, double end_point_factor=0.6, double mid_point_factor=0.9, int nlevel=3, double time_tick_cut=20, double charge_cut=2000);

        // track trajectory fitting // should fit all APA ...
        void trajectory_fit(std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec, int charge_div_method = 1, double div_sigma = 0.6*units::cm);

        bool skip_trajectory_point(WireCell::Point& p, std::pair<int,int>& apa_face, int i, std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec,  std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& fine_tracking_path); 

        // prepare for dQ/dx fitting
        double cal_gaus_integral(int tbin, int wbin, double t_center, double t_sigma, 
                                       double w_center, double w_sigma, int flag, double nsigma, int cur_ntime_ticks);

        double cal_gaus_integral_seg(int tbin, int wbin, std::vector<double>& t_centers, std::vector<double>& t_sigmas, std::vector<double>& w_centers, std::vector<double>& w_sigmas, std::vector<double>& weights, int flag, double nsigma, int cur_ntime_ticks);

        void update_dQ_dx_data();
        void recover_original_charge_data();

        /**
         * Calculate compact matrix analysis for wire plane sharing
         * 
         * This function analyzes the sharing patterns between 2D measurements and 3D positions
         * to compute overlap ratios and adjust weight matrix coefficients. It processes sparse
         * matrices representing the relationship between 2D wire measurements and 3D positions.
         * 
         * @param weight_matrix Reference to sparse weight matrix (MW, MV, or MU) to be modified
         * @param response_matrix_transpose Transposed response matrix (RWT, RVT, or RUT)  
         * @param n_2d_measurements Number of 2D measurements (wire/time points)
         * @param n_3d_positions Number of 3D positions
         * @param cut_position Threshold for wire sharing cut (default 2.0)
         * @return Vector of pairs containing overlap ratios for each 3D position
         *         Each pair contains (previous_neighbor_ratio, next_neighbor_ratio)
         */
        std::vector<std::pair<double, double>> calculate_compact_matrix(Eigen::SparseMatrix<double>& weight_matrix, const Eigen::SparseMatrix<double>& response_matrix_transpose, int n_2d_measurements, int n_3d_positions, double cut_position = 2.0);

        void dQ_dx_fill(double dis_end_point_ext=0.45*units::cm);

        void dQ_dx_fit(double dis_end_point_ext=0.45*units::cm, bool flag_dQ_dx_fit_reg=true);

        void do_single_tracking(std::shared_ptr<PR::Segment> segment, bool flag_dQ_dx_fit_reg= true, bool flag_dQ_dx_fit= true, bool flag_force_load_data = false, bool flag_hack = false);

        

        /**  
         * Get anode for a specific APA identifier
         * @param apa_ident APA identifier (typically same as APA number)
         * @return Pointer to IAnodePlane, or nullptr if not found
         */
        IAnodePlane::pointer get_anode(int apa_ident = 0) const;

        /**
         * Get all available anodes from the grouping
         * @return Map of APA identifier to anode pointer
         */
        std::map<int, IAnodePlane::pointer> get_all_anodes() const;

        /**
         * Get channel number for a specific wire location
         * Uses hybrid caching for optimal performance
         * @param apa APA number
         * @param face Face number (0 or 1)
         * @param plane Plane index (0=U, 1=V, 2=W typically)  
         * @param wire Wire index within the plane
         * @return Channel number, or -1 if invalid
         */
        int get_channel_for_wire(int apa, int face, int plane, int wire) const;

        /**
         * Get all wires that belong to a specific channel
         * @param apa APA number
         * @param channel_number Channel identifier
         * @return Vector of wire information (face, plane, wire_index)
         */
        std::vector<std::tuple<int, int, int>> get_wires_for_channel(int apa, int channel_number) const;

        /**
         * Clear all caches (useful for memory management)
         */
        void clear_cache() const;

        /**
         * Get cache statistics for monitoring/debugging
         */
        struct CacheStats {
            size_t hot_planes_count;
            size_t cold_entries_count;
            size_t total_lookups;
            size_t hot_hits;
            size_t cold_hits;
            double hit_rate() const { 
                return total_lookups > 0 ? (double)(hot_hits + cold_hits) / total_lookups : 0.0; 
            }
        };
        CacheStats get_cache_stats() const;

        /**
         * Set the detector volume for this TrackFitting instance
         * @param dv Pointer to IDetectorVolumes
         */
        void set_detector_volume(IDetectorVolumes::pointer dv) { m_dv = dv; }
        
        /**
         * Set the PCTransformSet for coordinate transformations
         * @param pcts Pointer to PCTransformSet interface
         */
        void set_pc_transforms(IPCTransformSet::pointer pcts) { m_pcts = pcts; }

        /**
         * Get the current detector volumes
         * @return Pointer to detector volumes interface
         */
        IDetectorVolumes::pointer get_detector_volume() const { return m_dv; }

        /**
         * Get the current PCTransformSet
         * @return Pointer to PCTransformSet interface
         */
        IPCTransformSet::pointer get_pc_transforms() const { return m_pcts; }

        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> get_fine_tracking_path() const { return fine_tracking_path; }
        std::vector<double> get_dQ() const { return dQ; }
        std::vector<double> get_dx() const { return dx; }
        std::vector<double> get_pu() const { return pu; }
        std::vector<double> get_pv() const { return pv; }
        std::vector<double> get_pw() const { return pw; }
        std::vector<double> get_pt() const { return pt; }
        std::vector<std::pair<int,int>> get_paf() const {return paf;}
        std::vector<double> get_reduced_chi2() const { return reduced_chi2; }

    private:
         // Core parameters - centralized storage
        Parameters m_params;

        // Helper method to get parameter value or default
        double get_param_or_default(double param_value, double default_value) const {
            return (param_value < 0) ? default_value : param_value;
        }

        FittingType m_fitting_type;
        IDetectorVolumes::pointer m_dv{nullptr};  
        IPCTransformSet::pointer m_pcts{nullptr};          // PC Transform Set
        
        // cluster and grouping, CTPC is from m_grouping ...
        Facade::Grouping* m_grouping{nullptr};
        std::set<Facade::Cluster*> m_clusters;

        std::set<Facade::Blob*> m_blobs;

        // input segment
        std::set<std::shared_ptr<PR::Segment> > m_segments;

        // input graph 
        std::shared_ptr<PR::Graph> m_graph{nullptr};

        // =====================================================================
        // HYBRID CACHE IMPLEMENTATION
        // =====================================================================
        
        // Key types for caching
        using PlaneKey = std::tuple<int, int, int>;    // (apa, face, plane)
        using WireKey = std::tuple<int, int, int, int>; // (apa, face, plane, wire)
        
        // Hot cache: frequently accessed plane mappings (full plane cached)
        mutable std::map<PlaneKey, std::vector<int>> m_hot_cache;
        
        // Cold cache: individual wire lookups
        mutable std::map<WireKey, int> m_cold_cache;
        
        // Access frequency tracking
        mutable std::map<PlaneKey, int> m_access_count;
        
        // Cache statistics
        mutable CacheStats m_cache_stats = {0, 0, 0, 0, 0};
        
        // Configuration
        static constexpr int HOT_THRESHOLD = 50; // Access count to promote to hot cache
        
        // Helper methods
        void cache_entire_plane(int apa, int face, int plane) const;
        int fetch_channel_from_anode(int apa, int face, int plane, int wire) const;
    

        // ----------------------------------------
        // Internal Storage
        // ----------------------------------------
        std::map<CoordReadout, ChargeMeasurement> m_charge_data;  ///< Internal charge data storage using ChargeMeasurement struct
        std::map<CoordReadout, ChargeMeasurement> m_orig_charge_data; // saved original charge measurement, if modified

        std::map<Coord2D, std::set<int>> m_2d_to_3d;  ///< Internal 2D→3D mapping
        std::map<int, Point3DInfo> m_3d_to_2d;               ///< Internal 3D→2D mapping
    
        // Global (apa, time, channel) to blobs
        std::map<CoordReadout, std::set<Facade::Blob* > > global_rb_map;

        // global geometry

        void BuildGeometry();

        std::map<WirePlaneId , std::tuple<WireCell::Point, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<WireCell::Point, double> > wpid_U_dir;
        std::map<WirePlaneId, std::pair<WireCell::Point, double> > wpid_V_dir;
        std::map<WirePlaneId, std::pair<WireCell::Point, double> > wpid_W_dir;
        std::set<int> apas;

        // Time_width, Pitch_u, pitch_v, pitch_w, for each apa/face
        std::map<WirePlaneId, std::tuple<double, double, double, double >> wpid_geoms;

        // geometry information T, U, V, W for each apa/face
        std::map<WirePlaneId, std::tuple<double, double, double, double >> wpid_offsets;
        // T, slope_yu slope_zu, slope_yv slope_zv, slope_yw slope_zw 
        std::map<WirePlaneId, std::tuple<double, std::pair<double, double>, std::pair<double, double>, std::pair<double, double> >> wpid_slopes;

        // result
        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
        std::vector<double> dQ;
        std::vector<double> dx;
        std::vector<double> pu;
        std::vector<double> pv;
        std::vector<double> pw;
        std::vector<double> pt;
        std::vector<std::pair<int, int>> paf;
        std::vector<double> reduced_chi2;
    };

} // namespace WireCell::Clus

#endif // WIRECELLCLUS_TRACKFITTING_H