#ifndef WIRECELLCLUS_TRACKFITTING_H
#define WIRECELLCLUS_TRACKFITTING_H

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Logging.h"
#include "WireCellClus/PRGraph.h"

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

        void add_segment(PR::Segment* segment); 

        // after the first round of track fitting, adjust the rough path ...
        WireCell::Point adjust_rough_path(PR::Segment& segment);

        // collect charge
        void collect_charge(double dis_cut, double range_cut);

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

    private:
        FittingType m_fitting_type;
        
        // cluster and grouping, CTPC is from m_grouping ...
        Facade::Grouping* m_grouping{nullptr};
        std::set<Facade::Cluster*> m_clusters;

        // input segment
        std::set<PR::Segment*> m_segments;

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

        // ----------------------------------------
        // Internal Storage
        // ----------------------------------------

        std::map<CoordReadout, ChargeMeasurement> m_charge_data;  ///< Internal charge data storage using ChargeMeasurement struct
        std::map<Coord2D, std::set<int>> m_2d_to_3d;  ///< Internal 2D→3D mapping
        std::map<int, Point3DInfo> m_3d_to_2d;               ///< Internal 3D→2D mapping
    
    
    };

} // namespace WireCell::Clus

#endif // WIRECELLCLUS_TRACKFITTING_H