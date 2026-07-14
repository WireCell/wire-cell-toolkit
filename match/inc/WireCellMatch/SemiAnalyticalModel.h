// SemiAnalyticalModel
//
// WCT port of larsim/PhotonPropagation/SemiAnalyticalModel{.h,.cxx}.
//
// Computes direct (VUV) and reflected (VIS) visibilities from a scintillation
// point to a list of optical detectors using the semi-analytical
// parameterisation. No larsoft / fhicl / art dependencies. Geometry and
// optical-detector arrays are supplied as plain structures by the caller
// (loaded from JSON for SBND in practice).
//
// Minimal scope (only what QLMatching exercises for SBND):
//   - dome PMTs (type=1) at anode/cathode orientation (0)
//   - flat PDs / (X)Arapucas (type=0) at anode/cathode orientation (0)
// The following branches are intentionally NOT ported:
//   - lateral PD VUV / VIS corrections
//   - anode reflections
//   - Xe absorption / vertical-border / field-cage corrections
//   - disk PMTs (type=2)
// If extra detector configurations are needed in the future, port the
// matching branches from larsim/PhotonPropagation/SemiAnalyticalModel.cxx.

#ifndef WIRECELL_MATCH_SEMIANALYTICALMODEL
#define WIRECELL_MATCH_SEMIANALYTICALMODEL

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Point.h"

#include <cstddef>
#include <vector>

namespace WireCell::Match {

    class SemiAnalyticalModel {
    public:
        struct Dims {
            double h;  // height [cm]
            double w;  // width  [cm]
        };

        struct OpticalDetector {
            double h;               // height [cm], -1 for dome
            double w;               // width  [cm], -1 for dome
            WireCell::Point center; // [cm]
            int type;               // 0=(X)Arapuca, 1=Dome PMT
            int orientation;        // 0=anode/cathode (only one supported)
        };

        // Detector-level geometry needed by the parameterisation.
        // All lengths in cm.
        struct Geometry {
            double active_center_y{0.0};
            double active_center_z{0.0};
            double active_size_y{0.0};
            double active_size_z{0.0};
            // Cathode plane sits at +/- |cathode_x|; the sign is taken from
            // the scintillation point's X. Cathode-centre Y/Z are taken from
            // active_center_*.
            double cathode_x{0.0};
            // VUV absorption length in LAr [cm]. (~85 cm @ 9.7 eV in larsoft
            // for SBND; tune per-detector or per-job.)
            double vuv_absorption_length{85.0};
        };

        // VUVHits and VISHits mirror the corresponding FHICL blocks from
        // larsim's semimodel_*.fcl. See sbnd/standalone-sample/
        // semi-analytical-sbnd.json (or equivalent) for the contents.
        SemiAnalyticalModel(const Configuration& VUVHits,
                            const Configuration& VISHits,
                            const Geometry& geom,
                            const std::vector<OpticalDetector>& opdets,
                            bool doReflectedLight = true);

        // opdet_mask (optional): when non-null, opdets with mask[i]==0 are skipped
        // (their visibility stays 0) instead of evaluating the per-opdet solid-angle
        // and Gaisser-Hillas corrections. The caller must size the mask to nOpDets();
        // a null mask (default) computes every opdet, leaving existing callers
        // bit-identical. Used by QLMatching to skip the masked PMTs / non-PMT opdets
        // whose predicted light it discards anyway.
        void detectedDirectVisibilities(std::vector<double>& vis,
                                        const WireCell::Point& scintPoint,
                                        const std::vector<unsigned int>* opdet_mask = nullptr) const;
        void detectedReflectedVisibilities(std::vector<double>& vis,
                                           const WireCell::Point& scintPoint,
                                           const std::vector<unsigned int>* opdet_mask = nullptr) const;

        std::size_t nOpDets() const { return m_opdets.size(); }

    private:
        // r is the precomputed border radius hypot(y-cy, z-cz) and distance the
        // precomputed point->opdet range; both are passed in so VUVVisibility does
        // not recompute them per opdet (bit-identical, same operands).
        double VUVVisibility(const WireCell::Point& scintPoint,
                             const OpticalDetector& opDet, double r, double distance) const;
        // rxy (cathode-plane border radius) and d_c (cathode distance) are passed in
        // (point-invariant across opdets) so VISVisibility does not recompute them
        // per opdet (bit-identical, same operands).
        double VISVisibility(const WireCell::Point& scintPoint,
                             const OpticalDetector& opDet,
                             double cathode_visibility,
                             const WireCell::Point& hotspot,
                             double rxy, double d_c) const;

        double Gaisser_Hillas(double x, const double* par) const;
        double Rectangle_SolidAngle(double a, double b, double d) const;
        double Rectangle_SolidAngle(const Dims& o,
                                    const WireCell::Vector& v,
                                    int OpDetOrientation) const;
        double Disk_SolidAngle(double d, double h, double b) const;
        double Omega_Dome_Model(double distance, double theta) const;

        // Geometry / opdets
        Geometry m_geom;
        std::vector<OpticalDetector> m_opdets;

        // VUV parameters
        bool m_isFlatPDCorr{false};
        bool m_isDomePDCorr{false};
        double m_delta_angulo_vuv{10.0};
        double m_pmt_radius{10.16}; // cm
        double m_maxPDDistance{1000.0};

        // flat PD (anode/cathode)
        std::vector<std::vector<double>> m_GHvuvpars_flat;
        std::vector<double> m_border_corr_angulo_flat;
        std::vector<std::vector<double>> m_border_corr_flat;
        // dome PD
        std::vector<std::vector<double>> m_GHvuvpars_dome;
        std::vector<double> m_border_corr_angulo_dome;
        std::vector<std::vector<double>> m_border_corr_dome;

        // VIS parameters
        bool m_doReflectedLight{true};
        double m_delta_angulo_vis{10.0};
        std::vector<double> m_vis_distances_x_flat;
        std::vector<double> m_vis_distances_r_flat;
        std::vector<std::vector<std::vector<double>>> m_vispars_flat;
        std::vector<double> m_vis_distances_x_dome;
        std::vector<double> m_vis_distances_r_dome;
        std::vector<std::vector<std::vector<double>>> m_vispars_dome;

        // Cathode plane (set when doReflectedLight)
        Dims m_cathode_plane{0., 0.};
        double m_plane_depth{0.0};
    };

} // namespace WireCell::Match

#endif
