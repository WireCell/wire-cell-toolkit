// SCECorrection.h
//
// A Clus::IPCTransform that performs combined T0 + SCE (Space Charge Effect)
// correction for SBND, delegating the SCE displacement field to an externally
// provided ISCEField.  This keeps clus/ ROOT-free; the field implementation
// (typically WireCell::Root::SCEFieldTH3) lives in the root/ subpackage and
// is looked up by jsonnet TypeName.
//
// Applies (all in WCT mm):
//   (1) T0:  x_t0  = x_raw - dirx * cluster_t0 * drift_speed
//   (2) SCE: x_sce = x_t0 + field->displacement_x(apa, x_t0, y, z)
//            y_sce = y    + field->displacement_y(apa, x_t0, y, z)
//            z_sce = z    + field->displacement_z(apa, x_t0, y, z)
//
// East/West (apa==0/1) routing is handled inside the ISCEField implementation.
// If no field is provided (nullptr), the SCE step is a no-op (T0 still applied
// to x; y, z pass through).
//
// Author: Avinay Bhat (UChicago), for SBND
// Modeled on T0Correction by Haiwang Yu.

#ifndef WIRECELLCLUS_SCECORRECTION_H
#define WIRECELLCLUS_SCECORRECTION_H

#include "WireCellClus/IPCTransform.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ISCEField.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace WireCell::Clus {

    class SCECorrection : public WireCell::Clus::IPCTransform {
    public:
        SCECorrection(IDetectorVolumes::pointer dv,
                      WireCell::ISCEField::pointer field = nullptr);
        virtual ~SCECorrection() = default;

        virtual Point forward(const Point& pos_raw,  double cluster_t0, int face, int apa) const override;
        virtual Point backward(const Point& pos_cor, double cluster_t0, int face, int apa) const override;
        virtual bool  filter(const Point& pos_cor,   double cluster_t0, int face, int apa) const override;

        virtual Dataset forward(const Dataset& pc_raw,
                                const std::vector<std::string>& arr_raw_names,
                                const std::vector<std::string>& arr_cor_names,
                                double cluster_t0, int face, int apa) const override;
        virtual Dataset backward(const Dataset& pc_cor,
                                 const std::vector<std::string>& arr_cor_names,
                                 const std::vector<std::string>& arr_raw_names,
                                 double cluster_t0, int face, int apa) const override;
        virtual Dataset filter(const Dataset& pc_cor,
                               const std::vector<std::string>& arr_cor_names,
                               double cluster_t0, int face, int apa) const override;

        virtual PointCloud::Tree::Scope output_scope() const override {
            return {"3d", {"x_sce", "y_sce", "z_sce"}};
        }
        virtual std::vector<std::string> stored_array_names() const override {
            return {"x_sce", "y_sce", "z_sce"};
        }

    private:
        IDetectorVolumes::pointer m_dv;
        WireCell::ISCEField::pointer m_field;
        std::map<int, std::map<int, double>> m_drift_speeds;

        // Combined T0 + full 3D SCE for one point.  All quantities in WCT mm.
        inline void corrected_xyz(double x_raw_mm, double y_mm, double z_mm,
                                  double cluster_t0, int face, int apa,
                                  double& xs, double& ys, double& zs) const {
            const double dirx = m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa));
            const double x_t0_mm = x_raw_mm - dirx * cluster_t0 * m_drift_speeds.at(apa).at(face);
            if (m_field) {
                xs = x_t0_mm + m_field->displacement_x(apa, x_t0_mm, y_mm, z_mm);
                ys = y_mm    + m_field->displacement_y(apa, x_t0_mm, y_mm, z_mm);
                zs = z_mm    + m_field->displacement_z(apa, x_t0_mm, y_mm, z_mm);
            } else {
                xs = x_t0_mm; ys = y_mm; zs = z_mm;
            }
        }
    };

    inline SCECorrection::SCECorrection(IDetectorVolumes::pointer dv,
                                        WireCell::ISCEField::pointer field)
        : m_dv(dv), m_field(field)
    {
        for (const auto& [wfid, _] : m_dv->wpident_faces()) {
            WirePlaneId wpid(wfid);
            const auto md = m_dv->metadata(wpid);
            m_drift_speeds[wpid.apa()][wpid.face()] = md["drift_speed"].asDouble();
        }
        // Use a named logger rather than the process-global default logger.
        auto log = Log::logger("clus.SCECorrection");
        if (m_field) {
            log->info("SCECorrection: ISCEField wired in (x,y,z)");
            // Sanity probe (once per construction): report the forward
            // (reco -> true) SCE displacement at a few fixed interior points so
            // the direction and magnitude (~cm) can be eyeballed in the log.
            // t0=0 -> pure SCE part (no drift term).
            const std::vector<Point> probes = {
                {-10 * units::cm,     100 * units::cm,        250 * units::cm},
                {-190 * units::cm,    100 * units::cm,        250 * units::cm},
                {10 * units::cm,      100 * units::cm,        250 * units::cm},
                {190 * units::cm,     100 * units::cm,        250 * units::cm},
            };
            for (const auto& p : probes) {
                WirePlaneId wpid = m_dv->contained_by(p);
                if (wpid.apa() < 0 || wpid.face() < 0) {
                    log->info("SCECorrection probe ({:.0f},{:.0f},{:.0f})cm: outside TPC",
                              p[0]/units::cm, p[1]/units::cm, p[2]/units::cm);
                    continue;
                }
                Point pc = forward(p, 0.0, wpid.face(), wpid.apa());
                double d = std::sqrt(std::pow(pc[0]-p[0],2) + std::pow(pc[1]-p[1],2) + std::pow(pc[2]-p[2],2));
                log->info("SCECorrection probe apa{} ({:.1f},{:.1f},{:.1f})cm -> "
                          "({:.1f},{:.1f},{:.1f})cm  |d|={:.2f}cm (forward = reco->true)",
                          wpid.apa(), p[0]/units::cm, p[1]/units::cm, p[2]/units::cm,
                          pc[0]/units::cm, pc[1]/units::cm, pc[2]/units::cm, d/units::cm);
            }
        }
        else         log->info("SCECorrection: no ISCEField -- SCE disabled (T0 still applied)");
    }

    inline Point SCECorrection::forward(const Point& pos_raw, double t0, int face, int apa) const {
        Point pc(pos_raw);
        corrected_xyz(pos_raw[0], pos_raw[1], pos_raw[2], t0, face, apa, pc[0], pc[1], pc[2]);
        return pc;
    }

    inline Point SCECorrection::backward(const Point& pos_cor, double t0, int face, int apa) const {
        Point pr(pos_cor);
        if (m_field) {
            // Approximate: evaluate field at the corrected position.
            pr[0] -= m_field->displacement_x(apa, pos_cor[0], pos_cor[1], pos_cor[2]);
            pr[1] -= m_field->displacement_y(apa, pos_cor[0], pos_cor[1], pos_cor[2]);
            pr[2] -= m_field->displacement_z(apa, pos_cor[0], pos_cor[1], pos_cor[2]);
        }
        const double dirx = m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa));
        pr[0] += dirx * t0 * m_drift_speeds.at(apa).at(face);
        return pr;
    }

    inline bool SCECorrection::filter(const Point& pos_cor, double, int face, int apa) const {
        auto wpid = m_dv->contained_by(pos_cor);
        return wpid.valid() && wpid.apa() == apa && wpid.face() == face;
    }

    inline PointCloud::Dataset SCECorrection::forward(const PointCloud::Dataset& pc_raw,
                                          const std::vector<std::string>& arr_raw_names,
                                          const std::vector<std::string>& arr_cor_names,
                                          double t0, int face, int apa) const
    {
        const auto& ax = pc_raw.get(arr_raw_names[0])->elements<double>();
        const auto& ay = pc_raw.get(arr_raw_names[1])->elements<double>();
        const auto& az = pc_raw.get(arr_raw_names[2])->elements<double>();
        std::vector<double> ax_c(ax.size()), ay_c(ax.size()), az_c(ax.size());
        for (size_t i = 0; i < ax.size(); ++i) {
            corrected_xyz(ax[i], ay[i], az[i], t0, face, apa,
                          ax_c[i], ay_c[i], az_c[i]);
        }
        Dataset ds;
        ds.add(arr_cor_names[0], Array(ax_c));   // x_sce
        ds.add(arr_cor_names[1], Array(ay_c));   // y_sce
        ds.add(arr_cor_names[2], Array(az_c));   // z_sce
        return ds;
    }

    inline PointCloud::Dataset SCECorrection::backward(const PointCloud::Dataset& pc_cor,
                                           const std::vector<std::string>& arr_cor_names,
                                           const std::vector<std::string>& arr_raw_names,
                                           double t0, int face, int apa) const
    {
        const auto& ax = pc_cor.get(arr_cor_names[0])->elements<double>();
        const auto& ay = pc_cor.get(arr_cor_names[1])->elements<double>();
        const auto& az = pc_cor.get(arr_cor_names[2])->elements<double>();
        const double dirx = m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa));
        const double v = m_drift_speeds.at(apa).at(face);
        std::vector<double> ax_r(ax.size()), ay_r(ax.size()), az_r(ax.size());
        for (size_t i = 0; i < ax.size(); ++i) {
            double x = ax[i], y = ay[i], z = az[i];
            if (m_field) {
                x -= m_field->displacement_x(apa, ax[i], ay[i], az[i]);
                y -= m_field->displacement_y(apa, ax[i], ay[i], az[i]);
                z -= m_field->displacement_z(apa, ax[i], ay[i], az[i]);
            }
            x += dirx * t0 * v;
            ax_r[i] = x; ay_r[i] = y; az_r[i] = z;
        }
        Dataset ds;
        ds.add(arr_raw_names[0], Array(ax_r));
        ds.add(arr_raw_names[1], Array(ay_r));
        ds.add(arr_raw_names[2], Array(az_r));
        return ds;
    }

    inline PointCloud::Dataset SCECorrection::filter(const PointCloud::Dataset& pc_cor,
                                         const std::vector<std::string>& arr_cor_names,
                                         double, int face, int apa) const
    {
        std::vector<int> flt(pc_cor.size_major());
        const auto& ax = pc_cor.get(arr_cor_names[0])->elements<double>();
        const auto& ay = pc_cor.get(arr_cor_names[1])->elements<double>();
        const auto& az = pc_cor.get(arr_cor_names[2])->elements<double>();
        for (size_t i = 0; i < ax.size(); ++i) {
            flt[i] = 0;
            auto wpid = m_dv->contained_by(Point(ax[i], ay[i], az[i]));
            if (wpid.valid() && wpid.apa() == apa && wpid.face() == face) flt[i] = 1;
        }
        Dataset ds;
        ds.add("filter", Array(flt));
        return ds;
    }

}  // namespace WireCell::Clus

#endif  // WIRECELLCLUS_SCECORRECTION_H
