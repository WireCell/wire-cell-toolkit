// SCECorrection.h
//
// A Clus::IPCTransform that performs combined T0 + SCE (Space Charge Effect)
// correction for SBND, delegating the SCE displacement field to an externally
// provided ISCEField.  This keeps clus/ ROOT-free; the field implementation
// (typically WireCell::Root::SCEFieldTH3) lives in the root/ subpackage and
// is looked up by jsonnet TypeName.
//
// Applies:
//   (1) T0:  x_t0  = x_raw - dirx * cluster_t0 * drift_speed
//   (2) SCE: x_sce = x_t0  + field->displacement_x(apa, x_t0, y, z)
//
// East/West (apa==0/1) routing is handled inside the ISCEField implementation.
// If no field is provided (nullptr), the SCE step is a no-op (T0 still applied).
//
// Author: Avinay Bhat (UChicago), for SBND
// Modeled on T0Correction by Haiwang Yu.

#ifndef WIRECELLCLUS_SCECORRECTION_H
#define WIRECELLCLUS_SCECORRECTION_H

#include "WireCellClus/IPCTransform.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ISCEField.h"
#include "WireCellUtil/Logging.h"

#include <map>
#include <memory>
#include <string>

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
            return {"3d", {"x_sce", "y", "z"}};
        }
        virtual std::vector<std::string> stored_array_names() const override {
            return {"x_sce"};
        }

    private:
        IDetectorVolumes::pointer m_dv;
        WireCell::ISCEField::pointer m_field;
        std::map<int, std::map<int, double>> m_drift_speeds;

        // Combined T0 + SCE for one X coordinate.  All quantities in WCT mm.
        inline double corrected_x(double x_raw_mm, double y_mm, double z_mm,
                                  double cluster_t0, int face, int apa) const {
            const double dirx = m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa));
            const double x_t0_mm = x_raw_mm - dirx * cluster_t0 * m_drift_speeds.at(apa).at(face);
            const double dx_mm = m_field ? m_field->displacement_x(apa, x_t0_mm, y_mm, z_mm) : 0.0;
            return x_t0_mm + dx_mm;
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
        if (m_field) spdlog::info("SCECorrection: ISCEField wired in");
        else         spdlog::info("SCECorrection: no ISCEField -- SCE disabled (T0 still applied)");
    }

    inline Point SCECorrection::forward(const Point& pos_raw, double t0, int face, int apa) const {
        Point pc(pos_raw);
        pc[0] = corrected_x(pos_raw[0], pos_raw[1], pos_raw[2], t0, face, apa);
        return pc;
    }

    inline Point SCECorrection::backward(const Point& pos_cor, double t0, int face, int apa) const {
        Point pr(pos_cor);
        if (m_field) {
            // Approximate: evaluate field at x_sce instead of x_t0 (small displacement OK).
            pr[0] -= m_field->displacement_x(apa, pos_cor[0], pos_cor[1], pos_cor[2]);
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
        std::vector<double> ax_c(ax.size());
        for (size_t i = 0; i < ax.size(); ++i) {
            ax_c[i] = corrected_x(ax[i], ay[i], az[i], t0, face, apa);
        }
        Dataset ds;
        ds.add(arr_cor_names[0], Array(ax_c));
        ds.add(arr_cor_names[1], Array(ay));
        ds.add(arr_cor_names[2], Array(az));
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
        std::vector<double> ax_r(ax.size());
        for (size_t i = 0; i < ax.size(); ++i) {
            double x = ax[i];
            if (m_field) x -= m_field->displacement_x(apa, ax[i], ay[i], az[i]);
            x += dirx * t0 * v;
            ax_r[i] = x;
        }
        Dataset ds;
        ds.add(arr_raw_names[0], Array(ax_r));
        ds.add(arr_raw_names[1], Array(ay));
        ds.add(arr_raw_names[2], Array(az));
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
