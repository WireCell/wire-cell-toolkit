// A Clus::IPCTransform that does T0 correction.
//
// It takes a single configuration parameter:
//
// detector_volumes which defaults to "DetectorVolumes".
//

#include "WireCellClus/IPCTransform.h"
#include "WireCellIface/IDetectorVolumes.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

#include <map>
#include <memory>
#include <string>

class PCTransformSet;

WIRECELL_FACTORY(PCTransformSet, PCTransformSet,
                 WireCell::Clus::IPCTransformSet,
                 WireCell::IConfigurable)


// Note, we do not register any of the individual IPCTransforms because the
// crazy clus code evolved to reinvent this pattern and it's too damn ugly at
// this point to try to refactor properly.


using namespace WireCell;
using namespace WireCell::Clus;

class T0Correction : public WireCell::Clus::IPCTransform
{
public:

    virtual ~T0Correction() = default;
    
    T0Correction(IDetectorVolumes::pointer dv, bool relax_containment_filter = false)
        : m_dv(dv)
        , m_relax_containment_filter(relax_containment_filter) {
        
        for (const auto& [wfid, _] : m_dv->wpident_faces()) {
            WirePlaneId wpid(wfid);
            const auto md = m_dv->metadata(wpid);
            m_time_global_offsets[wpid.apa()][wpid.face()] = md["time_offset"].asDouble();
            m_drift_speeds[wpid.apa()][wpid.face()] = md["drift_speed"].asDouble();
            // Per-event trigger offset applied ALONGSIDE cluster_t0 in the x
            // correction, for detectors that do NOT bake the readout-vs-trigger
            // offset into x_raw at imaging time (BlobSampler time_offset=0, e.g.
            // PDHD).  Distinct from "time_offset" above, which carries the value
            // baked into x_raw (and is intentionally NOT re-applied here).  Absent
            // key => 0 => production bit-identical.
            m_trigger_offsets[wpid.apa()][wpid.face()] =
                md.isMember("trigger_offset") ? md["trigger_offset"].asDouble() : 0.0;
            // Optional per-TPC transverse position offset (Y,Z); default zero.  The
            // x component (pos_offset[0]) is intentionally ignored -- the drift/x
            // correction stays with the t0/flash_x_offset term.  Presence of the
            // "pos_offset" key on ANY face flips the corrected scope to carry
            // y_cor/z_cor (so the post-QLMatching switch_scope materializes the
            // shift); absence => OFF => production bit-identical.  See
            // match/docs/cathode-offset-correction.md.
            double dy = 0.0, dz = 0.0;
            if (md.isMember("pos_offset") && md["pos_offset"].isArray() &&
                md["pos_offset"].size() >= 3) {
                dy = md["pos_offset"][1].asDouble();
                dz = md["pos_offset"][2].asDouble();
                m_has_pos_offset = true;
            }
            m_pos_offsets[wpid.apa()][wpid.face()] = {dy, dz};
        }
    }

    /**
     * From time2drift in Facade_Util.cxx
     * x_raw = xorig + face->dirx * (time_read_out + time_global_offset) * abs_drift_speed;
     * x_corr = xorig + face->dirx * (time_read_out - clustser_t0) * abs_drift_speed;
     *x_corr - x_raw = face->dirx * (- clustser_t0 - time_global_offset) * abs_drift_speed;
     *
     * Detectors that bake time_global_offset into x_raw (BlobSampler time_offset != 0,
     * e.g. SBND) carry the whole correction in cluster_t0, and trigger_offset = 0.
     * Detectors that do NOT bake it (BlobSampler time_offset = 0, e.g. PDHD) supply the
     * per-event readout-vs-trigger offset as "trigger_offset" in the DV metadata; it is
     * applied here ALONGSIDE cluster_t0:
     *   x_corr = x_raw - face->dirx * (cluster_t0 + trigger_offset) * abs_drift_speed;
     */

    // get x_corr from x_raw
    virtual Point forward(const Point &pos_raw, double cluster_t0, int face,
                          int apa) const override        {
        Point pos_corr(pos_raw);
        pos_corr[0] -= m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) *
            (cluster_t0 + m_trigger_offsets.at(apa).at(face)) *
            m_drift_speeds.at(apa).at(face);
        // Transverse (Y,Z) shift: a per-TPC constant, independent of cluster_t0.
        const auto& [dy, dz] = m_pos_offsets.at(apa).at(face);
        pos_corr[1] += dy;
        pos_corr[2] += dz;
        return pos_corr;
    }

    virtual Point backward(const Point &pos_corr, double cluster_t0, int face,
                           int apa) const override        {
        Point pos_raw(pos_corr);
        pos_raw[0] += m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) *
            (cluster_t0 + m_trigger_offsets.at(apa).at(face)) *
            m_drift_speeds.at(apa).at(face);
        // Invert the transverse (Y,Z) shift symmetrically with forward().
        const auto& [dy, dz] = m_pos_offsets.at(apa).at(face);
        pos_raw[1] -= dy;
        pos_raw[2] -= dz;
        return pos_raw;
    }

    virtual bool filter(const Point &pos_corr, double clustser_t0, int face,
                        int apa) const override        {
        // Relaxed mode (no-T0 detectors): accept every point.  Without an
        // event T0 the apparent x of out-of-time activity is unreliable in
        // both directions — early activity lands between the wire planes and
        // the sensitive-volume boundary, late activity past the cathode (gap
        // or opposite volume) — so any apparent-x containment test wrongly
        // excludes whole clusters from corrected-scope passes.
        if (m_relax_containment_filter) return true;
        auto wpid = m_dv->contained_by(pos_corr);
        if (!wpid.valid()) return false;
        if (wpid.apa() != apa || wpid.face() != face) return false;
        return true;
        //  return ().valid() ? true : false;
    }

    virtual Dataset forward(const Dataset &pc_raw, const std::vector<std::string>& arr_raw_names, const std::vector<std::string>& arr_cor_names, double cluster_t0, int face,
                            int apa) const override        {
        // std::cout << "Test: " << m_time_global_offsets.at(apa).at(face) << " " << cluster_t0 << std::endl;

        const auto &arr_x = pc_raw.get(arr_raw_names[0])->elements<double>();
        const auto &arr_y = pc_raw.get(arr_raw_names[1])->elements<double>();
        const auto &arr_z = pc_raw.get(arr_raw_names[2])->elements<double>();
        const auto& [dy, dz] = m_pos_offsets.at(apa).at(face);
        std::vector<double> arr_x_corr(arr_x.size());
        std::vector<double> arr_y_corr(arr_y.size());
        std::vector<double> arr_z_corr(arr_z.size());
        for (size_t i = 0; i < arr_x.size(); ++i) {
            arr_x_corr[i] = arr_x[i] - m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) *
                (cluster_t0 + m_trigger_offsets.at(apa).at(face)) *
                m_drift_speeds.at(apa).at(face);
            arr_y_corr[i] = arr_y[i] + dy;
            arr_z_corr[i] = arr_z[i] + dz;
        }
        Dataset ds_corr;
        ds_corr.add(arr_cor_names[0], Array(arr_x_corr));
        ds_corr.add(arr_cor_names[1], Array(arr_y_corr));
        ds_corr.add(arr_cor_names[2], Array(arr_z_corr));
        return ds_corr;
    }

    virtual Dataset backward(const Dataset &pc_corr, const std::vector<std::string>& arr_cor_names, const std::vector<std::string>& arr_raw_names, double cluster_t0, int face,
                             int apa) const override        {
        const auto &arr_x = pc_corr.get(arr_cor_names[0])->elements<double>();
        const auto &arr_y = pc_corr.get(arr_cor_names[1])->elements<double>();
        const auto &arr_z = pc_corr.get(arr_cor_names[2])->elements<double>();
        const auto& [dy, dz] = m_pos_offsets.at(apa).at(face);
        std::vector<double> arr_x_raw(arr_x.size());
        std::vector<double> arr_y_raw(arr_y.size());
        std::vector<double> arr_z_raw(arr_z.size());
        for (size_t i = 0; i < arr_x.size(); ++i) {
            arr_x_raw[i] = arr_x[i] + m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) *
                (cluster_t0 + m_trigger_offsets.at(apa).at(face)) *
                m_drift_speeds.at(apa).at(face);
            arr_y_raw[i] = arr_y[i] - dy;
            arr_z_raw[i] = arr_z[i] - dz;
        }
        Dataset ds_raw;
        ds_raw.add(arr_raw_names[0], Array(arr_x_raw));
        ds_raw.add(arr_raw_names[1], Array(arr_y_raw));
        ds_raw.add(arr_raw_names[2], Array(arr_z_raw));
        return ds_raw;
    }

    virtual Dataset filter(const Dataset &pc_corr, const std::vector<std::string>& arr_cor_names, double clustser_t0, int face,
                           int apa) const override        {
        std::vector<int> arr_filter(pc_corr.size_major());
        // See the Point overload: relaxed mode accepts every point.
        if (m_relax_containment_filter) {
            std::fill(arr_filter.begin(), arr_filter.end(), 1);
            Dataset ds;
            ds.add("filter", Array(arr_filter));
            return ds;
        }
        const auto &arr_x = pc_corr.get(arr_cor_names[0])->elements<double>();
        const auto &arr_y = pc_corr.get(arr_cor_names[1])->elements<double>();
        const auto &arr_z = pc_corr.get(arr_cor_names[2])->elements<double>();
        for (size_t i = 0; i < arr_x.size(); ++i) {
            arr_filter[i] = false;
            auto wpid = m_dv->contained_by(Point(arr_x[i], arr_y[i], arr_z[i]));
            if (wpid.valid()) {
                if (wpid.apa() == apa && wpid.face() == face) {
                    arr_filter[i] = true;
                }
            }
            // if (wpid.apa() != apa || wpid.face() != face) return false;
            //   ().valid() ? 1 : 0;
        }
        Dataset ds;
        ds.add("filter", Array(arr_filter));
        return ds;
    }
 
    virtual PointCloud::Tree::Scope output_scope() const override {
        // Without a transverse offset, the correction shifts x only; y and z are
        // unchanged and already exist in the blob's 3d PC under their original
        // names (production bit-identical).  With a per-TPC (Y,Z) offset, y and z
        // are shifted too, so the corrected scope names new y_cor/z_cor arrays.
        if (m_has_pos_offset) return {"3d", {"x_t0cor", "y_cor", "z_cor"}};
        return {"3d", {"x_t0cor", "y", "z"}};
    }

    virtual std::vector<std::string> stored_array_names() const override {
        // Persist only the arrays that actually changed.  Without a transverse
        // offset only x_t0cor is new; with one, y_cor/z_cor are new too.
        if (m_has_pos_offset) return {"x_t0cor", "y_cor", "z_cor"};
        return {"x_t0cor"};
    }

private:
    IDetectorVolumes::pointer m_dv; // do not own

    // m_time_global_offsets.at(apa).at(face) = time_global_offset
    std::map<int, std::map<int, double>> m_time_global_offsets;
    // m_trigger_offsets.at(apa).at(face) = per-event trigger offset applied with
    // cluster_t0 (un-baked detectors only; 0 => bit-identical).
    std::map<int, std::map<int, double>> m_trigger_offsets;
    std::map<int, std::map<int, double>> m_drift_speeds;
    // m_pos_offsets.at(apa).at(face) = {dy, dz} transverse position offset (cm-internal).
    std::map<int, std::map<int, std::pair<double, double>>> m_pos_offsets;
    // True iff any face declared a "pos_offset"; flips the corrected scope/stored
    // arrays to carry y_cor/z_cor.  False => OFF => bit-identical.
    bool m_has_pos_offset{false};
    // Disable filter() entirely — accept every point (for no-T0 detectors,
    // e.g. PDVD, where apparent x makes any containment test unreliable).
    // False => own-(apa,face) containment required => production bit-identical.
    bool m_relax_containment_filter{false};
};

class PCTransformSet : public WireCell::Clus::IPCTransformSet,
                       public WireCell::IConfigurable 
{
public:

    PCTransformSet() {}
    virtual ~PCTransformSet() {}
    
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["relax_containment_filter"] = false;
        return cfg;
    }
    virtual void configure(const Configuration& cfg) {
        std::string dvtn = get<std::string>(cfg, "detector_volumes", "DetectorVolumes");
        auto dv = Factory::find_tn<WireCell::IDetectorVolumes>(dvtn);
        const bool relax = get(cfg, "relax_containment_filter", false);
        m_pcts["T0Correction"] = std::make_shared<T0Correction>(dv, relax);
        // ...
    }

    virtual IPCTransform::pointer pc_transform(const std::string &name) const {
        auto it = m_pcts.find(name);
        if (it == m_pcts.end()) { return nullptr; }
        return it->second;
    }

private:
    std::map<std::string, IPCTransform::pointer> m_pcts;
};
