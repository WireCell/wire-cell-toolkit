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
    
    T0Correction(IDetectorVolumes::pointer dv)
        : m_dv(dv) {
        
        for (const auto& [wfid, _] : m_dv->wpident_faces()) {
            WirePlaneId wpid(wfid);
            m_time_global_offsets[wpid.apa()][wpid.face()] = m_dv->metadata(wpid)["time_offset"].asDouble();
            m_drift_speeds[wpid.apa()][wpid.face()] = m_dv->metadata(wpid)["drift_speed"].asDouble();
        }
    }

    /**
     * From time2drift in Facade_Util.cxx
     * x_raw = xorig + face->dirx * (time_read_out + time_global_offset) * abs_drift_speed;
     * x_corr = xorig + face->dirx * (time_read_out - clustser_t0) * abs_drift_speed;
     *x_corr - x_raw = face->dirx * (- clustser_t0 - time_global_offset) * abs_drift_speed;
     */

    // get x_corr from x_raw
    virtual Point forward(const Point &pos_raw, double cluster_t0, int face,
                          int apa) const override        {
        Point pos_corr(pos_raw);
        pos_corr[0] -= m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (cluster_t0 ) *
            m_drift_speeds.at(apa).at(face);
        return pos_corr;
    }

    virtual Point backward(const Point &pos_corr, double cluster_t0, int face,
                           int apa) const override        {
        Point pos_raw(pos_corr);
        pos_raw[0] += m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (cluster_t0 ) *
            m_drift_speeds.at(apa).at(face);
        return pos_raw;
    }

    virtual bool filter(const Point &pos_corr, double clustser_t0, int face,
                        int apa) const override        {
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
        std::vector<double> arr_x_corr(arr_x.size());
        for (size_t i = 0; i < arr_x.size(); ++i) {
            arr_x_corr[i] = arr_x[i] - m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (cluster_t0 ) *
                m_drift_speeds.at(apa).at(face);
        }
        Dataset ds_corr;
        ds_corr.add(arr_cor_names[0], Array(arr_x_corr));
        ds_corr.add(arr_cor_names[1], Array(arr_y));
        ds_corr.add(arr_cor_names[2], Array(arr_z));
         
        //  ds_corr.add("x_corr", Array(arr_x_corr));
        //  ds_corr.add("y_corr", Array(arr_y));
        //  ds_corr.add("z_corr", Array(arr_z));
        return ds_corr;
    }

    virtual Dataset backward(const Dataset &pc_corr, const std::vector<std::string>& arr_cor_names, const std::vector<std::string>& arr_raw_names, double cluster_t0, int face,
                             int apa) const override        {
        const auto &arr_x = pc_corr.get(arr_cor_names[0])->elements<double>();
        const auto &arr_y = pc_corr.get(arr_cor_names[1])->elements<double>();
        const auto &arr_z = pc_corr.get(arr_cor_names[2])->elements<double>();
        std::vector<double> arr_x_corr(arr_x.size());
        for (size_t i = 0; i < arr_x.size(); ++i) {
            arr_x_corr[i] = arr_x[i] + m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa)) * (cluster_t0 ) *
                m_drift_speeds.at(apa).at(face);
        }
        Dataset ds_raw;
        ds_raw.add("x", Array(arr_x_corr));
        ds_raw.add("y", Array(arr_y));
        ds_raw.add("z", Array(arr_z));
        return ds_raw;
    }

    virtual Dataset filter(const Dataset &pc_corr, const std::vector<std::string>& arr_cor_names, double clustser_t0, int face,
                           int apa) const override        {
        std::vector<int> arr_filter(pc_corr.size_major());
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
 
private:
    IDetectorVolumes::pointer m_dv; // do not own
 
    // // m_time_global_offsets.at(apa).at(face) = time_global_offset
    std::map<int, std::map<int, double>> m_time_global_offsets;
    std::map<int, std::map<int, double>> m_drift_speeds;
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
        return cfg;
    }
    virtual void configure(const Configuration& cfg) {
        std::string dvtn = get<std::string>(cfg, "detector_volumes", "DetectorVolumes");
        auto dv = Factory::find_tn<WireCell::IDetectorVolumes>(dvtn);
        m_pcts["T0Correction"] = std::make_shared<T0Correction>(dv);
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
