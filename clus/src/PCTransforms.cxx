// A Clus::IPCTransform that does T0 correction.
//
// It takes a single configuration parameter:
//
// detector_volumes which defaults to "DetectorVolumes".
//
#include <TH3F.h>
#include <TFile.h>

#include <cmath>
#include <stdexcept>
#include <vector>

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

class SCECorrection : public WireCell::Clus::IPCTransform
{
public:
  virtual ~SCECorrection() = default;

  SCECorrection(IDetectorVolumes::pointer dv,
		const std::string& sce_file,
		double cathode_eps = 2.5)
    : m_dv(dv)
    , m_cathode_eps(cathode_eps)
  {
    std::unique_ptr<TFile> infile(TFile::Open(sce_file.c_str(), "READ"));
    if (!infile || infile->IsZombie()) {
      throw std::runtime_error("SCECorrection: failed to open SCE file: " + sce_file);
    }

    m_bkwd_e[0] = must_clone_hist(infile.get(), "TrueBkwd_Displacement_X_E");
    m_bkwd_e[1] = must_clone_hist(infile.get(), "TrueBkwd_Displacement_Y_E");
    m_bkwd_e[2] = must_clone_hist(infile.get(), "TrueBkwd_Displacement_Z_E");

    m_bkwd_w[0] = must_clone_hist(infile.get(), "TrueBkwd_Displacement_X_W");
    m_bkwd_w[1] = must_clone_hist(infile.get(), "TrueBkwd_Displacement_Y_W");
    m_bkwd_w[2] = must_clone_hist(infile.get(), "TrueBkwd_Displacement_Z_W");
  }

  virtual Point forward(const Point& pos_in, double cluster_t0, int face, int apa) const override
  {
    auto dpos = cal_offset(pos_in, face, apa);
    Point pos_out(pos_in);
    pos_out[0] += dpos[0];
    pos_out[1] += dpos[1];
    pos_out[2] += dpos[2];
    return pos_out;
  }

  virtual Point backward(const Point& pos_in, double cluster_t0, int face, int apa) const override
  {
    // First-pass approximate inverse.
    auto dpos = cal_offset(pos_in, face, apa);
    Point pos_out(pos_in);
    pos_out[0] -= dpos[0];
    pos_out[1] -= dpos[1];
    pos_out[2] -= dpos[2];
    return pos_out;
  }

  virtual bool filter(const Point& pos_corr, double cluster_t0, int face, int apa) const override
  {
    auto wpid = m_dv->contained_by(pos_corr);
    if (!wpid.valid()) return false;
    return (wpid.apa() == apa && wpid.face() == face);
  }

  virtual Dataset forward(const Dataset& pc_in,
			  const std::vector<std::string>& arr_in_names,
			  const std::vector<std::string>& arr_out_names,
			  double cluster_t0, int face, int apa) const override
  {
    const auto& arr_x = pc_in.get(arr_in_names[0])->elements<double>();
    const auto& arr_y = pc_in.get(arr_in_names[1])->elements<double>();
    const auto& arr_z = pc_in.get(arr_in_names[2])->elements<double>();

    std::vector<double> arr_x_out(arr_x.size());
    std::vector<double> arr_y_out(arr_y.size());
    std::vector<double> arr_z_out(arr_z.size());

    for (size_t i = 0; i < arr_x.size(); ++i) {
      Point pin(arr_x[i], arr_y[i], arr_z[i]);
      auto dpos = cal_offset(pin, face, apa);
      arr_x_out[i] = arr_x[i] + dpos[0];
      arr_y_out[i] = arr_y[i] + dpos[1];
      arr_z_out[i] = arr_z[i] + dpos[2];
    }

    Dataset ds_out;
    ds_out.add(arr_out_names[0], Array(arr_x_out));
    ds_out.add(arr_out_names[1], Array(arr_y_out));
    ds_out.add(arr_out_names[2], Array(arr_z_out));
    return ds_out;
  }

  virtual Dataset backward(const Dataset& pc_in,
			   const std::vector<std::string>& arr_in_names,
			   const std::vector<std::string>& arr_out_names,
			   double cluster_t0, int face, int apa) const override
  {
    const auto& arr_x = pc_in.get(arr_in_names[0])->elements<double>();
    const auto& arr_y = pc_in.get(arr_in_names[1])->elements<double>();
    const auto& arr_z = pc_in.get(arr_in_names[2])->elements<double>();

    std::vector<double> arr_x_out(arr_x.size());
    std::vector<double> arr_y_out(arr_y.size());
    std::vector<double> arr_z_out(arr_z.size());

    for (size_t i = 0; i < arr_x.size(); ++i) {
      Point pin(arr_x[i], arr_y[i], arr_z[i]);
      auto dpos = cal_offset(pin, face, apa);
      arr_x_out[i] = arr_x[i] - dpos[0];
      arr_y_out[i] = arr_y[i] - dpos[1];
      arr_z_out[i] = arr_z[i] - dpos[2];
    }

    Dataset ds_out;
    ds_out.add(arr_out_names[0], Array(arr_x_out));
    ds_out.add(arr_out_names[1], Array(arr_y_out));
    ds_out.add(arr_out_names[2], Array(arr_z_out));
    return ds_out;
  }

  virtual Dataset filter(const Dataset& pc_in,
			 const std::vector<std::string>& arr_names,
			 double cluster_t0, int face, int apa) const override
  {
    const auto& arr_x = pc_in.get(arr_names[0])->elements<double>();
    const auto& arr_y = pc_in.get(arr_names[1])->elements<double>();
    const auto& arr_z = pc_in.get(arr_names[2])->elements<double>();

    std::vector<int> arr_filter(arr_x.size(), 0);
    for (size_t i = 0; i < arr_x.size(); ++i) {
      auto wpid = m_dv->contained_by(Point(arr_x[i], arr_y[i], arr_z[i]));
      if (wpid.valid() && wpid.apa() == apa && wpid.face() == face) {
	arr_filter[i] = 1;
      }
    }

    Dataset ds;
    ds.add("filter", Array(arr_filter));
    return ds;
  }

private:
  IDetectorVolumes::pointer m_dv;
  double m_cathode_eps{2.5};

  TH3F* m_bkwd_e[3] = {nullptr, nullptr, nullptr};
  TH3F* m_bkwd_w[3] = {nullptr, nullptr, nullptr};

  static TH3F* must_clone_hist(TFile* tf, const char* name)
  {
    auto* h = dynamic_cast<TH3F*>(tf->Get(name));
    if (!h) {
      throw std::runtime_error(std::string("SCECorrection: missing histogram: ") + name);
    }
    auto* hc = dynamic_cast<TH3F*>(h->Clone());
    hc->SetDirectory(nullptr);
    return hc;
  }

  Point cal_offset(const Point& pin, int face, int apa) const
  {
    // SCE map axes are in cm; WCT internal length is mm. Convert.
    const double mm_per_cm = 10.0;
    double xx = pin[0] / mm_per_cm;
    double yy = pin[1] / mm_per_cm;
    double zz = pin[2] / mm_per_cm;

    if (xx < -199.999) xx = -199.999;
    else if (xx > 199.999) xx = 199.999;

    if (yy < -199.999) yy = -199.999;
    else if (yy > 199.999) yy = 199.999;

    if (zz < 0.001) zz = 0.001;
    else if (zz > 499.999) zz = 499.999;

    if (std::abs(xx) < m_cathode_eps) {
      const auto dirx = m_dv->face_dirx(WirePlaneId(kAllLayers, face, apa));
      xx = (dirx < 0 ? -m_cathode_eps : m_cathode_eps);
    }

    TH3F** hs = (xx < 0.0) ? const_cast<TH3F**>(m_bkwd_e)
      : const_cast<TH3F**>(m_bkwd_w);

    return Point(hs[0]->Interpolate(xx, yy, zz) * mm_per_cm,
		 hs[1]->Interpolate(xx, yy, zz) * mm_per_cm,
		 hs[2]->Interpolate(xx, yy, zz) * mm_per_cm);
  }
};

class T0SCECorrection : public WireCell::Clus::IPCTransform
{
public:
  virtual ~T0SCECorrection() = default;

  T0SCECorrection(IDetectorVolumes::pointer dv,
		  const std::string& sce_file,
		  double cathode_eps = 2.5)
    : m_t0(dv)
    , m_sce(dv, sce_file, cathode_eps)
  {}

  virtual Point forward(const Point& pos_in, double cluster_t0, int face, int apa) const override
  {
    auto p1 = m_t0.forward(pos_in, cluster_t0, face, apa);
    return m_sce.forward(p1, cluster_t0, face, apa);
  }

  virtual Point backward(const Point& pos_in, double cluster_t0, int face, int apa) const override
  {
    auto p1 = m_sce.backward(pos_in, cluster_t0, face, apa);
    return m_t0.backward(p1, cluster_t0, face, apa);
  }

  virtual bool filter(const Point& pos_corr, double cluster_t0, int face, int apa) const override
  {
    return m_sce.filter(pos_corr, cluster_t0, face, apa);
  }

  virtual Dataset forward(const Dataset& pc_in,
			  const std::vector<std::string>& arr_in_names,
			  const std::vector<std::string>& arr_out_names,
			  double cluster_t0, int face, int apa) const override
  {
    auto tmp = m_t0.forward(pc_in, arr_in_names, arr_out_names, cluster_t0, face, apa);
    return m_sce.forward(tmp, arr_out_names, arr_out_names, cluster_t0, face, apa);
  }

  virtual Dataset backward(const Dataset& pc_in,
			   const std::vector<std::string>& arr_in_names,
			   const std::vector<std::string>& arr_out_names,
			   double cluster_t0, int face, int apa) const override
  {
    auto tmp = m_sce.backward(pc_in, arr_in_names, arr_out_names, cluster_t0, face, apa);
    return m_t0.backward(tmp, arr_out_names, arr_out_names, cluster_t0, face, apa);
  }

  virtual Dataset filter(const Dataset& pc_in,
			 const std::vector<std::string>& arr_names,
			 double cluster_t0, int face, int apa) const override
  {
    return m_sce.filter(pc_in, arr_names, cluster_t0, face, apa);
  }

private:
  T0Correction m_t0;
  SCECorrection m_sce;
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
	cfg["enable_sce_correction"] = false;
	cfg["sce_file"] = "";
	cfg["sce_cathode_eps"] = 2.5;
        return cfg;
    }
    virtual void configure(const Configuration& cfg) {
        std::string dvtn = get<std::string>(cfg, "detector_volumes", "DetectorVolumes");
        auto dv = Factory::find_tn<WireCell::IDetectorVolumes>(dvtn);
        m_pcts["T0Correction"] = std::make_shared<T0Correction>(dv);
	const bool enable_sce = get<bool>(cfg, "enable_sce_correction", false);
	const std::string sce_file = get<std::string>(cfg, "sce_file", "");
	const double sce_cathode_eps = get<double>(cfg, "sce_cathode_eps", 2.5);

	if (enable_sce) {
	  if (sce_file.empty()) {
            throw std::runtime_error("PCTransformSet: enable_sce_correction=true but sce_file is empty");
	  }
	  m_pcts["SCECorrection"] = std::make_shared<SCECorrection>(dv, sce_file, sce_cathode_eps);
	  m_pcts["T0SCECorrection"] = std::make_shared<T0SCECorrection>(dv, sce_file, sce_cathode_eps);
	}
    }

    virtual IPCTransform::pointer pc_transform(const std::string &name) const {
        auto it = m_pcts.find(name);
        if (it == m_pcts.end()) { return nullptr; }
        return it->second;
    }

private:
    std::map<std::string, IPCTransform::pointer> m_pcts;
};
