/* Comentary:

   This class is one element in a 2 by 2 matrix of classes which start
   very symmetric.

   (Depo, Impact) outer-product (Zipper, Transform)

   The "Depo" class is high level and uses the low level "Impact" class.

   The "Transform" version is newer than the "Zipper" version and
   while it started with a symmetric interface, it must grow out of
   it.

   The "Transform" version is faster than the "Zipper" for most
   realistically complicated event topology.  It will be made faster
   yet but this will break the symmetry:

   1) The "Tranform" version builds a dense array so sending its
   output through the narrow waveform() is silly.

   2) It can be optimized if the array has a channel basis instead of
   a wire one so the "zipper" interface is not fitting.

   3) It may be further optimized by performing its transforms on two
   planes from each face (for two-faced APAs) simultaneously.  This
   requires a different intrerface, possibly by (ab)using
   BinnedDiffusion in a different way than "zipper".

   4) It may be optimized yet further by breaking up the FFTs to do
   smaller ones (in time) on charge (X) FR and then doing the larger
   ones over ER(X)RC(X)RC.  This work also requires extension of PIR.
   This improvement may also be applicable to zipper.

   All of this long winded explantion, which nobody will ever read, is
   "justification" for the otherwise disgusting copy-paste which is
   going on with these two pairs of classes.

 */

#include "WireCellGen/DepoTransform.h"
#include "WireCellGen/ImpactTransform.h"
#include "WireCellGen/BinnedDiffusion_transform.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/DepoTools.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/NamedFactory.h"

#include <algorithm>
#include <vector>

WIRECELL_FACTORY(DepoTransform, WireCell::Gen::DepoTransform, WireCell::IDepoFramer, WireCell::IConfigurable)

using namespace WireCell;
using namespace std;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Gen::DepoTransform::DepoTransform()
  : Aux::Logger("DepoTransform", "gen")
  , m_start_time(0.0 * units::ns)
  , m_readout_time(5.0 * units::ms)
  , m_tick(0.5 * units::us)
  , m_drift_speed(1.0 * units::mm / units::us)
  , m_nsigma(3.0)
  , m_frame_count(0)
{
}

Gen::DepoTransform::~DepoTransform() {}

void Gen::DepoTransform::configure(const WireCell::Configuration& cfg)
{
    auto anode_tn = get<string>(cfg, "anode", "");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    m_nsigma = get<double>(cfg, "nsigma", m_nsigma);
    bool fluctuate = get<bool>(cfg, "fluctuate", false);
    m_rng = nullptr;
    if (fluctuate) {
        auto rng_tn = get<string>(cfg, "rng", "");
        m_rng = Factory::find_tn<IRandom>(rng_tn);
    }
    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_readout_time = get<double>(cfg, "readout_time", m_readout_time);
    m_tick = get<double>(cfg, "tick", m_tick);
    m_start_time = get<double>(cfg, "start_time", m_start_time);
    m_drift_speed = get<double>(cfg, "drift_speed", m_drift_speed);
    m_frame_count = get<int>(cfg, "first_frame_number", m_frame_count);

    m_process_planes = {0,1,2};

    if (cfg.isMember("process_planes")) {
      m_process_planes.clear();
      for (auto jplane : cfg["process_planes"]) {
	m_process_planes.push_back(jplane.asInt());
      }
    }

    log->debug("tick={} us, start={} us, readin={} us, drift_speed={} mm/us",
               m_tick/units::us, m_start_time/units::us,
               m_readout_time/units::us, m_drift_speed/(units::mm/units::us));

    auto jpirs = cfg["pirs"];
    if (jpirs.isNull() or jpirs.empty()) {
        std::string msg = "must configure with some plane impact response components";
        log->error(msg);
        THROW(ValueError() << errmsg{"Gen::Ductor: " + msg});
    }
    m_pirs.clear();
    for (auto jpir : jpirs) {
        auto tn = jpir.asString();
        auto pir = Factory::find_tn<IPlaneImpactResponse>(tn);
        m_pirs.push_back(pir);
    }
}
WireCell::Configuration Gen::DepoTransform::default_configuration() const
{
    Configuration cfg;

    /// How many Gaussian sigma due to diffusion to keep before truncating.
    put(cfg, "nsigma", m_nsigma);

    /// Whether to fluctuate the final Gaussian deposition.
    put(cfg, "fluctuate", false);

    /// The open a gate.  This is actually a "readin" time measured at
    /// the input ("response") plane.
    put(cfg, "start_time", m_start_time);

    /// The time span for each readout.  This is actually a "readin"
    /// time span measured at the input ("response") plane.
    put(cfg, "readout_time", m_readout_time);

    /// The sample period
    put(cfg, "tick", m_tick);

    /// The nominal speed of drifting electrons
    put(cfg, "drift_speed", m_drift_speed);

    /// Allow for a custom starting frame number
    put(cfg, "first_frame_number", m_frame_count);

    /// Name of component providing the anode plane.
    put(cfg, "anode", "");
    /// Name of component providing the anode pseudo random number generator.
    put(cfg, "rng", "");

    /// Plane impact responses
    cfg["pirs"] = Json::arrayValue;

    // type-name for the DFT to use
    cfg["dft"] = "FftwDFT";

    // Need to remove this otherwise [] will be the default
    // cfg["process_planes"] = Json::arrayValue;

    return cfg;
}

bool Gen::DepoTransform::operator()(const input_pointer& in, output_pointer& out)
{
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        out = nullptr;
        return true;
    }

    auto depos = in->depos();
    size_t ndepos_used=0;

    Binning tbins(m_readout_time / m_tick, m_start_time, m_start_time + m_readout_time);
    ITrace::vector traces;
    for (auto face : m_anode->faces()) {
        // Select the depos which are in this face's sensitive volume
        IDepo::vector face_depos = Aux::sensitive(*depos, face);
        ndepos_used += face_depos.size();

        int iplane = -1;
        for (auto plane : face->planes()) {
            ++iplane;

	    int plane_index = plane->planeid().index();

	    if (std::find(m_process_planes.begin(),  m_process_planes.end(), plane_index) == m_process_planes.end()) {
          log->debug("skip plane {}", plane_index);         
	      continue;
	    }

            const Pimpos* pimpos = plane->pimpos();

            Binning tbins(m_readout_time / m_tick, m_start_time, m_start_time + m_readout_time);

            Gen::BinnedDiffusion_transform bindiff(*pimpos, tbins, m_nsigma, m_rng);
            for (auto depo : face_depos) {
                depo = modify_depo(plane->planeid(), depo);
                bindiff.add(depo, depo->extent_long() / m_drift_speed, depo->extent_tran());
            }

            auto& wires = plane->wires();

            auto pir = m_pirs.at(iplane);
            Gen::ImpactTransform transform(pir, m_dft, bindiff);

            const int nwires = pimpos->region_binning().nbins();
            for (int iwire = 0; iwire < nwires; ++iwire) {
                auto wave = transform.waveform(iwire);

                auto mm = Waveform::edge(wave);
                if (mm.first == (int) wave.size()) {  // all zero
                    continue;
                }

                int chid = wires[iwire]->channel();
                int tbin = mm.first;

                ITrace::ChargeSequence charge(wave.begin() + mm.first, wave.begin() + mm.second);
                auto trace = make_shared<SimpleTrace>(chid, tbin, charge);
                traces.push_back(trace);
            }
            // fixme: use SPDLOG_LOGGER_DEBUG
            log->debug("plane={} face={} depos={} total traces={}",
                       iplane, face->ident(), face_depos.size(), traces.size());
        }
    }

    auto frame = make_shared<SimpleFrame>(m_frame_count, m_start_time, traces, m_tick);
    log->debug("call={} count={} ndepos_in={} ndepos_used={}",
               m_count, m_frame_count, depos->size(), ndepos_used);
    log->debug("output: {}", Aux::taginfo(frame));

    ++m_frame_count;
    ++m_count;
    out = frame;
    return true;
}
