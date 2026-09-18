#include "WireCellGen/DepoFluxSplat.h"

#include "WireCellIface/IFieldResponse.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellUtil/Range.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "WireCellGen/GaussianDiffusion.h"

#include "WireCellUtil/NumpyHelper.h"

#include "WireCellUtil/RayTiling.h"

#include <functional>           //  std::plus
#include <array>
#include <algorithm>
#include <unordered_map>
#include <limits>

WIRECELL_FACTORY(DepoFluxSplat,
                 WireCell::Gen::DepoFluxSplat,
                 WireCell::INamed,
                 WireCell::IConfigurable,
                 WireCell::IDepoFramer)

using namespace WireCell;
using namespace WireCell::Aux;
using WireCell::Range::irange;


Gen::DepoFluxSplat::DepoFluxSplat()
    : Aux::Logger("DepoFluxSplat", "gen")
{
}


Gen::DepoFluxSplat::~DepoFluxSplat()
{
}


WireCell::Configuration Gen::DepoFluxSplat::default_configuration() const
{
    Configuration cfg;

    // Accept array of strings or single string
    cfg["anode"] = Json::arrayValue;
    cfg["field_response"] = Json::stringValue;

    // time binning
    const double default_tick = 0.5 * units::us;
    cfg["tick"] = default_tick;
    cfg["window_start"] = 0;
    cfg["window_duration"] = 8096 * default_tick;

    cfg["sparse"] = m_sparse;

    cfg["nsigma"] = m_nsigma;

    cfg["reference_time"] = 0;
    cfg["time_offsets"] = Json::arrayValue;

    cfg["smear_long"] = 0.0;
    cfg["smear_tran"] = 0.0;

    // Depo-to-channel association output.  Empty disables it.
    cfg["trio_file"] = m_trio_file;
    cfg["trio_min_charge"] = m_trio_min_charge;

    return cfg;
}

namespace {

    // A "trio" is one channel in each plane, at its tick, known to share a
    // cause because a single deposition produced all three.  The frame sums
    // every depo's contribution and so discards which depo went where.
    //
    // Key is (tbin, chan_u, chan_v, chan_w) with tbin the NOMINAL tick, before
    // m_tick_offsets is applied.  The offsets are per plane, so the three
    // channels of a trio need not be hot at the same tick -- but the offsets
    // are constants of the job, so one nominal tick plus the offsets written
    // once carries the same information as three ticks per row, at a third of
    // the width.  Storing three would be exact but wasteful; ASSUMING they are
    // equal would be wrong the moment anyone sets time_offsets.
    using TrioKey = std::array<int, 4>;

    struct TrioKeyHash {
        size_t operator()(const TrioKey& k) const noexcept {
            size_t h = 1469598103934665603ull;          // FNV-1a
            for (int v : k) {
                h ^= static_cast<size_t>(static_cast<unsigned>(v));
                h *= 1099511628211ull;
            }
            return h;
        }
    };

    // unordered, not ordered: a full-drift depo contributes ~10^2 rows and an
    // event ~5*10^7 insertions, where std::map's log(n) compare chain costs
    // real time.  Sorted only once, at write.
    using TrioMap = std::unordered_map<TrioKey, double, TrioKeyHash>;

    // One plane's view of one depo: the patch this component already computed
    // for the frame, kept so the three planes can be combined after the plane
    // loop has run.
    struct PlanePatch {
        bool filled{false};
        int pbeg{0}, nwire{0};      // wire index range
        int tbeg{0}, ntick{0};      // tick range, BEFORE the per-plane offset
        int toff{0};                // m_tick_offsets for this plane
        std::vector<float> q;       // nwire x ntick, row major
        const IWirePlane* plane{nullptr};

        float charge(int iw, int it) const { return q[iw*ntick + it]; }
    };

    // Turn one depo's three patches into trios.
    //
    // The patch is a separable Gaussian in (pitch, time), so its wire set does
    // not vary with tick and the overlap test runs ONCE per depo rather than
    // once per tick.  RayGrid tiling is what decides which wire combinations
    // genuinely overlap; taking the cross product of the three planes' wire
    // ranges without it would assert correspondences that no electron realises.
    void add_trios(const IAnodeFace::pointer& face,
                   const std::array<PlanePatch, 3>& pp,
                   double min_charge,
                   TrioMap& trios,
                   size_t& nkept, size_t& ndropped)
    {
        using namespace WireCell::RayGrid;

        for (const auto& p : pp) {
            if (!p.filled) { ++ndropped; return; }
        }

        // Layers 0 and 1 are RayGrid's horizontal/vertical bounds, as in
        // img/src/GridTiling.cxx; the wire planes follow.
        activities_t activities;
        activities.push_back(Activity(0, 1, 1.0));
        activities.push_back(Activity(1, 1, 1.0));
        for (int i = 0; i < 3; ++i) {
            activities.push_back(Activity(2+i, (size_t)pp[i].nwire, 1.0, pp[i].pbeg));
        }

        auto blobs = make_blobs(face->raygrid(), activities);
        if (blobs.empty()) { ++ndropped; return; }

        // Ticks are shared across planes -- one drift time, seen by all three --
        // so the correspondence runs over the ticks they have in common, with
        // each plane's own offset applied only when forming the key.
        const int t0 = std::max({pp[0].tbeg, pp[1].tbeg, pp[2].tbeg});
        const int t1 = std::min({pp[0].tbeg + pp[0].ntick,
                                 pp[1].tbeg + pp[1].ntick,
                                 pp[2].tbeg + pp[2].ntick});
        if (t0 >= t1) { ++ndropped; return; }

        // Collect the wire triples the tiling admits.
        std::vector<std::array<int,3>> triples;
        for (const auto& blob : blobs) {
            std::array<grid_range_t, 3> wr;
            bool ok = true;
            for (int i = 0; i < 3; ++i) {
                const auto& strips = blob.strips();
                auto it = std::find_if(strips.begin(), strips.end(),
                                       [i](const Strip& s) { return s.layer == (layer_index_t)(2+i); });
                if (it == strips.end()) { ok = false; break; }
                wr[i] = it->bounds;
            }
            if (!ok) continue;
            // A blob's strips are bounds on the tiled region and can reach past
            // the depo's own patch.  Clip to it: outside, this depo deposits no
            // charge, and the wire index would not even be safe to look up.
            for (int i = 0; i < 3; ++i) {
                wr[i].first  = std::max(wr[i].first,  pp[i].pbeg);
                wr[i].second = std::min(wr[i].second, pp[i].pbeg + pp[i].nwire);
                if (wr[i].first >= wr[i].second) { ok = false; break; }
            }
            if (!ok) continue;
            for (int u = wr[0].first; u < wr[0].second; ++u)
                for (int v = wr[1].first; v < wr[1].second; ++v)
                    for (int w = wr[2].first; w < wr[2].second; ++w)
                        triples.push_back({u, v, w});
        }
        if (triples.empty()) { ++ndropped; return; }
        // Blobs may overlap at their corners, so the same triple can appear twice.
        std::sort(triples.begin(), triples.end());
        triples.erase(std::unique(triples.begin(), triples.end()), triples.end());

        for (int t = t0; t < t1; ++t) {
            // The depo's charge at this tick, and the per-triple weights.  The
            // three planes see the same electrons, so their tick totals agree
            // up to per-plane smearing; average them rather than privileging one.
            double qtick = 0;
            for (int i = 0; i < 3; ++i) {
                double s = 0;
                for (int iw = 0; iw < pp[i].nwire; ++iw) s += pp[i].charge(iw, t - pp[i].tbeg);
                qtick += s;
            }
            qtick /= 3.0;
            if (qtick <= 0) continue;

            // Weight a triple by the product of the three planes' charge on its
            // wires -- the independent-marginals estimate of the joint, which is
            // all the per-plane patches can support -- then renormalise so the
            // triples at this tick carry exactly the depo's charge.
            std::vector<double> wgt(triples.size(), 0.0);
            double norm = 0;
            for (size_t k = 0; k < triples.size(); ++k) {
                double prod = 1.0;
                for (int i = 0; i < 3; ++i) {
                    const int iw = triples[k][i] - pp[i].pbeg;
                    if (iw < 0 || iw >= pp[i].nwire) { prod = 0; break; }
                    prod *= pp[i].charge(iw, t - pp[i].tbeg);
                }
                wgt[k] = prod;
                norm += prod;
            }
            if (norm <= 0) continue;

            for (size_t k = 0; k < triples.size(); ++k) {
                const double q = qtick * wgt[k] / norm;
                // A triple the depo puts no charge on is not a trio.  Note
                // min_charge defaults to 0, so this must be its own test --
                // "q < 0" would let every zero-weight triple through.
                if (q <= 0 || q < min_charge) continue;
                TrioKey key;
                key[0] = t;             // nominal tick; offsets applied by the reader
                for (int i = 0; i < 3; ++i) {
                    key[1+i] = pp[i].plane->wires()[triples[k][i]]->channel();
                }
                trios[key] += q;
            }
        }
        ++nkept;
    }
}

static
std::vector<double> get_n(const Configuration& cfg, size_t n=3)
{
    if (cfg.isDouble()) {
        return std::vector<double>(n, cfg.asDouble());
    }
    if (cfg.size() == n) {
        std::vector<double> ret;
        for (auto const& one : cfg) {
            ret.push_back(one.asDouble());
        }
        return ret;
    }
    return std::vector<double>(n, 0);
}

void Gen::DepoFluxSplat::configure(const WireCell::Configuration& cfg)
{
    // For response plane info
    auto jfrname = cfg["field_response"];
    if (! jfrname.isNull()) {
        auto ifr = Factory::find_tn<IFieldResponse>(jfrname.asString());
        const auto& fr = ifr->field_response();
        m_speed = fr.speed;
        m_origin = fr.origin;
    }
    auto jspeed = cfg["drift_speed"];
    if (! jspeed.isNull()) {
        m_speed = jspeed.asDouble();
    }
    auto jorigin = cfg["response_plane"];
    if (! jorigin.isNull()) {
        m_origin = jorigin.asDouble();
    }

    // Anode plane for down-selecting depos.
    std::string anode_tn = cfg["anode"].asString();
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    // Acceptance window 
    const double wtick = get(cfg, "tick", 0.5 * units::us);
    const double wstart = cfg["window_start"].asDouble();
    const double wduration = cfg["window_duration"].asDouble();
    const int nwbins = wduration / wtick;
    m_tbins = Binning(nwbins, wstart, wstart + nwbins * wtick);

    // Frame form.
    m_sparse = get(cfg, "sparse", m_sparse);

    // Gaussian cut-off.
    m_nsigma = get(cfg, "nsigma", m_nsigma);

    // Depo-to-channel association output.
    m_trio_file = get<std::string>(cfg, "trio_file", m_trio_file);
    m_trio_min_charge = get(cfg, "trio_min_charge", m_trio_min_charge);

    //Check which plane to work on
    m_process_planes = {0,1,2};

    if (cfg["process_planes"].isArray()) {
	m_process_planes.clear();
	for (auto jplane : cfg["process_planes"]) {
	    m_process_planes.push_back(jplane.asInt());
	}
    }

    // Additional smearing.
    m_smear_long = get_n(cfg["smear_long"]);
    m_smear_tran = get_n(cfg["smear_tran"]);

    // Arbitrary time subtracted from window_start when setting frame time
    m_reftime = get(cfg, "reference_time", m_reftime);

    // Arbitrary time added to tbin of traces on per plane basis.
    m_tick_offsets.clear();
    m_tick_offsets.resize(3,0);
    auto jto = cfg["time_offsets"];
    if (jto.isArray()) {
        if (jto.size() == 3) {
            for (int ind = 0; ind < 3; ++ind) {
                m_tick_offsets[ind] = jto[ind].asDouble() / wtick;
            }
        }
        else if (!jto.empty()) {
            THROW(ValueError() << errmsg{"DepoFluxSplat: time_offsets must be empty or be a 3-array"});
        }
    }
    log->debug("speed={} mm/us, origin={} mm, tbins: {} {}us ticks: [{},{}]us, reftime={} us, tick offsets=({},{},{})",
               m_speed / (units::mm/units::us), m_origin/units::mm,
               m_tbins.nbins(), m_tbins.binsize()/units::us,
               m_tbins.min()/units::us, m_tbins.max()/units::us,
               m_reftime/units::us,
               m_tick_offsets[0],m_tick_offsets[1],m_tick_offsets[2]);

}

IAnodeFace::pointer Gen::DepoFluxSplat::find_face(const IDepo::pointer& depo)
{
    for (auto face : m_anode->faces()) {
        auto bb = face->sensitive();
        if (bb.inside(depo->pos())) { return face; }
    }
    return nullptr;
}

// Return intersection of two half open ranges.
using intrange_t = std::pair<int,int>;

// Return intersection of half-open ranges r1 and r2.  If r1 is fully before r2,
// return [r2.first,r2.first] and if fully after return [r2.second.r2.second].
static intrange_t intersect(intrange_t const& r1, intrange_t const& r2)
{
    if (r1.second <= r2.first) {
        return std::make_pair(r2.first, r2.first);
    }
    if (r1.first >= r2.second) {
        return std::make_pair(r2.second, r2.second);
    }
    return std::make_pair(std::max(r1.first, r2.first),
                          std::min(r1.second, r2.second));
}

// A base class for sparse vs dense accumulation
struct Accumulator {

    virtual ~Accumulator() {};
    // Accumulate one depo's patch
    virtual void add(int chid, int tbin, const std::vector<float>& charge) = 0;
    virtual IFrame::pointer frame(int ident, double time, double tick) = 0;
    virtual size_t ntraces() const = 0;
};


// Accumulate in a sparse way.
struct SparseAccumulator : public Accumulator {
    ITrace::vector itraces;

    virtual ~SparseAccumulator() {}
    virtual void add(int chid, int tbin, const std::vector<float>& charge) {
        itraces.push_back(std::make_shared<SimpleTrace>(chid,tbin,charge));
    }
    virtual IFrame::pointer frame(int ident, double time, double tick) {
        return std::make_shared<SimpleFrame>(ident, time, itraces, tick);
    }
    virtual size_t ntraces() const { return itraces.size(); }
};

static intrange_t make_intrange(int beg, size_t siz)
{
    return std::make_pair(beg, beg+siz);
}
static intrange_t union_intrange(intrange_t const& r1, intrange_t const& r2)
{
    return std::make_pair(std::min(r1.first, r2.first),
                          std::max(r1.second, r2.second));
}

// Accumulate in a dense way.
struct DenseAccumulator : public Accumulator {
    using SharedSimpleTrace = std::shared_ptr<SimpleTrace>;
    std::map<int, SharedSimpleTrace> traces; // map chid->trace
    virtual void add(int chid, int tbin, const std::vector<float>& charge) {
        auto tp = traces[chid];
        if (!tp) {              // first seen
            traces[chid] = std::make_shared<SimpleTrace>(chid,tbin,charge);
            return;
        }
        const auto& oldcharge = tp->charge();
        auto have = make_intrange(tp->tbin(), oldcharge.size());
        auto want = make_intrange(tbin, charge.size());
        auto need = union_intrange(have, want);
        if (need.first < have.first || need.second > have.second) {
            std::vector<float> newcharge(need.second-need.first, 0);
            std::copy(oldcharge.begin(), oldcharge.end(),
                      newcharge.begin() + have.first-need.first);
            std::transform(charge.begin(), charge.end(),
                           newcharge.begin() + want.first-need.first,
                           newcharge.begin() + want.first-need.first,
                           std::plus<float>());
            tp->charge() = newcharge;
        }
    }

    virtual IFrame::pointer frame(int ident, double time, double tick) {
        ITrace::vector itraces;
        for (const auto& [chid, tp] : traces) {
            itraces.push_back(tp);
        }
        return std::make_shared<SimpleFrame>(ident, time, itraces, tick); 
    }        
    virtual size_t ntraces() const { return traces.size(); }
};



bool Gen::DepoFluxSplat::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at {}", m_count);
        ++m_count;
        return true;            // EOS
    }

    std::unique_ptr<Accumulator> accum;
    if (m_sparse) {
        accum = std::make_unique<SparseAccumulator>();
    }
    else {
        accum = std::make_unique<DenseAccumulator>();
    }

    size_t ndepos_seen=0;
    // size_t ndepos_skipped=0;
    size_t nplanes_skipped=0;

    // Depo-to-channel association, collected only when asked for.
    const bool want_trios = !m_trio_file.empty();
    TrioMap trios;
    size_t ntrio_depos=0, ntrio_incomplete=0;

    for (const auto& depo : *in->depos()) {
        if (!depo) {
            // ++ndepos_skipped;
            continue;
        }

        auto face = find_face(depo);
        if (!face) {
            // ++ndepos_skipped;
            continue;
        }

        ++ndepos_seen;

        // Depo is at response plane.  Find its time at the collection
        // plane assuming it were to continue along a uniform field.
        // After this, all times are nominal up until we add arbitrary
        // time offsets at output.
        const double nominal_depo_time = depo->time() + m_origin / m_speed;

        // This depo's patch in each plane, gathered as the plane loop runs and
        // consumed after it.
        std::array<PlanePatch, 3> patches;

        // Tabulate depo flux for wire regions from each plane
        for (auto plane : face->planes()) {

            int iplane = plane->planeid().index();
            if (iplane < 0) {
                ++nplanes_skipped;
                continue;
            }
	    
	    if (std::find(m_process_planes.begin(),  m_process_planes.end(), iplane) == m_process_planes.end()) {
		continue;
            }

            // Allow for extra smear in time unique to each plane.
            double sigma_L = depo->extent_long(); // [length]
            const double smear_long = m_smear_long[iplane];
            if (smear_long > 0) {
                const double extra = smear_long * m_tbins.binsize() * m_speed;
                sigma_L = sqrt(sigma_L * sigma_L + extra * extra);
            }
            Gen::GausDesc time_desc(nominal_depo_time, sigma_L / m_speed);

            // Check if patch is fully outside time binning
            {
                double nmin_sigma = time_desc.distance(m_tbins.min());
                double nmax_sigma = time_desc.distance(m_tbins.max());

                double eff_nsigma = depo->extent_long() > 0 ? m_nsigma : 0;
                if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) {
                    // ++ndepos_skipped;
                    continue;
                }
            }

            const Pimpos* pimpos = plane->pimpos();
            auto& wires = plane->wires();
            auto wbins = pimpos->region_binning(); // wire binning

            double sigma_T = depo->extent_tran();
            const double smear_tran = m_smear_tran[iplane];
            if (smear_tran > 0) {
                const double extra = smear_tran * wbins.binsize();
                sigma_T = sqrt(sigma_T * sigma_T + extra * extra);
            }

            const double center_pitch = pimpos->distance(depo->pos());
            Gen::GausDesc pitch_desc(center_pitch, sigma_T);
            {
                double nmin_sigma = pitch_desc.distance(wbins.min());
                double nmax_sigma = pitch_desc.distance(wbins.max());

                double eff_nsigma = depo->extent_tran() > 0 ? m_nsigma : 0;
                if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) {
                    // We are more than "N-sigma" outside the wire plane.
                    break;
                }
            }

            // The heavy lifting
            Gen::GaussianDiffusion gd(depo, time_desc, pitch_desc);
            gd.set_sampling(m_tbins, wbins, m_nsigma, 0, 1);

            // Transfer depo's patch to itraces
            const auto patch = gd.patch(); // 2D array


            // The absolute pitch bin for the first row of the patch array.
            const int pbin0 = gd.poffset_bin();

            // Absolute pitch bin range truncated to fit wire plane.
            const auto p_range = intersect({pbin0, pbin0 + patch.rows()},
                                           {0, wbins.nbins()});
            if (p_range.first == p_range.second) {
                ++nplanes_skipped;
                continue;
            }

            // The absolute tick bin for the first column of the patch array
            const int tbin0 = gd.toffset_bin();

            // Absolute tick bin range truncated to fit time window.
            const auto t_range = intersect({tbin0, tbin0 + patch.cols()},
                                           {0, m_tbins.nbins()});
            if (t_range.first == t_range.second) {
                ++nplanes_skipped;
                continue;
            }

            // Iterate over the valid wires covered by the patch
            for (int pbin : irange(p_range)) {
                auto iwire = wires[pbin];
                const int chid = iwire->channel();

                const int ncharges = t_range.second - t_range.first;
                std::vector<float> charge(ncharges);
                
                const int prel = pbin          - pbin0;
                const int trel = t_range.first - tbin0;
                
                // truly cursed
                Eigen::VectorXf::Map(&charge[0], ncharges) = patch.row(prel).segment(trel, charge.size());

                // Differing conventions exist for the sign of the charge of
                // "number of electrons" in the depo.  Force positive signal.
                std::transform(charge.begin(), charge.end(), charge.begin(),
                               static_cast<float (*)(float)>(&std::abs));
                
                accum->add(chid, t_range.first + m_tick_offsets[iplane], charge);

                if (want_trios) {
                    // Keep this plane's patch so the three can be combined
                    // after the plane loop.  tbeg is stored BEFORE the
                    // per-plane tick offset: the ticks correspond across
                    // planes in nominal time, and the offset is applied only
                    // when the key is formed.
                    auto& pp = patches[iplane];
                    if (!pp.filled) {
                        pp.filled = true;
                        pp.plane = plane.get();
                        pp.pbeg = p_range.first;
                        pp.nwire = p_range.second - p_range.first;
                        pp.tbeg = t_range.first;
                        pp.ntick = ncharges;
                        pp.toff = m_tick_offsets[iplane];
                        pp.q.assign((size_t)pp.nwire * ncharges, 0.0f);
                    }
                    const int iw = pbin - pp.pbeg;
                    if (iw >= 0 && iw < pp.nwire) {
                        std::copy(charge.begin(), charge.end(),
                                  pp.q.begin() + (size_t)iw * pp.ntick);
                    }
                }
            } // wires
        }     // plane

        if (want_trios) {
            add_trios(face, patches, m_trio_min_charge,
                      trios, ntrio_depos, ntrio_incomplete);
        }
    }         // depos

    out = accum->frame(in->ident(), m_tbins.min() - m_reftime, m_tbins.binsize());
    log->debug("splat {} ndepos={}/{}/[{}] ntraces={}",
               out->ident(), ndepos_seen, in->depos()->size(), nplanes_skipped, accum->ntraces());

    if (want_trios) {
        // Sorted so the file is reproducible; the accumulator is unordered for
        // the ~10^7 insertions, not for the one pass taken here.
        std::vector<const TrioMap::value_type*> rows;
        rows.reserve(trios.size());
        for (const auto& kv : trios) rows.push_back(&kv);
        std::sort(rows.begin(), rows.end(),
                  [](auto a, auto b) { return a->first < b->first; });

        // int16 holds a tick or a channel for any plausible readout, and these
        // arrays are the bulk of the output, so the narrower type is worth it --
        // but silently wrapping would corrupt the association, so check.
        const int kmax = std::numeric_limits<int16_t>::max();
        for (const auto* r : rows) {
            for (int icol = 0; icol < 4; ++icol) {
                if (r->first[icol] < 0 || r->first[icol] > kmax) {
                    THROW(ValueError() << errmsg{
                        "DepoFluxSplat: trio index " + std::to_string(r->first[icol])
                        + " does not fit in int16; readout too long or too many channels"});
                }
            }
        }

        // A given (u,v,w) is typically hot for a RUN of consecutive ticks --
        // measured at ~24 of them, since that is the depo's diffused time
        // extent.  Emitting one row per tick repeats the three channels 24
        // times over.  Run-length encoding the tick axis removes exactly that
        // redundancy and nothing else: the charges stay one per tick, in run
        // order, so the file still says what it said before.
        std::sort(rows.begin(), rows.end(), [](auto a, auto b) {
            const auto& x = a->first;
            const auto& y = b->first;
            if (x[1] != y[1]) return x[1] < y[1];   // u
            if (x[2] != y[2]) return x[2] < y[2];   // v
            if (x[3] != y[3]) return x[3] < y[3];   // w
            return x[0] < y[0];                     // then tick
        });

        std::vector<size_t> run_start;
        for (size_t i = 0; i < rows.size(); ++i) {
            bool fresh = (i == 0);
            if (!fresh) {
                const auto& p = rows[i-1]->first;
                const auto& c = rows[i]->first;
                // A run breaks on a different triple or a gap in tick.  Not
                // every triple is one run: ~7% are interrupted, and those
                // simply emit another row.
                fresh = (c[1] != p[1] || c[2] != p[2] || c[3] != p[3]
                         || c[0] != p[0] + 1);
            }
            if (fresh) run_start.push_back(i);
        }

        Array::array_xxs runs(run_start.size(), 5);
        Array::array_xxf value(rows.size(), 1);
        for (size_t r = 0; r < run_start.size(); ++r) {
            const size_t b = run_start[r];
            const size_t e = (r+1 < run_start.size()) ? run_start[r+1] : rows.size();
            const auto& k = rows[b]->first;
            if ((e - b) > (size_t)kmax) {
                THROW(ValueError() << errmsg{"DepoFluxSplat: trio run too long for int16"});
            }
            runs(r, 0) = (short)k[1];        // chan U
            runs(r, 1) = (short)k[2];        // chan V
            runs(r, 2) = (short)k[3];        // chan W
            runs(r, 3) = (short)k[0];        // first nominal tick
            runs(r, 4) = (short)(e - b);     // number of ticks
        }
        for (size_t i = 0; i < rows.size(); ++i) {
            value(i, 0) = (float)rows[i]->second;
        }

        // Only the very first write may truncate.  Everything after must
        // append, or it truncates away the array just written.
        const std::string sid = std::to_string(out->ident());
        Numpy::save2d(runs, "trio_runs_" + sid, m_trio_file,
                      m_count ? "a" : "w");
        Numpy::save2d(value, "trio_value_" + sid, m_trio_file, "a");

        if (!m_count) {
            // Written once: the reader adds these to the nominal tick to get
            // each plane's own tick.  Zero unless time_offsets is configured.
            Array::array_xxs offs(1, 3);
            for (int i = 0; i < 3; ++i) offs(0, i) = (short)m_tick_offsets[i];
            Numpy::save2d(offs, "trio_tick_offsets", m_trio_file, "a");
        }
        log->debug("splat {} trio runs={} depos={}/[{}] -> {}",
                   out->ident(), run_start.size(), ntrio_depos, ntrio_incomplete,
                   m_trio_file);
    }

    ++m_count;

    // (void)ndepos_skipped;
    return true;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
