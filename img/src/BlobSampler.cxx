#include "WireCellImg/BlobSampler.h"

#include "WireCellUtil/Range.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

#include <iostream>             // debug

WIRECELL_FACTORY(BlobSampler, WireCell::Img::BlobSampler,
                 WireCell::INamed,
                 WireCell::IBlobSampler,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Img;
using namespace WireCell::Range;
using namespace WireCell::String;
using namespace WireCell::RayGrid;
using namespace WireCell::PointCloud;

BlobSampler::BlobSampler()
    : Aux::Logger("BlobSampler", "img")
{

}


Configuration cpp2cfg(const BlobSampler::CommonConfig& cc)
{
    Configuration cfg;
    cfg["time_offset"] = cc.time_offset;
    cfg["drift_speed"] = cc.drift_speed;
    cfg["prefix"] = cc.prefix;
    cfg["tbins"] = cc.tbinning.nbins();
    cfg["tmin"] = cc.tbinning.min();
    cfg["tmax"] = cc.tbinning.max();
    cfg["extra"] = Json::arrayValue;
    for (const auto& res : cc.extra) {
        cfg["extra"].append(res);
    }
    return cfg;
}
BlobSampler::CommonConfig cfg2cpp(const Configuration& cfg,
                                  const BlobSampler::CommonConfig& def = {})
{
    BlobSampler::CommonConfig cc;
    cc.time_offset = get(cfg, "time_offset", def.time_offset);
    cc.drift_speed = get(cfg, "drift_speed", def.drift_speed);
    cc.prefix = get(cfg, "prefix", def.prefix);
    cc.tbinning = Binning(
        get(cfg, "tbins", def.tbinning.nbins()),
        get(cfg, "tmin", def.tbinning.min()),
        get(cfg, "tmax", def.tbinning.max()));
    cc.extra = get(cfg, "extra", cc.extra);
    cc.extra_re.clear();
    for (const auto& res : get(cfg, "extra", def.extra)) {
        cc.extra_re.push_back(std::regex(res));
    }

    return cc;
}

WireCell::Configuration BlobSampler::default_configuration() const
{
    auto cfg = cpp2cfg(m_cc);
    cfg["strategy"] = "center";
    return cfg;
}

void BlobSampler::configure(const WireCell::Configuration& cfg)
{
    m_cc = cfg2cpp(cfg, m_cc);

    add_strategy(cfg["strategy"]);
}


// Local base class providing common functionality to all samplers.
//
// Subclass should provide sample() which should call intern() with
// points.
struct BlobSampler::Sampler : public Aux::Logger
{
    // Little helper to reduct setting the same value over n points and
    // only doing so if the suffix is configured.
    struct npts_dup {
        Sampler& samp;
        Dataset& ds;
        size_t npts;

        template<typename Num>
        void operator()(const std::string& suffix, Num val)
        {
            if (samp.is_extra(suffix)) {
                std::vector<Num> vals(npts, val);
                ds.add(samp.cc.prefix+suffix, Array(vals));
            }
        }
    };
    struct pts_vec {
        Sampler& samp;
        Dataset& ds;
        std::string letter;

        template<typename Num>
        void operator()(const std::string& suffix, const std::vector<Num>& vals)
        {
            std::string lsuffix = letter + suffix;
            if (samp.is_extra(lsuffix)) {
                ds.add(samp.cc.prefix+lsuffix, Array(vals));
            }
        }
    };

    // Hard-wired, common subset of full config.
    BlobSampler::CommonConfig cc;

    // The identity of this sampler
    size_t my_ident;

    // Current blob index in the iterated IBlob::vector
    size_t blob_index;
    IBlob::pointer iblob;
    // size_t points_added{0};

    explicit Sampler(const Configuration& cfg, size_t ident)
        : Aux::Logger("BlobSampler", "img")
        , cc(cfg2cpp(cfg)), my_ident(ident)
    {
        // this->configure(cfg); can not call virtual methods from ctro!  Outer
        // context must call our configure() after we have been constructed.
    }

    virtual ~Sampler() {}

    // Context manager around sample()
    IAnodeFace::pointer anodeface;
    void begin_sample(size_t bind, IBlob::pointer fresh_iblob)
    {
        // points_added = 0;
        blob_index = bind;
        iblob = fresh_iblob;
        anodeface = fresh_iblob->face();
    }
    void end_sample()
    {
        // points_added=0;
        anodeface = nullptr;
        blob_index=0;
        iblob = nullptr;
    }

    // Entry point to subclasses
    virtual void sample(Dataset& ds) = 0;

    // subclass may want to config self.
    virtual void configure(const Configuration& cfg) { }

    // Return pimpos for a given plane index
    const Pimpos* pimpos(int plane_index=2) const
    {
        return anodeface->planes()[plane_index]->pimpos();
    }
    double plane_x(int plane_index=2) const
    {
        return anodeface->planes()[plane_index]->wires().front()->center().x();
    }

    // Convert a signal time to its location along global drift
    // coordinates.
    double time2drift(double time) const
    {
        const Pimpos* colpimpos = pimpos(2);
        const double drift = (time + cc.time_offset)*cc.drift_speed;
        double xorig = plane_x(2); // colpimpos->origin()[0];
        /// TODO: how to determine xsign?
        double xsign = colpimpos->axis(0)[0];
        return xorig + xsign*drift;
    }

    // Return a wire crossing as a 2D point held in a 3D point.  The
    // the drift coordinate is unset.  Use time2drift().
    Point crossing_point(const RayGrid::crossing_t& crossing)
    {
        const auto& coords = anodeface->raygrid();
        const auto& [one, two] = crossing;
        auto pt = coords.ray_crossing(one, two);
        return pt;
    }

    Point center_point(const crossings_t& corners)
    {
        Point avg;
        for (const auto& corner : corners) {
            avg += crossing_point(corner);
        }
        avg = avg / corners.size();
        return avg;
    }

    // Match array name suffix against regexp
    bool is_extra(const std::string& suffix) {
        std::smatch smatch;
        for (const auto& re : cc.extra_re) {
            if (std::regex_match(suffix, smatch, re)) {
                return true;
            }
        }
        return false;
    }
    
    // Fixme: this cache is not thread safe if a BlobSampler is shared
    // to multiple BlobSamplings.  To fix, have BlobSampler be
    // configured for IAnodePlanes and pre-fill the lookup.
    struct ChanInfo {
        int ident, index;
    };
    using ident2index_t = std::unordered_map<int, int>;
    using plane_ident2index_t = std::unordered_map<IWirePlane::pointer, ident2index_t>;
    plane_ident2index_t plane_ident2index;

    bool want_extra(const std::string& letter, const std::vector<std::string>& names) {
        for (const auto& name : names) {
            if (is_extra(letter + name)) { return true; }
        }
        return false;
    }

    // Return a dataset covering multiple points related to a blob
    Dataset make_dataset(const std::vector<Point>& pts, double time)
    {
        size_t npts = pts.size();

        std::vector<float> times(npts, time);
        std::vector<Point::coordinate_t> x(npts),y(npts),z(npts);

        for (size_t ind=0; ind<npts; ++ind) {
            const auto& pt = pts[ind];
            x[ind] = pt.x();
            y[ind] = pt.y();
            z[ind] = pt.z();
        }

        Dataset ds({
                {cc.prefix + "x", Array(x)},
                {cc.prefix + "y", Array(y)},
                {cc.prefix + "z", Array(z)},
                {cc.prefix + "t", Array(times)}});

        // Extra values
        const std::vector<std::string> uvw = {"u","v","w"};
        auto islice = iblob->slice();

        // Per blob, duplicated over all pts.
        {
            npts_dup nd{*this, ds, npts};
            nd("sample_strategy", my_ident);
            nd("blob_ident", iblob->ident());
            nd("blob_index", blob_index);
            nd("slice_ident", islice->ident());
            nd("slice_start", islice->start());
            nd("slice_span", islice->span());
        }        

        if (cc.extra_re.empty()) {
            return ds;
        }

        // Per point arrays

        const auto& activity = islice->activity();
        auto iface = iblob->face();
        for (const auto& iplane : iface->planes()) {

            const auto* pimpos = iplane->pimpos();
            const int pind = iplane->planeid().index();
            const std::string letter = uvw[pind];

            // Not sure if pre-checking all the regex is faster than
            // simply calculating everything first and checking
            // one-by-one at nv/add() time....
            if (! want_extra(letter, {"wire_index", "channel_ident", "channel_attach",
                                      "pitch_coord", "wire_coord", "charge_val", "charge_unc"})) {
                continue;
            }

            pts_vec nv{*this, ds, letter};
            
            const IChannel::vector& channels = iplane->channels();

            // Hit cache for ch ident -> ch index
            auto& p_chi2i = plane_ident2index[iplane];
            if (p_chi2i.empty()) {
                const size_t nchannels = channels.size();
                for (size_t chind=0; chind<nchannels; ++chind) {
                    auto ich = channels[chind];
                    p_chi2i[ich->ident()] = chind;
                }
            }

            std::vector<int> wire_index(npts, -1), channel_ident(npts, -1), channel_attach(npts, -1);
            std::vector<double> pitch_coord(npts,0), wire_coord(npts,0), charge_val(npts,0), charge_unc(npts,0);

            const IWire::vector& iwires = iplane->wires();

            for (size_t ipt=0; ipt<npts; ++ipt) {
                const Point xwp = pimpos->transform(pts[ipt]);
                wire_coord[ipt] = xwp[1];
                const double pitch = xwp[2];
                pitch_coord[ipt] = pitch;
                int wind = pimpos->closest(pitch).first;
                if (wind < 0) {
                    log->debug("sampler={}, point={} cartesian={} pimpos={}", my_ident, ipt, pts[ipt], xwp);
                    log->error("Negative wire index: {}, will segfault soon", wind);

                }
                wire_index[ipt] = wind;
        
                IWire::pointer iwire = iwires[wire_index[ipt]];
                channel_ident[ipt] = iwire->channel();
                channel_attach[ipt] = p_chi2i[channel_ident[ipt]];
                auto ich = channels[channel_attach[ipt]];

                auto ait = activity.find(ich);
                if (ait != activity.end()) {
                    auto act = ait->second;
                    charge_val[ipt] = act.value();
                    charge_unc[ipt] = act.uncertainty();
                }
            }

            nv("wire_index", wire_index);
            nv("pitch_coord", pitch_coord);
            nv("wire_coord", wire_coord);
            nv("channel_ident", channel_ident);
            nv("channel_attach", channel_attach);
            nv("charge_val", charge_val);
            nv("charge_unc", charge_unc);

        } // over planes

        return ds;
    }

    // Append points to PC.  
    void intern(Dataset& ds,
                std::vector<Point> points)
    {
        auto islice = iblob->slice();
        const double t0 = islice->start();
        const double dt = islice->span();

        const Binning bins(cc.tbinning.nbins(),
                           cc.tbinning.min()*dt + t0,
                           cc.tbinning.max()*dt + t0);
        // log->debug("t0 {} dt {} bins {}", t0, dt, bins);
        const size_t npts = points.size();
        for (int tbin : irange(bins.nbins())) {
            const double time = bins.edge(tbin);
            const double x = time2drift(time);
            // points_added += npts;
            for (size_t ind=0; ind<npts; ++ind) {
                points[ind].x(x);
            }
            auto tail = make_dataset(points, time);
            const size_t before = ds.size_major();
            ds.append(tail);
            const size_t after = ds.size_major();
            // log->debug("sampler {} iblob {} intern {} points, ds size {}, tail size {} with binning {}",
            //            ident, iblob->ident(), npts, ds.size_major(), tail.size_major(), bins);
            if (after != before + npts) {
                THROW(LogicError() << errmsg{"PointCloud append() is broken"});
            }
        }
        
    }
};

PointCloud::Dataset BlobSampler::sample_blob(const IBlob::pointer& iblob,
                                             int blob_index)
{
    if (!iblob) {
        THROW(ValueError() << errmsg{"can not sample null blob"});
    }

    PointCloud::Dataset ret;
    // size_t points_added = 0;

    for (auto& sampler : m_samplers) {
        sampler->begin_sample(blob_index, iblob);
        sampler->sample(ret);
        // points_added += sampler->points_added;
        sampler->end_sample();
    }
    // log->debug("got {} blobs, sampled {} points with {} samplers, returning {}",
    //            nblobs, points_added, m_samplers.size(), ret.size_major());
    return ret;
}

// PointCloud::Dataset BlobSampler::sample_blobs(const IBlob::vector& iblobs)
// {
//     PointCloud::Dataset ret;
//     size_t nblobs = iblobs.size();
//     size_t points_added = 0;
//     for (size_t bind=0; bind<nblobs; ++bind) {
//         auto fresh_iblob = iblobs[bind];
//         for (auto& sampler : m_samplers) {
//             if (!fresh_iblob) {
//                 THROW(ValueError() << errmsg{"can not sample null blob"});
//             }
//             sampler->begin_sample(bind, fresh_iblob);
//             sampler->sample(ret);
//             points_added += sampler->points_added;
//             sampler->end_sample();
//         }
//     }
//     log->debug("got {} blobs, sampled {} points with {} samplers, returning {}",
//                nblobs, points_added, m_samplers.size(), ret.size_major());
//     return ret;
// }
        

struct Center : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Center(const Center&) = default;
    Center& operator=(const Center&) = default;

    void sample(Dataset& ds)
    {
        auto corners = iblob->shape().corners();
        std::vector<Point> points(1);
        points[0] = center_point(corners);
        intern(ds, points);
    }
};


struct Corner : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Corner(const Corner&) = default;
    Corner& operator=(const Corner&) = default;
    int span{-1};
    virtual void configure(const Configuration& cfg)
    {
        span = get(cfg, "span", span);
    }
    void sample(Dataset& ds)
    {
        const auto& corners = iblob->shape().corners();
        const size_t npts = corners.size();
        if (! npts) return;
        std::vector<Point> points;
        points.reserve(npts);
        for (const auto& corner : corners) {
            points.emplace_back(crossing_point(corner));
        }
        intern(ds, points);
    }
};


struct Edge : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Edge(const Edge&) = default;
    Edge& operator=(const Edge&) = default;

    void sample(Dataset& ds)
    {
        const auto& coords = anodeface->raygrid();
        auto pts = coords.ring_points(iblob->shape().corners());
        const size_t npts = pts.size();

        // walk around the ring of points, find midpoint of each edge.
        std::vector<Point> points;
        for (size_t ind1=0; ind1<npts; ++ind1) {
            size_t ind2 = (1+ind1)%npts;
            const auto& origin = pts[ind1];
            const auto& egress = pts[ind2];
            const auto mid = 0.5*(egress + origin);
            points.push_back(mid);
        }
        intern(ds, points);
    }
};


struct Grid : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Grid(const Grid&) = default;
    Grid& operator=(const Grid&) = default;
    
    double step{1.0};
    // default planes to use
    std::vector<size_t> planes = {0,1,2}; 
    
    virtual void configure(const Configuration& cfg)
    {
        step = get(cfg, "step", step);
        planes = get(cfg, "planes", planes);
        if (planes.size() != 2) {
            THROW(ValueError() << errmsg{"illegal size Grid.planes"});
        }
        size_t tot = planes[0] + planes[1];
        //                                  x,0+1,0+2, 1+2
        const std::vector<size_t> other = {42,  2,  1, 0};
        planes.push_back(other[tot]);
    }

    void sample(Dataset& ds)
    {
        if (step == 1.0) {
            aligned(ds, iblob);
        }
        else {
            unaligned(ds, iblob);
        }
    }

    // Special case where points are aligned to grid
    void aligned(Dataset& ds, IBlob::pointer iblob)
    {
        const auto& coords = anodeface->raygrid();
        const auto& strips = iblob->shape().strips();

        // chosen layers
        const layer_index_t l1 = 2 + planes[0];
        const layer_index_t l2 = 2 + planes[1];
        const layer_index_t l3 = 2 + planes[2];

        // strips
        const auto& s1 = strips[l1];
        const auto& s2 = strips[l2];
        const auto& s3 = strips[l3];

        std::vector<Point> points;
        for (auto gi1 : irange(s1.bounds)) {
            coordinate_t c1{l1, gi1};
            for (auto gi2 : irange(s2.bounds)) {
                coordinate_t c2{l2, gi2};
                const double pitch = coords.pitch_location(c1, c2, l3);
                auto gi3 = coords.pitch_index(pitch, l3);
                if (s3.in(gi3)) {
                    auto pt = coords.ray_crossing(c1, c2);
                    points.push_back(pt);
                }
            }
        }
        intern(ds, points);
    }

    void unaligned(Dataset& ds, IBlob::pointer iblob)
    {
        std::vector<Point> pts;
        const auto& strips = iblob->shape().strips();

        // chosen layers
        const layer_index_t l1 = 2 + planes[0];
        const layer_index_t l2 = 2 + planes[1];
        const layer_index_t l3 = 2 + planes[2];

        const auto* pimpos3 = pimpos(planes[2]);

        const auto& coords = anodeface->raygrid();

        // fixme: the following code would be nice to transition to a
        // visitor pattern in the RayGrid::Coordinates class.

        // A jump along one axis between neighboring ray crossings from another axis, 
               // expressed as relative 3-vectors in Cartesian space.
               const auto& jumps = coords.ray_jumps();

        const auto& j1 = jumps(l1,l2);
        const double m1 = j1.magnitude();
        const auto& u1 = j1/m1; // unit vector

        const auto& j2 = jumps(l2,l1);
        const double m2 = j2.magnitude();
        const auto& u2 = j2/m2; // unit vector

        // The strips bound the blob in terms of ray indices.
        const auto& s1 = strips[l1];
        const auto& s2 = strips[l2];
        const auto& s3 = strips[l3];

        // Work out the index spans of first two planes
        auto b1 = s1.bounds.first;
        auto e1 = s1.bounds.second;
        if (e1 < b1) std::swap(b1,e1);
        auto b2 = s2.bounds.first;
        auto e2 = s2.bounds.second;
        if (e2 < b2) std::swap(b2,e2);

        // The maximum distance the blob spans along one direction is
        // that direction's jump size times number crossing rays from
        // the strip in the other direction.
        const int max1 = m1*(e2-b2);
        const int max2 = m2*(e1-b1);
        
        // Physical lowest point of two crossing strips.
        const auto origin = coords.ray_crossing({l1,b1}, {l2, b2});

        // Iterate by jumping along the direction of each of the two
        // planes and check if the result is inside the third strip.
        std::vector<Point> points;
        for (double step1 = m1; step1 < max1; step1 += m1) {
            const auto pt1 = origin + u1 * step1;
            for (double step2 = m2; step2 < max2; step2 += m2) {
                const auto pt2 = pt1 + u2 * step2;

                const double pitch = pimpos3->distance(pt2);
                const int pi = coords.pitch_index(pitch, l3);
                if (s3.in(pi)) {
                    points.push_back(pt2);
                }
            }
        }
        intern(ds, points);
    }
};


struct Bounds : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Bounds(const Bounds&) = default;
    Bounds& operator=(const Bounds&) = default;

    double step{1.0};
    virtual void configure(const Configuration& cfg)
    {
        step = get(cfg, "step", step);
    }

    void sample(Dataset& ds)
    {
        const auto& coords = anodeface->raygrid();
        auto pts = coords.ring_points(iblob->shape().corners());
        const size_t npts = pts.size();

        // walk around the ring of points, taking first step from a
        // corner until we step off the boundary.
        std::vector<Point> points;
        for (size_t ind1=0; ind1<npts; ++ind1) {
            size_t ind2 = (1+ind1)%npts;
            const auto& origin = pts[ind1];
            const auto& egress = pts[ind2];
            const auto diff = egress-origin;
            const double mag = diff.magnitude();
            const auto unit = diff/mag;
            for (double loc = step; loc < mag; loc += step) {
                auto pt = origin + unit*loc;
                points.push_back(pt);
            }
        }
        intern(ds, points);
    }
};


// Implement the "stepped" sampling.
//
// Outline:
//
// 1. Find N_{1,2} wires for plain p_{1,2} with {minimum,maximum} number of wires in blob.
// 2. Find S_{1,2} = max(3, N_{1,2}/12)
// 3. Accept points on sub-grid steps (S_1, S_2)
//
// Note, combine with "bounds" strategy to fully reproduce the
// sampling used in the WC prototype imaging.
struct Stepped : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Stepped(const Stepped&) = default;
    Stepped& operator=(const Stepped&) = default;
    virtual ~Stepped() {}

    // The minimium number of wires over which a step will be made.
    double min_step_size{3};
    // The maximum fraction of a blob a step may take.  If
    // non-positive, then all steps are min_step_size.
    double max_step_fraction{1.0/12.0};

    // Distance along each ray from a crossing to place a point.  Value is in
    // units of pitch.  Default value of 0.5 picks points at the wire crossing
    // point instead of the ray crossing point.
    double offset{0.5};

    virtual void configure(const Configuration& cfg)
    {
        min_step_size = get(cfg, "min_step_size", min_step_size);
        max_step_fraction = get(cfg, "max_step_fraction", max_step_fraction);
        offset = get(cfg, "offset", offset);
    }


    void sample(Dataset& ds) {
        const auto& coords = anodeface->raygrid();
        auto strips = iblob->shape().strips();

        auto swidth = [](const Strip& s) -> int {
            return s.bounds.second - s.bounds.first;
        };
        std::sort(strips.begin()+2, strips.end(),
                  [&](const Strip& a, const Strip& b) -> bool {
                      return swidth(a) < swidth(b);
                  });
        const Strip& smin = strips[2];
        const Strip& smid = strips[3];
        const Strip& smax = strips[4];
        
        int nmin = std::max(min_step_size, max_step_fraction*swidth(smin));
        int nmax = std::max(min_step_size, max_step_fraction*swidth(smax));

        std::vector<Point> points;

        const Vector adjust = offset * (
            coords.ray_crossing({smin.layer, 0}, {smax.layer, 0}) -
            coords.ray_crossing({smin.layer, 1}, {smax.layer, 1}));

        // log->debug("offset={} adjust={},{},{}", offset, adjust.x(), adjust.y(), adjust.z());

        for (auto gmin=smin.bounds.first; gmin < smin.bounds.second; gmin += nmin) {
            coordinate_t cmin{smin.layer, gmin};
            for (auto gmax=smax.bounds.first; gmax < smax.bounds.second; gmax += nmax) {
                coordinate_t cmax{smax.layer, gmax};
                
                const double pitch = coords.pitch_location(cmin, cmax, smid.layer);
                auto gmid = coords.pitch_index(pitch, smid.layer);
                if (smid.in(gmid)) {
                    const auto pt = coords.ray_crossing(cmin, cmax);
                    points.push_back(pt + adjust);
                }
            }
        }
        intern(ds, points);
    }
};

void BlobSampler::add_strategy(Configuration strategy)
{
    if (strategy.isNull()) {
        strategy = cpp2cfg(m_cc);
        strategy["name"] = "center";
        add_strategy(strategy);
        return;
    }
    if (strategy.isString()) {
        std::string name = strategy.asString();
        strategy = cpp2cfg(m_cc);
        strategy["name"] = name;
        add_strategy(strategy);
        return;
    }
    if (strategy.isArray()) {
        for (auto one : strategy) {
            add_strategy(one);
        }
        return;
    }
    if (! strategy.isObject()) {
        THROW(ValueError() << errmsg{"unsupported strategy type"});
    }

    auto full = cpp2cfg(m_cc);
    for (auto key : strategy.getMemberNames()) {
        full[key] = strategy[key];
    }
    // log->debug("making strategy: {}", full);

    std::string name = full["name"].asString();
    // use startswith() to be a little friendly to the user
    // w.r.t. spelling.  eg "centers", "bounds", "boundaries" are all
    // accepted.
    if (startswith(name, "center")) {
        m_samplers.push_back(std::make_unique<Center>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }
    if (startswith(name, "corner")) {
        m_samplers.push_back(std::make_unique<Corner>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }
    if (startswith(name, "edge")) {
        m_samplers.push_back(std::make_unique<Edge>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }
    if (startswith(name, "grid")) {
        m_samplers.push_back(std::make_unique<Grid>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }
    if (startswith(name, "bound")) {
        m_samplers.push_back(std::make_unique<Bounds>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }
    if (startswith(name, "stepped")) {
        m_samplers.push_back(std::make_unique<Stepped>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }

    THROW(ValueError() << errmsg{"unknown strategy: " + name});
}

        
