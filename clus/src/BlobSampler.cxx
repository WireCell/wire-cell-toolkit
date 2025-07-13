#include "WireCellClus/BlobSampler.h"

#include "WireCellUtil/Range.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

#include <iostream>             // debug

WIRECELL_FACTORY(BlobSampler, WireCell::Clus::BlobSampler,
                 WireCell::INamed,
                 WireCell::IBlobSampler,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Clus;
using namespace WireCell::Range;
using namespace WireCell::String;
using namespace WireCell::RayGrid;
using namespace WireCell::PointCloud;

BlobSampler::BlobSampler()
    : Aux::Logger("BlobSampler", "clus")
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
                // std::cout << "test: " << samp.cc.prefix << " " << lsuffix << " " << vals.size() << std::endl;
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
        : Aux::Logger("BlobSampler", "clus")
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
    virtual void sample(Dataset& ds, Dataset& aux) = 0;

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
        // const Pimpos* colpimpos = pimpos(2);
        const double drift = (time + cc.time_offset)*cc.drift_speed;
        double xorig = plane_x(2); // colpimpos->origin()[0];
        /// TODO: how to determine xsign?
        double xsign = anodeface->dirx();
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
        WirePlaneId wpid_blob{0}; // 0 is invalid, assign when we get it. duplicate over all pts later.
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

                // auto temp = pimpos->closest(pitch);
                int wind =  pimpos->closest(pitch + 0.1*units::mm).first; // shift to the higher wires in case of a tie ... 
                // std::cout << temp.second << " " << wind << " " << temp1.second << " " << temp1.first << std::endl;
                
                if (wind < 0) {
                    log->debug("sampler={}, point={} cartesian={} pimpos={}", my_ident, ipt, pts[ipt], xwp);
                    log->error("Negative wire index: {}, will segfault soon", wind);

                }
                wire_index[ipt] = wind;
        
                IWire::pointer iwire = iwires[wire_index[ipt]];
                const auto& wpid_wire = iwire->planeid();
                wpid_blob = WireCell::WirePlaneId(kAllLayers, wpid_wire.face(), wpid_wire.apa());
                channel_ident[ipt] = iwire->channel();
                channel_attach[ipt] = p_chi2i[channel_ident[ipt]];
                auto ich = channels[channel_attach[ipt]];

                

                auto ait = activity.find(ich);
                if (ait != activity.end()) {
                    auto act = ait->second;
                    charge_val[ipt] = act.value();
                    charge_unc[ipt] = act.uncertainty();
                }

                // std::cout << "Test: wire_index " << wire_index[ipt]
                //           << " pitch_coord " << pitch_coord[ipt]
                //           << " wire_coord " << wire_coord[ipt]
                //           << " channel_ident " << channel_ident[ipt]
                //           << " channel_attach " << channel_attach[ipt]
                //           << " charge_val " << charge_val[ipt]
                //           << " charge_unc " << charge_unc[ipt] 
                //           << std::endl;
            }

            nv("wire_index", wire_index);
            nv("pitch_coord", pitch_coord);
            nv("wire_coord", wire_coord);
            nv("channel_ident", channel_ident);
            nv("channel_attach", channel_attach);
            nv("charge_val", charge_val);
            nv("charge_unc", charge_unc);

        } // over planes

        {
            npts_dup nd{*this, ds, npts};
            nd("wpid", wpid_blob.ident());
        }        

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

std::tuple<PointCloud::Dataset, PointCloud::Dataset> BlobSampler::sample_blob(const IBlob::pointer& iblob,
                                             int blob_index)
{
    PointCloud::Dataset ret_main;
    PointCloud::Dataset ret_aux;
    // size_t points_added = 0;

    for (auto& sampler : m_samplers) {
        sampler->begin_sample(blob_index, iblob);
        sampler->sample(ret_main, ret_aux);
        // points_added += sampler->points_added;
        sampler->end_sample();
    }
    // log->debug("got {} blobs, sampled {} points with {} samplers, returning {}",
    //            nblobs, points_added, m_samplers.size(), ret_main.size_major());
    return {ret_main, ret_aux};
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

    void sample(Dataset& ds, Dataset& aux)
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
    void sample(Dataset& ds, Dataset& aux)
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

    void sample(Dataset& ds, Dataset& aux)
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
    std::vector<size_t> planes = {0,1}; 
    
    virtual void configure(const Configuration& cfg)
    {
        step = get(cfg, "step", step);
        planes = get(cfg, "planes", planes);
        if (planes.size() != 2) {
            raise<ValueError>("illegal size for Grid.planes: %d", planes.size());
        }
        size_t tot = planes[0] + planes[1];
        //                                  x,0+1,0+2, 1+2
        const std::vector<size_t> other = {42,  2,  1, 0};
        planes.push_back(other[tot]);
    }

    void sample(Dataset& ds, Dataset& aux)
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

    void sample(Dataset& ds, Dataset& aux)
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

    double tolerance{0.03};

    virtual void configure(const Configuration& cfg)
    {
        min_step_size = get(cfg, "min_step_size", min_step_size);
        max_step_fraction = get(cfg, "max_step_fraction", max_step_fraction);
        offset = get(cfg, "offset", offset);
    }


    void sample(Dataset& ds, Dataset& aux) {
        const auto& coords = anodeface->raygrid();
        auto strips = iblob->shape().strips();
        const int ndummy_index = strips.size() == 5 ? 2 : 0; // use this to skip dummy planes
        const int li[3] = {ndummy_index+0, ndummy_index+1, ndummy_index+2}; // layer index
        // std::cout << "DEBUG strips.size() " << strips.size() << std::endl;
        // for (const auto& strip : strips) {
        //     std::cout << "DEBUG strip " << strip.layer << " " << strip.bounds.first << " " << strip.bounds.second << std::endl;
        // }

        // This returns the number of wire regions covered by the strip s.
        auto swidth = [](const Strip& s) -> int {
            return s.bounds.second - s.bounds.first;
        };

        // Find the strip with largest coverage.
        Strip smax = strips[li[0]]; int max_id = li[0];
        if (swidth(strips[li[1]]) > swidth(smax)){
            smax = strips[li[1]]; max_id = li[1];
        }
        if(swidth(strips[li[2]]) > swidth(smax)){
            smax = strips[li[2]]; max_id = li[2];
        }

        // Find the strip with least coverage.
        Strip smin = strips[li[1]]; int min_id = li[1];
        if (swidth(strips[li[0]]) < swidth(smin)){
            smin = strips[li[0]]; min_id = li[0];
        }
        if(swidth(strips[li[2]]) < swidth(smin)){
            smin = strips[li[2]]; min_id = li[2];
        }

        // Find the other strip.
        Strip smid = strips[li[2]]; /*int mid_id = li[2];*/
        for (int i = li[0];i!=li[2]+1;i++){
            if (i != max_id && i != min_id){
                smid = strips[i];        
                // mid_id = i;
            }
        }
        
        // Step sizes for the min/max directions.
        int nmin = std::max(min_step_size, max_step_fraction*swidth(smin));
        int nmax = std::max(min_step_size, max_step_fraction*swidth(smax));

        std::vector<Point> points;

        //XQ: is the order of 0 vs. 1 correct for the wire center???
        //
        // BV: this gives a diagonal vector from the crossing point of the
        // 0-rays to the crossing point of the 0-wires (0-wire = half-way from
        // 0-ray to 1-ray, half-way assuming offset is default 0.5).
        //
        //        /(0-ray) 
        //       / /(0-wire)  
        //      +-+-+ b
        //     /   /
        //    / c +-------(0-wire)
        //   /   /
        //  +---+----(0-ray)
        //  a
        // 
        // if "a" is the crossing of two 0-rays and "b" is the crossing of two
        // 1-rays, "adjust is the vector "ac" which is coincident with the
        // crossing of the 0-wires.
        const Vector adjust = offset * (
            coords.ray_crossing({smin.layer, 1}, {smax.layer, 1}) -
            coords.ray_crossing({smin.layer, 0}, {smax.layer, 0}));

        // This gives a relative pitch distance measured in the "mid" view that
        // is half the distance between crossing point of the 0-rays and the
        // 1-rays in the other two views.  In general, this is NOT the same as
        // the magnitude of "adjust" / "ac" vector above as that diagonal of the
        // min/max parallelogram is not necessarily parallel to the pitch
        // direction in the third, "mid" view.  The two directions are
        // accidentally coincident for symmetric wire patterns like in
        // MicroBooNE.
        const double pitch_adjust = offset * (
            coords.pitch_location({smin.layer, 1}, {smax.layer, 1}, smid.layer) -
            coords.pitch_location({smin.layer, 0}, {smax.layer, 0}, smid.layer) ); 

        // log->debug("offset={} adjust={},{},{}", offset, adjust.x(), adjust.y(), adjust.z());

        std::set<decltype(smin.bounds.first)> min_wires_set;
        std::set<decltype(smin.bounds.first)> max_wires_set;

        // Load up wires for the min/max dimensions including first, along each
        // step size and last wire.
        for (auto gmin=smin.bounds.first; gmin < smin.bounds.second; gmin += nmin) { 
            min_wires_set.insert(gmin);
        }
        min_wires_set.insert(smin.bounds.second-1);
        for (auto gmax=smax.bounds.first; gmax < smax.bounds.second; gmax += nmax) {
            max_wires_set.insert(gmax);
        }
        max_wires_set.insert(smax.bounds.second-1);

        // size_t npre_missed=0;
        // size_t nrel_missed=0;
        for (auto it_gmin = min_wires_set.begin(); it_gmin != min_wires_set.end(); it_gmin++){
            coordinate_t cmin{smin.layer, *it_gmin};
            for (auto it_gmax = max_wires_set.begin(); it_gmax != max_wires_set.end(); it_gmax++){
                coordinate_t cmax{smax.layer, *it_gmax};

                // Added by hyu, dunno why
                const double ploc0 = coords.pitch_location(cmin, cmax, 0);
                const double prel0 = coords.pitch_relative(ploc0, 0);
                const double ploc1 = coords.pitch_location(cmin, cmax, 1);
                const double prel1 = coords.pitch_relative(ploc1, 1);
                if (prel0 > 1 or prel0 < 0 or prel1 > 1 or prel1 < 0) {
                    // ++npre_missed;
                    continue;
                }

                // This is the mid-view pitch of the crossings of the two *wires* associated with cmin/cmax *rays*.
                const double pitch = coords.pitch_location(cmin, cmax, smid.layer) + pitch_adjust;

                // This is the location from the mid view 0-ray to the point of
                // the wire-crossing measured in units of mid view pitches.
                const double pitch_relative = coords.pitch_relative(pitch, smid.layer); 

                if (pitch_relative > smid.bounds.first - tolerance && pitch_relative < smid.bounds.second + tolerance){
                    const auto pt = coords.ray_crossing(cmin, cmax);
                    points.push_back(pt + adjust);
                //     log->warn("Blob {} adding point {} prel={}, bounds=[{},{}], tol={}, off={} padj={} adj={}",
                //               iblob->ident(), pt, pitch_relative, smid.bounds.first, smid.bounds.second, tolerance, offset, pitch_adjust, adjust);
                // }
                // else {
                //     log->warn("Blob {} not adding point prel={}, bounds=[{},{}], tol={}, off={} padj={} adj={}",
                //               iblob->ident(), pitch_relative, smid.bounds.first, smid.bounds.second, tolerance, offset, pitch_adjust, adjust);
                //     ++nrel_missed;
                }
            }
        }

        // if (points.empty()) {
        //     int ident = iblob->ident();
        //     log->warn("Blob {} unsampled: minsiz={} maxsiz={} nwiresmin={} nwiresmax={} nrel={} npre={}.", 
        //               ident, nmin, nmax, min_wires_set.size(), max_wires_set.size(), nrel_missed, npre_missed);
        //     for (const auto& strip : strips) {
        //         log->warn("Blob {} strip: {}", ident, strip);
        //     }
        // }

        intern(ds, points);

        // make aux dataset
        /// TODO: hard coded for 5 planes, i.e., wire_type is id - "2"
        aux.add("max_wire_interval", Array({(int)nmax}));
        aux.add("min_wire_interval", Array({(int)nmin}));
        aux.add("max_wire_type", Array({(int)(max_id-ndummy_index)}));
        aux.add("min_wire_type", Array({(int)(min_id-ndummy_index)}));
    }
};

// ===========================================================================
// ChargeStepped Strategy Implementation for BlobSampler
// ===========================================================================

// Implement the "charge_stepped" sampling.
//
// This is an enhanced version of "stepped" sampling that includes:
// 1. Charge-based filtering of sampling points
// 2. Bad plane handling with configurable thresholds
// 3. Conditional use of all wires vs stepped sets based on wire product
// 4. Enhanced validation logic for charge requirements
// 5. Runtime configuration override capability
//
// Based on WCPPID::calc_sampling_points() from wire-cell-pid
struct ChargeStepped : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    ChargeStepped(const ChargeStepped&) = default;
    ChargeStepped& operator=(const ChargeStepped&) = default;
    virtual ~ChargeStepped() {}

    // Configuration parameters
    double min_step_size{3};
    double max_step_fraction{1.0/12.0};
    double offset{0.5};
    double tolerance{0.03};
    
    // Charge threshold parameters
    double charge_threshold_max{4000};
    double charge_threshold_min{4000};
    double charge_threshold_other{4000};
    
    // Control parameters
    int max_wire_product_threshold{2500};
    bool disable_mix_dead_cell{true};
    
    // Dead plane detection threshold (same as PointTreeBuilding default)
    double dead_threshold{1e10};

    virtual void configure(const Configuration& cfg)
    {
        min_step_size = get(cfg, "min_step_size", min_step_size);
        max_step_fraction = get(cfg, "max_step_fraction", max_step_fraction);
        offset = get(cfg, "offset", offset);
        tolerance = get(cfg, "tolerance", tolerance);
        
        charge_threshold_max = get(cfg, "charge_threshold_max", charge_threshold_max);
        charge_threshold_min = get(cfg, "charge_threshold_min", charge_threshold_min);
        charge_threshold_other = get(cfg, "charge_threshold_other", charge_threshold_other);
        
        max_wire_product_threshold = get(cfg, "max_wire_product_threshold", max_wire_product_threshold);
        disable_mix_dead_cell = get(cfg, "disable_mix_dead_cell", disable_mix_dead_cell);
        
        dead_threshold = get(cfg, "dead_threshold", dead_threshold);
    }

    // Runtime configuration override
    void apply_runtime_config(const Configuration& runtime_cfg)
    {
        if (runtime_cfg.isMember("charge_threshold_max")) {
            charge_threshold_max = runtime_cfg["charge_threshold_max"].asDouble();
        }
        if (runtime_cfg.isMember("charge_threshold_min")) {
            charge_threshold_min = runtime_cfg["charge_threshold_min"].asDouble();
        }
        if (runtime_cfg.isMember("charge_threshold_other")) {
            charge_threshold_other = runtime_cfg["charge_threshold_other"].asDouble();
        }
        if (runtime_cfg.isMember("disable_mix_dead_cell")) {
            disable_mix_dead_cell = runtime_cfg["disable_mix_dead_cell"].asBool();
        }
        if (runtime_cfg.isMember("dead_threshold")) {
            dead_threshold = runtime_cfg["dead_threshold"].asDouble();
        }
    }

    void sample(Dataset& ds, Dataset& aux) {
        sample_with_config(ds, aux, Configuration());
    }

    void sample_with_config(Dataset& ds, Dataset& aux, const Configuration& runtime_cfg) {
        // Apply runtime configuration if provided
        if (!runtime_cfg.isNull()) {
            apply_runtime_config(runtime_cfg);
        }

        const auto& coords = anodeface->raygrid();
        auto strips = iblob->shape().strips();
        const int ndummy_index = strips.size() == 5 ? 2 : 0;
        const int li[3] = {ndummy_index+0, ndummy_index+1, ndummy_index+2};

        auto swidth = [](const Strip& s) -> int {
            return s.bounds.second - s.bounds.first;
        };

        // Find strips with max, min, and middle coverage
        Strip smax = strips[li[0]]; int max_id = li[0];
        if (swidth(strips[li[1]]) > swidth(smax)){
            smax = strips[li[1]]; max_id = li[1];
        }
        if(swidth(strips[li[2]]) > swidth(smax)){
            smax = strips[li[2]]; max_id = li[2];
        }

        Strip smin = strips[li[1]]; int min_id = li[1];
        if (swidth(strips[li[0]]) < swidth(smin)){
            smin = strips[li[0]]; min_id = li[0];
        }
        if(swidth(strips[li[2]]) < swidth(smin)){
            smin = strips[li[2]]; min_id = li[2];
        }

        Strip smid = strips[li[2]]; int mid_id = li[2];
        for (int i = li[0]; i <= li[2]; i++){
            if (i != max_id && i != min_id){
                smid = strips[i];
                mid_id = i;
            }
        }

        // Step sizes for the min/max directions
        int nmin = std::max(min_step_size, max_step_fraction*swidth(smin));
        int nmax = std::max(min_step_size, max_step_fraction*swidth(smax));

        // Pre-cache activity data for charge lookup
        auto islice = iblob->slice();
        const auto& activity = islice->activity();
        auto iface = anodeface;

        // Detect bad planes dynamically based on charge uncertainty
        std::vector<bool> plane_is_bad(3, false);
        if (disable_mix_dead_cell){
            plane_is_bad[max_id - ndummy_index] = is_plane_bad(max_id, activity, iface);
            plane_is_bad[min_id - ndummy_index] = is_plane_bad(min_id, activity, iface);
            plane_is_bad[mid_id - ndummy_index] = is_plane_bad(mid_id, activity, iface);
        }

        

        // Adjust charge thresholds based on detected bad planes
        double thresh_max = plane_is_bad[max_id - ndummy_index] ? 0.0 : charge_threshold_max;
        double thresh_min = plane_is_bad[min_id - ndummy_index] ? 0.0 : charge_threshold_min;
        double thresh_other = plane_is_bad[mid_id - ndummy_index] ? 0.0 : charge_threshold_other;

        // if (!disable_mix_dead_cell)
        //     std::cout << islice->start()/islice->span()*4<< " " << (islice->start() + islice->span())/islice->span()*4 << " "
        //               << strips[2].bounds.first << " " << strips[2].bounds.second << " "
        //               << strips[3].bounds.first << " " << strips[3].bounds.second << " "
        //               << strips[4].bounds.first << " " << strips[4].bounds.second << " " 
        //               << thresh_max << " " << thresh_min << " " << thresh_other << " ";// << std::endl;

        // Create stepped wire sets (mandatory wires)
        std::set<decltype(smin.bounds.first)> min_wires_set;
        std::set<decltype(smax.bounds.first)> max_wires_set;

        for (auto gmin = smin.bounds.first; gmin < smin.bounds.second; gmin += nmin) {
            min_wires_set.insert(gmin);
        }
        min_wires_set.insert(smin.bounds.second-1);
        
        for (auto gmax = smax.bounds.first; gmax < smax.bounds.second; gmax += nmax) {
            max_wires_set.insert(gmax);
        }
        max_wires_set.insert(smax.bounds.second-1);

        // Determine which wire sets to use
        bool use_all_wires = (swidth(smax) * swidth(smin) <= max_wire_product_threshold);
        
        std::set<decltype(smin.bounds.first)> actual_min_wires;
        std::set<decltype(smax.bounds.first)> actual_max_wires;
        
        if (use_all_wires) {
            // Use all wires
            for (auto i = smin.bounds.first; i < smin.bounds.second; i++) {
                actual_min_wires.insert(i);
            }
            for (auto i = smax.bounds.first; i < smax.bounds.second; i++) {
                actual_max_wires.insert(i);
            }
        } else {
            // Use stepped sets
            actual_min_wires = min_wires_set;
            actual_max_wires = max_wires_set;
        }

        // Offset adjustment for wire crossing points
        const Vector adjust = offset * (
            coords.ray_crossing({smin.layer, 1}, {smax.layer, 1}) -
            coords.ray_crossing({smin.layer, 0}, {smax.layer, 0}));

        const double pitch_adjust = offset * (
            coords.pitch_location({smin.layer, 1}, {smax.layer, 1}, smid.layer) -
            coords.pitch_location({smin.layer, 0}, {smax.layer, 0}, smid.layer));

        std::vector<Point> points;


        //  // Collect wire indices for each plane //debug ...
        // std::vector<int> wires_u, wires_v, wires_w;
        // for (auto idx : actual_min_wires) {
        //     if (smin.layer - ndummy_index == 0) wires_u.push_back(idx);
        //     if (smin.layer - ndummy_index == 1) wires_v.push_back(idx);
        //     if (smin.layer - ndummy_index == 2) wires_w.push_back(idx);
        // }
        // for (auto idx : actual_max_wires) {
        //     if (smax.layer - ndummy_index == 0) wires_u.push_back(idx);
        //     if (smax.layer - ndummy_index == 1) wires_v.push_back(idx);
        //     if (smax.layer - ndummy_index == 2) wires_w.push_back(idx);
        // }
        // // Add actual_mid_wires collection
        // std::set<decltype(smid.bounds.first)> actual_mid_wires;
        // if (use_all_wires) {
        //     for (auto i = smid.bounds.first; i < smid.bounds.second; i++) {
        //         actual_mid_wires.insert(i);
        //     }
        // } else {
        //     // Only use the bounds (start and end-1)
        //     actual_mid_wires.insert(smid.bounds.first);
        //     actual_mid_wires.insert(smid.bounds.second-1);
        // }
        // for (auto idx : actual_mid_wires) {
        //     if (smid.layer - ndummy_index == 0) wires_u.push_back(idx);
        //     if (smid.layer - ndummy_index == 1) {wires_v.push_back(idx); }
        //     if (smid.layer - ndummy_index == 2) {wires_w.push_back(idx); }
        // }
        // // Remove duplicates
        // std::sort(wires_u.begin(), wires_u.end());
        // wires_u.erase(std::unique(wires_u.begin(), wires_u.end()), wires_u.end());
        // std::sort(wires_v.begin(), wires_v.end());
        // wires_v.erase(std::unique(wires_v.begin(), wires_v.end()), wires_v.end());
        // std::sort(wires_w.begin(), wires_w.end());
        // wires_w.erase(std::unique(wires_w.begin(), wires_w.end()), wires_w.end());
        bool flag_print = false;
        // Print debug info if specific sizes are matched
        // if (wires_u.size() == 10 && wires_v.size() == 10 && wires_w.size() == 3) {
        //     std::cout << "Xin1: " << points.size() << " " << wires_u.size() << " " << wires_v.size() << " " << wires_w.size()
        //           << " " << actual_max_wires.size() << " " << actual_min_wires.size() << " " << disable_mix_dead_cell << " " << use_all_wires
        //           << " bad_plane_max=" << plane_is_bad[max_id - ndummy_index] << " " << max_id 
        //           << " bad_plane_min=" << plane_is_bad[min_id - ndummy_index] << " " << min_id
        //           << " bad_plane_other=" << plane_is_bad[mid_id - ndummy_index] << " " << mid_id
        //           << std::endl;
        //     flag_print = true;
        // }
        // plane_is_bad[mid_id - ndummy_index] = is_plane_bad(mid_id, activity, iface, flag_print);
        // debug code 


        for (auto it_gmax = actual_max_wires.begin(); it_gmax != actual_max_wires.end(); it_gmax++) {
            coordinate_t cmax{smax.layer, *it_gmax};
            
            bool flag_must2 = max_wires_set.find(*it_gmax) != max_wires_set.end();
            double charge2 = get_wire_charge(cmax, activity, iface, ndummy_index, flag_print);

            // if (!disable_mix_dead_cell){
            //    std::cout << "max: " << charge2 << " " << *it_gmax << std::endl;
            // }
            
            if ((!flag_must2) && (charge2 < thresh_max) && (charge2 != 0 || disable_mix_dead_cell)) {
                continue;
            }
            // if (flag_print) {
            //     std::cout << "wire: " << *it_gmax << " " << charge2 << " " << flag_must2 << " " << (charge2 < thresh_max) << " " << (charge2 != 0 || disable_mix_dead_cell) << " " << thresh_max << std::endl;
            // }
           

            for (auto it_gmin = actual_min_wires.begin(); it_gmin != actual_min_wires.end(); it_gmin++) {
                coordinate_t cmin{smin.layer, *it_gmin};
                
                // Check if this is a "must" wire (from stepped set)
                bool flag_must1 = min_wires_set.find(*it_gmin) != min_wires_set.end();
                
                // Get charge for this wire
                double charge1 = get_wire_charge(cmin, activity, iface, ndummy_index);
                
                // if(!disable_mix_dead_cell){
                //     std::cout << "min: " << charge1 << " " << *it_gmin << std::endl;
                // }

                // Apply charge filtering for non-must wires
                if ((!flag_must1) && (charge1 < thresh_min) && (charge1 != 0 || disable_mix_dead_cell)) {
                    continue;
                }
                // if (flag_print) {
                //     std::cout << "min wire: " << *it_gmin << " " << charge1 << " " << flag_must1 << " " << (charge1 < thresh_min) << " " << (charge1 != 0 || disable_mix_dead_cell) << " " << thresh_min << std::endl;
                // }
           

                // Check basic bounds
                const double ploc0 = coords.pitch_location(cmin, cmax, 0);
                const double prel0 = coords.pitch_relative(ploc0, 0);
                const double ploc1 = coords.pitch_location(cmin, cmax, 1);
                const double prel1 = coords.pitch_relative(ploc1, 1);
                if (prel0 > 1 || prel0 < 0 || prel1 > 1 || prel1 < 0) {
                    continue;
                }

                // Check third plane bounds
                const double pitch = coords.pitch_location(cmin, cmax, smid.layer)+ pitch_adjust;
                const double pitch_relative = coords.pitch_relative(pitch, smid.layer);
                
                if (pitch_relative > smid.bounds.first - tolerance && 
                    pitch_relative < smid.bounds.second + tolerance) {
                    
                    // Get charge for third plane
                    coordinate_t cother{smid.layer, static_cast<int>(std::round(pitch_relative))};  
                    double charge3 = get_wire_charge(cother, activity, iface, ndummy_index);
                    
                    // Apply charge validation logic
                    if (flag_must1 && flag_must2) {
                        // Both wires are mandatory, no additional charge filtering
                    } else {
                        // At least one wire is not mandatory, apply charge filtering
                        if ((charge2 < thresh_max && (charge2 != 0 || disable_mix_dead_cell)) || // 2 is max ...
                            (charge1 < thresh_min && (charge1 != 0 || disable_mix_dead_cell)) || // 1 is min ...
                            (charge3 < thresh_other && (charge3 != 0 || disable_mix_dead_cell))) {
                       
                        
                            continue;
                        }
                        
                        // Skip if all charges are zero
                        if (charge1 == 0 && charge2 == 0 && charge3 == 0) {
                            continue;
                        }
                    }
                    
                    // if (flag_print) {
                    //     std::cout <<  cmax << " " << cmin 
                    //               << " " << cother 
                    //               << " (" << charge2 << ", " << charge1 << ", " << charge3 << ")" << " " << pitch << " " << pitch_relative << " " 
                    //               << " " << thresh_max << " " << thresh_min << " " << thresh_other << " " << std::endl;
                    // }

                    const auto pt = coords.ray_crossing(cmin, cmax);
                    points.push_back(pt + adjust);
                }
            }
        }

        // if (!disable_mix_dead_cell) std::cout << points.size() << std::endl;

        //de bug ...
        // if (wires_u.size() == 10 && wires_v.size() == 10 && wires_w.size() == 3) {
        //      std::cout << "Xin2: " << points.size() << " " << wires_u.front()  << " " << wires_v.front() << " " << wires_w.front() << " " << wires_u.size() << " " << wires_v.size() << " " << wires_w.size()
        //           << " " << actual_max_wires.size() << " " << actual_min_wires.size() << " " << disable_mix_dead_cell << " " << use_all_wires << std::endl;
        // }
        // debug ...

        intern(ds, points);

        // Add auxiliary data
        aux.add("max_wire_interval", Array({(int)nmax}));
        aux.add("min_wire_interval", Array({(int)nmin}));
        aux.add("max_wire_type", Array({(int)(max_id-ndummy_index)}));
        aux.add("min_wire_type", Array({(int)(min_id-ndummy_index)}));
        // aux.add("charge_threshold_max", Array({thresh_max}));
        // aux.add("charge_threshold_min", Array({thresh_min}));
        // aux.add("charge_threshold_other", Array({thresh_other}));
        // aux.add("use_all_wires", Array({use_all_wires}));
        // aux.add("bad_plane_max", Array({plane_is_bad[max_id - ndummy_index]}));
        // aux.add("bad_plane_min", Array({plane_is_bad[min_id - ndummy_index]}));
        // aux.add("bad_plane_other", Array({plane_is_bad[mid_id - ndummy_index]}));
    }

private:
    // Helper function to detect if a plane is bad based on charge uncertainty
    // Check the first wire of the blob for each plane
    bool is_plane_bad(int plane_layer, 
                      const ISlice::map_t& activity, 
                      IAnodeFace::pointer iface, bool flag_print = false) {
        
        // Get the blob strips to find the first wire index for this plane
        auto strips = iblob->shape().strips();
        const int ndummy_index = strips.size() == 5 ? 2 : 0;
        

        if (plane_layer -ndummy_index < 0 || plane_layer - ndummy_index >= (int)iface->planes().size()) {
            return false;
        }
        
        // Find the strip for this plane layer
        const Strip* target_strip = nullptr;
        for (const auto& strip : strips) {
            if (strip.layer == plane_layer) {
                target_strip = &strip;
                break;
            }
        }
        
        if (!target_strip) {
            return false;
        }

        // if (flag_print) std::cout << plane_layer << " " <<  (int)iface->planes().size() << target_strip << std::endl;

        
        // Get the first and last wire indices from the strip bounds
        int first_wire_index = target_strip->bounds.first;
        int last_wire_index = target_strip->bounds.second - 1;

        // Lambda to check if a wire is bad based on uncertainty
        auto is_wire_bad = [&](int wire_index) -> bool {
            // Create coordinate for the wire
            coordinate_t coord{plane_layer, wire_index};

           

            // Get the appropriate plane
            int plane_index = plane_layer - ndummy_index;

            // if (flag_print) {
            //     std::cout << "Plane check: " << plane_index << " " <<  (int)iface->planes().size() << std::endl;
            // }

            if (plane_index < 0 || plane_index >= (int)iface->planes().size()) {
            return false;
            }

            auto iplane = iface->planes()[plane_index];
            const IWire::vector& iwires = iplane->wires();
            const IChannel::vector& channels = iplane->channels();

            // if (flag_print) {
            //     std::cout << "Wire check: " << wire_index << " " << (int)iwires.size() << std::endl;
            // }

            // Bounds check for wire index
            if (wire_index < 0 || wire_index >= (int)iwires.size()) {
            return false;
            }

            // Get the wire and its channel
            IWire::pointer iwire = iwires[wire_index];
            int channel_ident = iwire->channel();

            // Build/use cache for channel ident to index mapping (same as get_wire_charge)
            auto& p_chi2i = plane_ident2index[iplane];
            if (p_chi2i.empty()) {
            const size_t nchannels = channels.size();
            for (size_t chind=0; chind<nchannels; ++chind) {
                auto ich = channels[chind];
                p_chi2i[ich->ident()] = chind;
            }
            }

            // if (flag_print){
            //     std::cout << "Channel ident: " << channel_ident << " " << p_chi2i.size() << std::endl;
            // }

            // Look up channel index using the cache
            auto chi2i_it = p_chi2i.find(channel_ident);
            if (chi2i_it == p_chi2i.end()) {
            return false;
            }

            // if (flag_print) {
            //     std::cout << "Found channel index: " << chi2i_it->second << " " << (int)channels.size() << std::endl;
            // }

            int channel_attach = chi2i_it->second;
            if (channel_attach < 0 || channel_attach >= (int)channels.size()) {
            return false;
            }

            auto ich = channels[channel_attach];

            // if (flag_print){
            //     std::cout << "Checking channel: " << ich << " " << activity.size() << std::endl;
            // }

            // Look up charge in activity map and check uncertainty
            auto ait = activity.find(ich);
            if (ait != activity.end()) {
                auto act = ait->second;
                double uncertainty = act.uncertainty();

                // if (flag_print) {
                //     std::cout << "Checking wire " << wire_index << " in plane " << plane_layer 
                //             << ": charge=" << act.value() << ", uncertainty=" << uncertainty 
                //             << ", threshold=" << dead_threshold << std::endl;
                // }

                // Plane is considered bad if uncertainty exceeds threshold
                return uncertainty > dead_threshold;
            }

            return false;
        };

        // Check both first and last wire
        return is_wire_bad(first_wire_index) || is_wire_bad(last_wire_index);
    }


    // Helper function to get wire charge using the same pattern as make_dataset
    double get_wire_charge(const coordinate_t& coord, 
                          const ISlice::map_t& activity, 
                          IAnodeFace::pointer iface, int ndummy_index = 2, bool flag_print = false) {
        
        // Get the appropriate plane
        int plane_index = coord.layer-ndummy_index;
        if (plane_index < 0 || plane_index >= (int)iface->planes().size()) {
            return 0.0;
        }
        
        auto iplane = iface->planes()[plane_index];
        const IWire::vector& iwires = iplane->wires();
        const IChannel::vector& channels = iplane->channels();
        
        // Bounds check for wire index
        if (coord.grid < 0 || coord.grid >= (int)iwires.size()) {
            return 0.0;
        }
        
        // Get the wire and its channel
        IWire::pointer iwire = iwires[coord.grid];
        int channel_ident = iwire->channel();
        
        // Build/use cache for channel ident to index mapping (same as make_dataset)
        auto& p_chi2i = plane_ident2index[iplane];
        if (p_chi2i.empty()) {
            const size_t nchannels = channels.size();
            for (size_t chind=0; chind<nchannels; ++chind) {
                auto ich = channels[chind];
                p_chi2i[ich->ident()] = chind;
            }
        }
        
        // Look up channel index using the cache
        auto chi2i_it = p_chi2i.find(channel_ident);
        if (chi2i_it == p_chi2i.end()) {
            return 0.0;
        }
        
        int channel_attach = chi2i_it->second;
        if (channel_attach < 0 || channel_attach >= (int)channels.size()) {
            return 0.0;
        }
        
        auto ich = channels[channel_attach];
        
        // Look up charge in activity map
        auto ait = activity.find(ich);

        if (ait != activity.end()) {
            auto act = ait->second;
            return act.value();
        }
        
        return 0.0;
    }
};

// ===========================================================================
// BlobSampler Extension for Runtime Configuration
// ===========================================================================

std::tuple<PointCloud::Dataset, PointCloud::Dataset> 
BlobSampler::sample_blob_with_config(const IBlob::pointer& iblob, 
                       int blob_index, 
                       const Configuration& runtime_config) {
    
    PointCloud::Dataset ret_main;
    PointCloud::Dataset ret_aux;
    
    for (auto& sampler : m_samplers) {
        sampler->begin_sample(blob_index, iblob);
        
        // Check if this is a ChargeStepped sampler
        auto charge_stepped = dynamic_cast<ChargeStepped*>(sampler.get());
        if (charge_stepped) {
            // Use the configuration-aware sampling
            charge_stepped->sample_with_config(ret_main, ret_aux, runtime_config);
        } else {
            // Use regular sampling
            sampler->sample(ret_main, ret_aux);
        }
        
        sampler->end_sample();
    }
    
    return {ret_main, ret_aux};
}



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
    if (startswith(name, "charge_stepped")) {
        m_samplers.push_back(std::make_unique<ChargeStepped>(full, m_samplers.size()));
        m_samplers.back()->configure(full);
        return;
    }

    THROW(ValueError() << errmsg{"unknown strategy: " + name});
}


