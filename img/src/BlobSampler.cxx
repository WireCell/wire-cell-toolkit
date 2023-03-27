#include "WireCellImg/BlobSampler.h"

#include "WireCellUtil/Range.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

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


// Base class providing common functionality to all samplers.
struct BlobSampler::Sampler
{

    // Hard-wired, common subset of full config.
    BlobSampler::CommonConfig cc;

    // The identity of the sampler
    size_t ident;

    explicit Sampler(const Configuration& cfg, size_t ident)
        : cc(cfg2cpp(cfg)), ident(ident)
    {
        this->configure(cfg);
    }

    virtual ~Sampler() {}

    // Context manager around sample()
    IAnodeFace::pointer anodeface;
    void begin_sample(IBlob::pointer iblob)
    {
        anodeface = iblob->face();
    }
    void end_sample()
    {
        anodeface = nullptr;
    }

    // Entry point to subclasses
    virtual void sample(Dataset& ds, IBlob::pointer iblob) = 0;

    // subclass may want to config self.
    virtual void configure(const Configuration& cfg) { }

    // Return pimpos for a given plane index
    const Pimpos* pimpos(int plane_index=2) const
    {
        return anodeface->planes()[plane_index]->pimpos();
    }

    // Convert a signal time to its location along global drift
    // coordinates.
    double time2drift(double time) const
    {
        const Pimpos* colpimpos = pimpos(2);
        const double drift = (time + cc.time_offset)/cc.drift_speed;
        double xorig = colpimpos->origin()[0];
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
        if (corners.empty()) return avg;
        for (const auto& corner : corners) {
            avg += crossing_point(corner);
        }
        avg = avg / corners.size();
        return avg;
    }

    // Return a dataset with a single point data.
    // Dataset make_dataset(const Point& p, int blobid)
    // {
    //     return Dataset({
    //             {cc.prefix + "x", Array({p.x()})},
    //             {cc.prefix + "y", Array({p.y()})},
    //             {cc.prefix + "z", Array({p.z()})},
    //             {cc.prefix + "strategy", Array({ident})},
    //             {cc.prefix + "blobid", Array({blobid})}});
    // }
    // Return a dataset covering multiple points related to a blob
    Dataset make_dataset(const std::vector<Point>& pts, int blobid, double time)
    {
        size_t npts = pts.size();
        std::vector<size_t> idents(npts, ident);
        std::vector<int> blobids(npts, blobid);
        std::vector<float> times(npts, time);
        std::vector<Point::coordinate_t> x(npts),y(npts),z(npts);
        for (size_t ind=0; ind<npts; ++ind) {
            const auto& pt = pts[ind];
            x[ind] = pt.x();
            y[ind] = pt.y();
            z[ind] = pt.z();
        }

        return Dataset({
                {cc.prefix + "x", Array(x)},
                {cc.prefix + "y", Array(y)},
                {cc.prefix + "z", Array(z)},
                {cc.prefix + "t", Array(times)},
                {cc.prefix + "strategy", Array(idents)},
                {cc.prefix + "blobid", Array(blobids)}});

        // fixme: include:
        // - pitch in each view
        // - index of nearest wire in each view
        // - tick
        // - activity

    }

    // Append points to PC.  
    void intern(Dataset& ds,
                std::vector<Point> points,
                IBlob::pointer iblob)
    {
        const int blobid = iblob->ident();
        auto islice = iblob->slice();
        const double t0 = islice->start();
        const double dt = islice->span();

        const Binning bins(cc.tbinning.nbins(),
                           cc.tbinning.min() * t0,
                           cc.tbinning.max() * (t0 + dt));

        const size_t npts = points.size();
        for (int tbin : irange(bins.nbins())) {
            const double time = bins.edge(tbin);
            const double x = time2drift(time);
            for (size_t ind=0; ind<npts; ++ind) {
                points[ind].x(x);
            }
            auto tail = make_dataset(points, blobid, time);
            ds.append(tail);
        }
    }
};


PointCloud::Dataset BlobSampler::sample_blobs(const IBlob::vector& iblobs)
{
    PointCloud::Dataset ret;
    for (const auto& iblob : iblobs) {
        for (auto& sampler : m_samplers) {
            sampler->begin_sample(iblob);
            sampler->sample(ret, iblob);
            sampler->end_sample();
        }
    }
    return ret;
}
        

struct Center : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Center(const Center&) = default;
    Center& operator=(const Center&) = default;

    void sample(Dataset& ds, IBlob::pointer iblob)
    {
        std::vector<Point> points(1);
        points[0] = center_point(iblob->shape().corners());
        intern(ds, points, iblob);
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
    void sample(Dataset& ds, IBlob::pointer iblob)
    {
        const auto& corners = iblob->shape().corners();
        const size_t npts = corners.size();
        if (! npts) return;
        std::vector<Point> points;
        points.reserve(npts);
        for (const auto& corner : corners) {
            points.emplace_back(crossing_point(corner));
        }
        intern(ds, points, iblob);
    }
};


struct Edge : public BlobSampler::Sampler
{
    using BlobSampler::Sampler::Sampler;
    Edge(const Edge&) = default;
    Edge& operator=(const Edge&) = default;

    void sample(Dataset& ds, IBlob::pointer iblob)
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
        intern(ds, points, iblob);
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

    void sample(Dataset& ds, IBlob::pointer iblob)
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
        intern(ds, points, iblob);
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
        intern(ds, points, iblob);
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

    void sample(Dataset& ds, IBlob::pointer iblob)
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
        intern(ds, points, iblob);
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

    std::string name = full["name"].asString();
    // use startswith() to be a little friendly to the user
    // w.r.t. spelling.  eg "centers", "bounds", "boundaries" are all
    // accepted.
    if (startswith(name, "center")) {
        m_samplers.push_back(std::make_unique<Center>(full, m_samplers.size()));
        return;
    }
    if (startswith(name, "corner")) {
        m_samplers.push_back(std::make_unique<Corner>(full, m_samplers.size()));
        return;
    }
    if (startswith(name, "edge")) {
        m_samplers.push_back(std::make_unique<Edge>(full, m_samplers.size()));
        return;
    }
    if (startswith(name, "grid")) {
        m_samplers.push_back(std::make_unique<Grid>(full, m_samplers.size()));
        return;
    }
    if (startswith(name, "bound")) {
        m_samplers.push_back(std::make_unique<Bounds>(full, m_samplers.size()));
        return;
    }

    THROW(ValueError() << errmsg{"unknown strategy: " + name});
}

        
