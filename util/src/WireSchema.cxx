#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Intersection.h"

#include <numeric>

using namespace WireCell;
using namespace WireCell::WireSchema;

static // for now
Ray ray_pitch_approx(const Ray& r1, const Ray& r2)
{
    // return ray_pitch(r1, r2);

    const auto c1 = 0.5*(r1.first + r1.second);
    const auto c2 = 0.5*(r2.first + r2.second);
    const auto d21 = c2 - c1;

    const auto v1 = ray_vector(r1);
    const auto ecks = v1.cross(d21).norm();
    const auto pdir = ecks.cross(v1).norm();
    const double pitch = d21.dot(pdir);
    return Ray(c1, c1+pitch*pdir);
}


static
void load_file(const std::string& path, StoreDB& store)
{
    Json::Value jtop = WireCell::Persist::load(path);
    Json::Value jstore = jtop["Store"];

    std::vector<Point> points;
    {
        Json::Value jpoints = jstore["points"];
        const int npoints = jpoints.size();
        points.resize(npoints);
        for (int ipoint = 0; ipoint < npoints; ++ipoint) {
            Json::Value jp = jpoints[ipoint]["Point"];
            points[ipoint].set(get<double>(jp, "x"), get<double>(jp, "y"), get<double>(jp, "z"));
        }
    }

    {  // wires
        Json::Value jwires = jstore["wires"];
        const int nwires = jwires.size();
        std::vector<Wire>& wires = store.wires;
        wires.resize(nwires);
        for (int iwire = 0; iwire < nwires; ++iwire) {
            Json::Value jwire = jwires[iwire]["Wire"];
            Wire& wire = wires[iwire];
            wire.ident = get<int>(jwire, "ident");
            wire.channel = get<int>(jwire, "channel");
            wire.segment = get<int>(jwire, "segment");

            int itail = get<int>(jwire, "tail");
            int ihead = get<int>(jwire, "head");
            wire.tail = points[itail];
            wire.head = points[ihead];
        }
    }

    {  // planes
        Json::Value jplanes = jstore["planes"];
        const int nplanes = jplanes.size();
        std::vector<Plane>& planes = store.planes;
        planes.resize(nplanes);
        for (int iplane = 0; iplane < nplanes; ++iplane) {
            Json::Value jplane = jplanes[iplane]["Plane"];
            Plane& plane = planes[iplane];
            plane.ident = get<int>(jplane, "ident");
            Json::Value jwires = jplane["wires"];
            const int nwires = jwires.size();
            plane.wires.resize(nwires);
            for (int iwire = 0; iwire < nwires; ++iwire) {
                plane.wires[iwire] = convert<int>(jwires[iwire]);
            }
        }
    }

    {  // faces
        Json::Value jfaces = jstore["faces"];
        const int nfaces = jfaces.size();
        std::vector<Face>& faces = store.faces;
        faces.resize(nfaces);
        for (int iface = 0; iface < nfaces; ++iface) {
            Json::Value jface = jfaces[iface]["Face"];
            Face& face = faces[iface];
            face.ident = get<int>(jface, "ident");
            Json::Value jplanes = jface["planes"];
            const int nplanes = jplanes.size();
            face.planes.resize(nplanes);
            for (int iplane = 0; iplane < nplanes; ++iplane) {
                face.planes[iplane] = convert<int>(jplanes[iplane]);
            }
        }
    }

    {  // anodes
        Json::Value janodes = jstore["anodes"];
        const int nanodes = janodes.size();
        std::vector<Anode>& anodes = store.anodes;
        anodes.resize(nanodes);
        for (int ianode = 0; ianode < nanodes; ++ianode) {
            Json::Value janode = janodes[ianode]["Anode"];
            Anode& anode = anodes[ianode];
            anode.ident = get<int>(janode, "ident");
            Json::Value jfaces = janode["faces"];
            const int nfaces = jfaces.size();
            anode.faces.resize(nfaces);
            for (int iface = 0; iface < nfaces; ++iface) {
                anode.faces[iface] = convert<int>(jfaces[iface]);
            }
        }
    }

    {  // detectors
        std::vector<Detector>& dets = store.detectors;

        Json::Value jdets = jstore["detectors"];
        const int ndets = jdets.size();

        if (!ndets) {           // some do not give us this section
            dets.resize(1);
            auto& det = dets[0];
            det.ident = 0;
            int nanodes = store.anodes.size();
            det.anodes.resize(nanodes);
            for (int ianode = 0; ianode < nanodes; ++ianode) {
                det.anodes[ianode] = ianode;
            }
        }
        else {
            dets.resize(ndets);
            for (int idet = 0; idet < ndets; ++idet) {
                Json::Value jdet = jdets[idet]["Detector"];
                Detector& det = dets[idet];
                det.ident = get<int>(jdet, "ident");
                Json::Value janodes = jdet["anodes"];
                const int nanodes = janodes.size();
                det.anodes.resize(nanodes);
                for (int ianode = 0; ianode < nanodes; ++ianode) {
                    det.anodes[ianode] = convert<int>(janodes[ianode]);
                }
            }
        }
    }
}

// Return axis along which wire centers are ascending when in proper
// wire-in-plane order.
static int wire_order_axis(const StoreDB& store, const Plane& plane)
{
    const auto& w = store.wires[plane.wires[0]];
    Vector wdir = ray_unit(Ray(w.tail, w.head));
    const double near_unity = 0.9999;
    if (std::abs(wdir.z()) > near_unity) { 
        return 1;             // less than ~1 deg from Zthen sort by Y
    }
    return 2;                   // sort by Z
}

// Sort wires according to their center coordinates
struct wip_order {
    const StoreDB &store;
    size_t axis{2};

    bool operator()(int a, int b) const {
        const auto& wa = store.wires[a];
        const auto& wb = store.wires[b];        
        return wa.head[axis] + wa.tail[axis] < wb.head[axis] + wb.tail[axis];
    }
};

using plane_fixer_f = std::function<void(StoreDB& store, Plane& plane)>;
using plane_fixers_t = std::vector<plane_fixer_f>;

// Fix order of plane's wires-in-plane and wire endpoints.
static void plane_fixer_order(StoreDB& store, Plane& plane)
{    

    const size_t nwires = plane.wires.size();

    std::vector<double> ycen(nwires, 0.0), zcen(nwires, 0.0);

    // Capture wire centers in original wire order
    for (size_t ind=0; ind<nwires; ++ind) {
        Wire& wire = store.wires[plane.wires[ind]];
        ycen[ind] = 0.5*(wire.tail.y() + wire.head.y());
        zcen[ind] = 0.5*(wire.tail.z() + wire.head.z());
    }

    // Find an origin as the center of centers and an approximate
    // pitch direction going from minimum to maximum center.
    double ymin, ymax, zmin, zmax;
    const int axis = wire_order_axis(store, plane);
    if (1 == axis) {            // minmax by Y
        const auto [emin, emax] = std::minmax_element(ycen.begin(), ycen.end());        
        ymin = *emin;
        ymax = *emax;
        zmin = zcen[std::distance(ycen.begin(), emin)];
        zmax = zcen[std::distance(ycen.begin(), emax)];
    }
    else {                      // minmax by Z
        const auto [emin, emax] = std::minmax_element(zcen.begin(), zcen.end());
        ymin = ycen[std::distance(zcen.begin(), emin)];
        ymax = ycen[std::distance(zcen.begin(), emax)];
        zmin = *emin;
        zmax = *emax;
    }
    Point origin(0, 0.5*(ymin+ymax), 0.5*(zmin+zmax));
    Vector dir(0, ymax-ymin, zmax-zmin);
    dir = dir.norm();

    // Find a position of wire centers from origin along direction.
    std::vector<double> pos(nwires, 0.0);
    for (size_t ind=0; ind<nwires; ++ind) {
        Point pt(0, ycen[ind], zcen[ind]);
        pos[ind] = dir.dot(pt - origin);
    }
    
    // This holds original indices into the plane.wires vector.
    std::vector<size_t> indices(nwires);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Sort those indices according to position.
    std::sort(indices.begin(), indices.end(),
              [&] (size_t a, size_t b) -> bool {
                  return pos[a] < pos[b];
              });

    // Replace the ordering to plane.wires.
    std::vector<int> new_wires(nwires);
    for (size_t ind=0; ind<nwires; ++ind) {
        new_wires[ind] = plane.wires[indices[ind]];
    }
    plane.wires = new_wires;

    // Maybe flip endpoints to get correct wire direction.
    for (int iwire : plane.wires) {
        Wire& wire = store.wires[iwire];

        if (axis == 1) {
            if (wire.head.z() > wire.tail.z()) {
                std::swap(wire.head, wire.tail);
            }
        }
        else {
            if (wire.head.y() < wire.tail.y()) {
                std::swap(wire.head, wire.tail);
            }
        }
    }
}

// Fix wire directions.  Assumes wire endpoint order is correct.
static void plane_fixer_direction(StoreDB& store, Plane& plane)
{
    Vector wdir;
    const size_t nwires = plane.wires.size();
    std::vector<double> half(nwires); // retain wire lengths

    for (size_t wind=0; wind<nwires; ++wind) {
        const Wire wire = store.wires[plane.wires[wind]];
        const Ray ray(wire.tail, wire.head);
        const auto rv = ray_vector(ray);
        wdir += rv;
        half[wind] = 0.5*rv.magnitude();
    }

    wdir.x(0);                  // bring into Y-Z plane
    wdir = wdir.norm();

    // Correct wire endpoints.
    for (size_t wind=0; wind<nwires; ++wind) {
        auto& wire = store.wires[plane.wires[wind]];
        const auto c = 0.5*(wire.tail + wire.head);
        const Vector whalf = wdir * half[wind];
        wire.head = c + whalf;
        wire.tail = c - whalf;
    }        
}

// uniform pitch, coplanar wires.  Assumes directions are parallel.
static void plane_fixer_pitch(StoreDB& store, Plane& plane)
{
    // Find the mean X position of wire centers
    double xmean = 0;
    const size_t nwires = plane.wires.size();
    for (size_t wind=0; wind<nwires; ++wind) {
        const auto& wire = store.wires[plane.wires[wind]];
        xmean += 0.5*(wire.tail.x() + wire.head.x());
    }
    xmean /= nwires;

    // mean pitch.  average all pitches.  Not sure best strategy.
    // Imprecise endoints lead to less precision on short wires but
    // including short wires means averaging over more.  Do simplest
    // thing and average over all of them.

    const size_t nhalf = nwires/2;
    Ray midway;                 // pick up central wire 

    Vector ptot;
    Ray prev;

    for (size_t wind=0; wind<nwires; ++wind) {
        auto& wire = store.wires[plane.wires[wind]];
        wire.tail.x(xmean);     // move wire to
        wire.head.x(xmean);     // average X location
        Ray next(wire.tail, wire.head);
        if (wind == nhalf) {
            midway = next;
        }
        if (wind) {            // wait until 2nd to calculate diff
            ptot += ray_vector(ray_pitch(prev, next));
            // ptot += ray_vector(ray_pitch_approx(prev, next));
        }
        prev = next;
    }
    const Vector pmean = ptot / (nwires-1);
    const double pmag = pmean.magnitude();
    const Vector pdir = pmean/pmag; // ptot.norm();

    // Center point of midway wire
    const Vector origin = 0.5*(midway.first + midway.second);

    for (size_t wind=0; wind<nwires; ++wind) {
        auto& wire = store.wires[plane.wires[wind]];

        // Center of wire relative to origin
        const auto wcen = 0.5*(wire.tail + wire.head) - origin;

        if (wind == nhalf) {
            // don't correct fixed midway wire
            continue;
        }

        const double have_pitch = pdir.dot(wcen);
        const double want_pitch = ((int)wind - (int)nhalf)*pmag;
        const double delta_pitch = want_pitch - have_pitch;
        const Vector diff = delta_pitch*pdir;

        wire.tail += diff;
        wire.head += diff;
    }

}

static void fix_planes(StoreDB& store, plane_fixers_t& fixers)
{
    if (fixers.empty()) return;

    for (auto& detector : store.detectors) {

        for (int ianode : detector.anodes) {
            auto& anode = store.anodes[ianode];

            for (int iface : anode.faces) {
                auto& face = store.faces[iface];

                for (int iplane : face.planes) {
                    auto& plane = store.planes[iplane];

                    for (auto& fixer : fixers) {
                        fixer(store, plane);
                    }
                }
            }
        }
    }
}



Store WireCell::WireSchema::load(const char* filename, Correction correction)
{
    using cache_key_t = std::pair<std::string, WireSchema::Correction>;
    static std::map<cache_key_t, StoreDBPtr> cache;

    // turn into absolute real path
    std::string realpath = WireCell::Persist::resolve(filename);

    // Make mutable shared pointer to shart with.
    StoreDBPtr cached;
    Correction have_level = Correction::empty;

    // See if we have already loaded this file at current or lower
    // correction level
    for (auto level = correction; level > Correction::empty; --level) {
        const cache_key_t key(realpath, level);
        auto maybe = cache.find(key);

        if (maybe == cache.end()) {
            continue;
        }

        have_level = level;
        cached = maybe->second;
        break;
    }
    if (have_level == correction) {
        return Store(cached);   // nothing to do
    }
    
    auto store = std::make_shared<StoreDB>();
    if (cached) {
        store = std::make_shared<StoreDB>(*cached);
    }

    const cache_key_t ckey(realpath, correction);
    cache[ckey] = store;

    // always load if we start empty
    if (have_level == Correction::empty) {
        load_file(realpath, *store);
        have_level = Correction::load;
    };

    // Pile on fixers to reach correction level
    plane_fixers_t fixers;
    while (have_level < correction) {
        if (have_level == Correction::load)
            fixers.push_back(plane_fixer_order);
        if (have_level == Correction::order)
            fixers.push_back(plane_fixer_direction);
        if (have_level == Correction::direction)
            fixers.push_back(plane_fixer_pitch);
        ++have_level;
    };

    fix_planes(*store, fixers);

    return Store(store);
}

static Json::Value to_json(const std::vector<int>& arr)
{
    Json::Value ret = Json::arrayValue;
    size_t siz = arr.size();
    ret.resize(siz);
    for (size_t ind=0; ind<siz; ++ind) {
        ret[(int)ind] = arr[ind];
    }
    return ret;    
}

static Json::Value jo_wrap(const char* name, Json::Value& j)
{
    Json::Value ret = Json::objectValue;
    ret[name] = j;
    return ret;
}

void WireCell::WireSchema::dump(const char* filename, const Store& store)
{
    Json::Value jstore = Json::objectValue;

    Json::Value jpoints = Json::arrayValue; // we don't bother deduplicating
    Json::Value jwires = Json::arrayValue;
    for (const auto& wire : store.wires()) {
        Json::Value jwire = Json::objectValue;
        jwire["ident"] = wire.ident;
        jwire["channel"] = wire.channel;
        jwire["segment"] = wire.segment;

        jwire["tail"] = jpoints.size();
        Json::Value tpoint = Json::objectValue;
        tpoint["x"] = wire.tail.x();
        tpoint["y"] = wire.tail.y();
        tpoint["z"] = wire.tail.z();
        jpoints.append(jo_wrap("Point", tpoint));

        jwire["head"] = jpoints.size();
        Json::Value hpoint = Json::objectValue;
        hpoint["x"] = wire.head.x();
        hpoint["y"] = wire.head.y();
        hpoint["z"] = wire.head.z();
        jpoints.append(jo_wrap("Point", hpoint));

        jwires.append(jo_wrap("Wire", jwire));
    }
    jstore["points"] = jpoints;
    jstore["wires"] = jwires;

    Json::Value jplanes = Json::arrayValue;
    for (const auto& plane : store.planes()) {
        Json::Value jplane = Json::objectValue;
        jplane["ident"] = plane.ident;
        jplane["wires"] = to_json(plane.wires);
        jplanes.append(jo_wrap("Plane", jplane));
    }
    jstore["planes"] = jplanes;

    Json::Value jfaces = Json::arrayValue;
    for (const auto& face : store.faces()) {
        Json::Value jface = Json::objectValue;
        jface["ident"] = face.ident;
        jface["planes"] = to_json(face.planes);
        jfaces.append(jo_wrap("Face",jface));
    }
    jstore["faces"] = jfaces;

    Json::Value janodes = Json::arrayValue;
    for (const auto& anode : store.anodes()) {
        Json::Value janode = Json::objectValue;
        janode["ident"] = anode.ident;
        janode["faces"] = to_json(anode.faces);
        janodes.append(jo_wrap("Anode", janode));
    }
    jstore["anodes"] = janodes;

    const auto& detectors = store.detectors();
    if (! detectors.empty()) {
        Json::Value jdetectors = Json::arrayValue;
        for (const auto& det: detectors) {
            Json::Value jdet = Json::objectValue;
            jdet["ident"] = det.ident;
            jdet["anodes"] = to_json(det.anodes);
            jdetectors.append(jo_wrap("Detector", jdet));
        }
        jstore["detectors"] = jdetectors;
    }

    Json::Value jtop = Json::objectValue;
    jtop["Store"] = jstore;
    Persist::dump(filename, jtop, true);
}


Store::Store()
  : m_db(nullptr)
{
}

Store::Store(StoreDBPtr db)
  : m_db(db)
{
}

Store::Store(const Store& other)
  : m_db(other.db())
{
}
Store& Store::operator=(const Store& other)
{
    m_db = other.db();
    return *this;
}

StoreDBPtr Store::db() const { return m_db; }

const std::vector<Detector>& Store::detectors() const { return m_db->detectors; }
const std::vector<Anode>& Store::anodes() const { return m_db->anodes; }
const std::vector<Face>& Store::faces() const { return m_db->faces; }
const std::vector<Plane>& Store::planes() const { return m_db->planes; }
const std::vector<Wire>& Store::wires() const { return m_db->wires; }

const Anode& Store::anode(int ident) const
{
    for (auto& a : m_db->anodes) {
        if (a.ident == ident) {
            return a;
        }
    }
    THROW(KeyError() << errmsg{String::format("Unknown anode: %d", ident)});
}

std::vector<Anode> Store::anodes(const Detector& detector) const
{
    std::vector<Anode> ret;
    for (auto ind : detector.anodes) {
        ret.push_back(m_db->anodes[ind]);
    }
    return ret;
}

std::vector<Face> Store::faces(const Anode& anode) const
{
    std::vector<Face> ret;
    for (auto ind : anode.faces) {
        ret.push_back(m_db->faces[ind]);
    }
    return ret;
}

std::vector<Plane> Store::planes(const Face& face) const
{
    std::vector<Plane> ret;
    for (auto ind : face.planes) {
        ret.push_back(m_db->planes[ind]);
    }
    return ret;
}

std::vector<Wire> Store::wires(const Plane& plane) const
{
    std::vector<Wire> ret;
    for (auto ind : plane.wires) {
        const auto& w = m_db->wires[ind];
        ret.push_back(w);
    }
    return ret;
}

BoundingBox Store::bounding_box(const Anode& anode, bool region) const
{
    BoundingBox bb;
    for (const auto& face : faces(anode)) {
        bb(bounding_box(face, region).bounds());
    }
    return bb;
}
BoundingBox Store::bounding_box(const Face& face, bool region) const
{
    BoundingBox bb;
    for (const auto& plane : planes(face)) {
        bb(bounding_box(plane, region).bounds());
    }
    return bb;
}
BoundingBox Store::bounding_box(const Plane& plane, bool region) const
{
    BoundingBox bb;
    
    Vector phalf;
    if (region) phalf = 0.5*mean_pitch(plane);

    const auto ws = wires(plane);
    for (const auto& wire : ws) {
        bb(Ray(wire.tail+phalf, wire.head+phalf));
    }
    if (region) {
        bb(Ray(ws[0].tail-phalf, ws[0].head-phalf));
    }
    return bb;
}

// Ray Store::wire_pitch_fast(const Plane& plane) const
// {
//     // Use likely longer wires from center of plane to attenuate any
//     // endpoint imprecision.
//     const size_t nhalf = plane.wires.size()/2;

//     const Wire& w1 = m_db->wires[plane.wires[nhalf]];
//     const Wire& w2 = m_db->wires[plane.wires[nhalf+1]];

//     Ray r1(w1.tail, w1.head);
//     Ray r2(w2.tail, w2.head);

//     const Vector W = ray_unit(r1);
//     const Vector P = ray_unit(ray_pitch(r1, r2));

//     return Ray(W, P);
// }

Vector Store::mean_pitch(const Plane& plane) const
{
    Vector ptot;
    const auto ws = wires(plane);
    const size_t nwires = ws.size();
    Ray prev;
    for (size_t wind=0; wind<nwires; ++wind) {
        const Wire& wire = ws[wind];
        Ray next(wire.tail, wire.head);
        if (wind) {             // wait until 2nd to calculate diff
            //ptot = ptot + ray_vector(ray_pitch(prev, next));
            ptot = ptot + ray_vector(ray_pitch_approx(prev, next));
        }
        prev = next;
    }
    return ptot / (nwires-1);
}

Vector Store::mean_wire(const Plane& plane) const
{
    Vector wtot;
    const auto ws = wires(plane);
    for (const auto& wire : ws) {
        Ray ray(wire.tail, wire.head);
        wtot += ray_vector(ray);
    }
    return wtot / ws.size();
}

Ray Store::wire_pitch(const Plane& plane) const
{
    auto wmean = mean_wire(plane);
    auto pmean = mean_pitch(plane);
    return Ray(wmean.norm(), pmean.norm());
}

std::vector<int> Store::channels(const Plane& plane) const
{
    std::vector<int> ret;
    for (const auto& wire : wires(plane)) {
        ret.push_back(wire.channel);
    }
    return ret;
}


//// Validation ////

using spdlog::error;
using fmt::format;
struct ValidationContext {

    bool fail_fast{false};
    size_t nerrors{0};
    std::string context{"top level"};

    void operator()(const std::string& ctx) {
        context = ctx;
    }

    void bail(const std::string& tested, const std::string& reason)
    {
        ++nerrors;
        const std::string full = format("{}: {} failed: {}", context, tested, reason);
        error(full);
        if (fail_fast) {
            THROW(ValueError() << errmsg{full});
        }
    }
    template<typename T>
    void positive(T val, const std::string& tested)
    {
        if (val > 0) return;
        bail(tested, format("non positive ({})", val));
    }
    void nonneg(int ind, const std::string& tested)
    {
        if (ind >= 0) return;
        bail(tested, format("negative ({})", ind));
    }
    void ubound(int ind, int end, const std::string& tested)
    {
        if (ind < end) return;
        bail(tested, format("{} not less than {}", ind, end));
    }
    void near(double val, double targ, double err, const std::string& tested)
    {
        const double diff = val-targ;
        if (std::abs(diff) < err) return;
        double per = (diff)/err;
        bail(tested, format("{} is {} from {} by more than {} ({} * allowed error)", val, diff, targ, err, per));
    }
    template<typename T>
    void equal(const T& a, const T& b, const std::string& tested)
    {
        if (a == b) return;
        bail(tested, format("{} and {} are not equal", a, b));
    }

    template<typename T>
    const T& element(const std::vector<T>& arr, int ind, const std::string& tested)
    {
        nonneg(ind, tested);
        ubound(ind, arr.size(), tested);
        return arr[ind];
    }
};

void WireCell::WireSchema::validate(const Store& store, double repsilon, bool fail_fast)
{
    ValidationContext v{fail_fast};
    v.positive(store.detectors().size(), "detector count");

    const double near_unity = 0.9999;

    for (const auto& detector : store.detectors()) {
        v.nonneg(detector.ident, "detector ident");
        const std::string dctx = format("detID={}: ", detector.ident);
        v(dctx);

        v.positive(detector.anodes.size(), "anode count");
        for (int ianode : detector.anodes) {
            const auto& anode = v.element(store.anodes(), ianode, "anodes array access");
            v.nonneg(anode.ident, "anode ident");
            const std::string actx = dctx + format("anodeID={}: ", anode.ident);
            v(actx);

            v.positive(anode.faces.size(), "face count");
            for (int iface : anode.faces) {
                const auto& face = v.element(store.faces(), iface, "faces array access");
                v.nonneg(face.ident, "face ident");
                const std::string fctx = actx + format("faceID={}: ", face.ident);
                v(fctx);

                v.positive(face.planes.size(), "plane count");
                std::set<int> face_plane_idents;
                for (int iplane : face.planes) {
                    const auto& plane = v.element(store.planes(), iplane, "planes array access");
                    v.nonneg(plane.ident, "plane ident");
                    const std::string pctx = fctx + format("planeID={}: ", plane.ident);
                    v(pctx);

                    face_plane_idents.insert(plane.ident);

                    const auto pmean = store.mean_pitch(plane);
                    const auto pmmag = pmean.magnitude();
                    const auto pmdir = pmean.norm();
                    const auto wmean = store.mean_wire(plane);
                    const auto wmdir = wmean.norm();
                    v.near(wmdir.dot(pmdir), 0, repsilon, "wire-pitch orthogonality");

                    v.positive(plane.wires.size(), "wire count");
                    Wire wlast;
                    bool seen=false;
                    for (int iwire : plane.wires) {
                        const auto& wire = v.element(store.wires(), iwire, "wires array access");
                        v.nonneg(wire.ident, "wire ident");
                        const std::string wctx = pctx + format("wireID={}: ", wire.ident);
                        v(wctx);

                        // Check wire direction convention
                        Ray wray(wire.tail, wire.head);
                        const auto wvec = ray_vector(wray);
                        const auto wdir = wvec.norm();
                        // const double wlen = wvec.magnitude();

                        if (std::abs(wdir.z()) > near_unity) {
                            v.near(wdir.y(), -1, repsilon, "wire direction Y component");
                        }
                        else {
                            v.positive(wdir.y(), "wire direction Y component");
                        }

                        v.near(wdir.dot(wmdir), 1, repsilon, "local-mean wire parallelism");

                        if (seen) { // validate wire neighbors
                            const auto lray = Ray(wlast.tail, wlast.head);
                            //const auto pray = ray_pitch(lray, wray);
                            const auto pray = ray_pitch_approx(lray, wray);
                            const auto pvec = ray_vector(pray);
                            const auto pdir = pvec.norm();
                            const double pmag = pvec.magnitude();
                            auto ldir = ray_unit(lray);

                            v.near(pmag, pmmag, repsilon*pmmag, "local pitch magnitude");
                            v.near(ldir.dot(pdir), 0, repsilon, "local wire-pitch orthogonality");
                            
                            if (std::abs(wdir.z()) > near_unity) {
                                v.near(pdir.z(), 1, repsilon, "wire pitch Z component");
                            }
                            else {
                                v.positive(pdir.z(), "wire pitch Z component");
                            }
                        }
                        wlast = wire;
                        seen = true;
                    } // wires
                } // planes
                v(fctx);
                v.equal(face_plane_idents.size(), face.planes.size(), "unique plane ident counts");

            } // faces
            
        } // anodes

    } // detectors

    // can't throw in v's destructor
    if (v.nerrors) {
        THROW(ValueError() << errmsg{format("{} wire validation errors", v.nerrors)});
    }

}


WireSchema::Plane& WireSchema::get_append(StoreDB& store, int pid, int fid, int aid, int did)
{
    auto& dobj = get_append<Detector>(store.detectors, did);
    auto& aobj = get_append<Anode>(store.anodes, dobj.anodes, aid);
    auto& fobj = get_append<Face>(store.faces, aobj.faces, fid);
    auto& pobj = get_append<Plane>(store.planes, fobj.planes, pid);
    return pobj;
}


int WireSchema::generate(StoreDB& store, Plane& plane,
                         const Ray& pitch, const Ray& bounds, int wid0)
{
    if (std::abs(pitch.first.x() - pitch.second.x()) > 1e-6) {
        THROW(ValueError() << errmsg{"pitch not orthogonal to drift"});
    }

    const Vector pvec = ray_vector(pitch);
    const Vector pdir = pvec.norm();

    // Obey WCT convention: W = P x X
    const Vector wdir = pdir.cross(Vector(1,0,0));

    Vector cen = pitch.first;
    int wip = 0;
    while (true) {

        Ray wray(cen, cen+wdir), hits;
        int hitmask = box_intersection(bounds, cen, wdir, hits);
        if (hitmask != 3) {
            break;
        }

        plane.wires.push_back(store.wires.size());
        store.wires.emplace_back(Wire{wip+wid0,0,0,hits.first, hits.second});

        ++wip;

        // Choose point on next wire.
        // cen = 0.5 * (hits.first + hits.second) + pvec;
        // cen = cen + pvec;
        cen = pitch.first + wip * pvec;
    }
        
    return wip;
}

/*
  When considering the "ray pairs" code, think of looking in
  the direction of postive X-axis with

  Y
  ^
  |
  |
 (X)----> Z
 
  WCT required conventions:

  Wire direction W: unit vector pointing from wire tail to wire head
  end points and is uniform across all wires in a plane.

  Pitch direction P: unit vector perpendicular to W, coplanar with and
  uniform across all wires in a plane.

  Pitch magnitude "pitch": perpendicular separation between any two
  neighboring wires in a plane.

  Cross product rule: X x W = P

  Pitch vector sign: if P||Z then dot(P,Y) > 0 else dot(P,Z) > 0. 

*/

ray_pair_vector_t WireSchema::ray_pairs_active(const Store& store,
                                               const Face& face,
                                               bool region) 
{
    ray_pair_vector_t raypairs;

    const auto bb = store.bounding_box(face, region).bounds();
    double bbyl=bb.first.y(),  bbzl=bb.first.z();
    double bbyh=bb.second.y(), bbzh=bb.second.z();

    // Corners of a box in X=0 plane, pretend Y points up, Z to right
    Point ll(0.0, bbyl, bbzl);
    Point lr(0.0, bbyl, bbzh);
    Point ul(0.0, bbyh, bbzl);
    Point ur(0.0, bbyh, bbzh);

    // Horizontal bounds layer. Rays point in +Y, pitch in +Z)
    //
    //   h1  h2
    // ul^   ^ur
    //   |   |
    //   |   |
    // ll+--->lr  pitch = X x h1
    Ray h1(ll, ul);
    Ray h2(lr, ur);
    raypairs.push_back(ray_pair_t(h1, h2));

    // Vertical bounds layer. Rays point in -Z, pitch in +Y)
    //
    // ul<--^--ur  v2
    //      |
    //      | pitch = X x v1
    //      |
    // ll<--+--lr  v1
    Ray v1(lr, ll);
    Ray v2(ur, ul);
    raypairs.push_back(ray_pair_t(v1, v2));
    return raypairs;
}

ray_pair_vector_t WireSchema::ray_pairs(const Store& store, 
                                        const Face& face, 
                                        bool region)
{
    auto raypairs = ray_pairs_active(store, face, region);

    // Each wire plane.
    for (const auto& plane : store.planes(face)) {
        Vector phalf;
        if (region) {
            phalf = 0.5*store.mean_pitch(plane);
        }
        const Wire& w = store.wires()[plane.wires.front()];
        const Ray r1(w.tail - phalf, w.head - phalf);
        const Ray r2(w.tail + phalf, w.head + phalf);
        raypairs.push_back(ray_pair_t(r1, r2));
    }

    return raypairs;
}

