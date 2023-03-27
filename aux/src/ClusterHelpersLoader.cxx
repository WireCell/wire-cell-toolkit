#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/SimpleWire.h"
#include "WireCellAux/SimpleChannel.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleMeasure.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/GraphTools.h"

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::GraphTools;

static Point to_point(Json::Value jpoint)
{
    return Point(jpoint[0].asDouble(),
                 jpoint[1].asDouble(),
                 jpoint[2].asDouble());
}

// Wire.  No dependencies.
static cluster_node_t to_wire(Json::Value jwire)
{
    const int ident = jwire["ident"].asInt();
    const int index = jwire["index"].asInt();
    const auto wpid =  WirePlaneId(jwire["wpid"].asInt());
    const int chid =  jwire["chid"].asInt();
    const int seg  =  jwire["seg"].asInt();
    const auto tail = to_point(jwire["tail"]);
    const auto head = to_point(jwire["head"]);
    const auto ray = Ray(tail, head);
    IWire::pointer iwire = std::make_shared<SimpleWire>(wpid, ident, index, chid, ray, seg);
    return cluster_node_t(iwire);
}

// Channel.  Needs edges to wires.
static cluster_node_t to_channel(Json::Value jchan, const IWire::vector& wires)
{
    const int ident = jchan["ident"].asInt();
    const int index = jchan["index"].asInt();
    IChannel::pointer ich = std::make_shared<SimpleChannel>(ident, index, wires);
    return cluster_node_t(ich);
}

// Blob. Needs edges to slice and wires
static cluster_node_t to_blob(Json::Value jblob,
                              ISlice::pointer islice,
                              IAnodeFace::pointer iface,
                              const IWire::vector& wires)
{
    const int ident = jblob["ident"].asInt();
    const double value = jblob["value"].asDouble();
    const double error = jblob["error"].asDouble();

    std::map<int, RayGrid::Strip> strips;
    for (const auto& iwire : wires) {
        auto wpid = iwire->planeid();
        const int wip = iwire->index();
        const int layer = 2 + wpid.index(); // fixme: how to set bounding box?
        auto sit = strips.find(layer);
        if (sit == strips.end()) {
            strips[layer] = RayGrid::Strip{layer, {wip,wip}};
            continue;
        }
        auto& strip = sit->second;
        strip.bounds.first = std::min(strip.bounds.first, wip);
        strip.bounds.second = std::min(strip.bounds.second, wip+1);
    }

    const RayGrid::Coordinates& coords = iface->raygrid();
    RayGrid::Blob bshape;
    for (const auto& [layer, strip] : strips) {
        bshape.add(coords, strip);
    }

    IBlob::pointer iblob = std::make_shared<SimpleBlob>(ident, value, error, bshape, islice, iface);
    return cluster_node_t(iblob);
}

struct Slice : public ISlice {
    IFrame::pointer m_frame{nullptr};
    int m_ident{0};
    double m_start{0}, m_span{0};
    ISlice::map_t m_activity{};

    Slice(IFrame::pointer frame, int ident, double start, double span, const ISlice::map_t& activity)
        : m_frame(frame), m_ident(ident), m_start(start), m_span(span), m_activity(activity) { }

    virtual ~Slice() {}

    // ISlice interface
    IFrame::pointer frame() const { return m_frame; }
    int ident() const { return m_ident; }
    double start() const { return m_start; }
    double span() const { return m_span; }
    map_t activity() const { return m_activity; }
};
struct Frame : public IFrame {
    virtual ~Frame() {};
    Frame(int ident) : m_ident(ident) {}

    int m_ident{0};
    virtual int ident() const { return m_ident; }
    virtual double time() const { return 0; }
    virtual double tick() const { return 0; }

    virtual ITrace::shared_vector traces() const {
        return std::make_shared<ITrace::vector>();
    }
    virtual Waveform::ChannelMaskMap masks() const {
        return Waveform::ChannelMaskMap{};
    }
    virtual const tag_list_t& frame_tags() const {
        static tag_list_t dummy;
        return dummy;
    }
    virtual const tag_list_t& trace_tags() const {
        static tag_list_t dummy;
        return dummy;
    }
    virtual const trace_list_t& tagged_traces(const tag_t& tag) const {
        static trace_list_t dummy;
        return dummy;
    }
    virtual const trace_summary_t& trace_summary(const tag_t& tag) const {
        static trace_summary_t dummy;
        return dummy;
    }
};

using chid2ichannel_t = std::unordered_map<int, IChannel::pointer>;
using frid2iframe_t = std::unordered_map<int, IFrame::pointer>;
using slid2islice_t = std::unordered_map<int, ISlice::pointer>;

// Slice. Needs edges to channels.  Needs iframe.
static cluster_node_t to_slice(Json::Value jslice,  IFrame::pointer iframe,
                               const chid2ichannel_t& channels)
{
    const int ident = jslice["ident"].asInt();
    const double start = jslice["start"].asDouble();
    const double span = jslice["span"].asDouble();


    ISlice::map_t activity;
    auto jsig = jslice["signal"]; // map from chid as string to val/unc object.
    for (auto jsig_it = jsig.begin(); jsig_it != jsig.end(); ++jsig_it) {
        const auto& jchid = jsig_it.key();
        const int ichid = atoi(jchid.asString().c_str());
        IChannel::pointer ich = nullptr;
        auto chit = channels.find(ichid);
        if (chit != channels.end()) {
            ich = chit->second;
        }
        const auto& jmeas = *jsig_it;
        activity[ich] = ISlice::value_t(jmeas["val"].asDouble(),
                                        jmeas["unc"].asDouble());
    }
    ISlice::pointer islice = std::make_shared<Slice>(iframe, ident, start, span, activity);
    return cluster_node_t(islice);
}

// Measure.  Needs edges to channels
static cluster_node_t to_measure(Json::Value jmeas,
                                 const IChannel::vector& chans)
{
    const int ident = jmeas["ident"].asInt();
    const auto wpid =  WirePlaneId(jmeas["wpid"].asInt());

    const auto sig = IMeasure::value_t(jmeas["val"].asDouble(),
                                       jmeas["unc"].asDouble());

    IMeasure::pointer imeas = std::make_shared<SimpleMeasure>(ident, wpid, sig, chans);
    return cluster_node_t(imeas);
}

ClusterLoader::ClusterLoader(const anodes_t& anodes)
{
    for (const auto& anode : anodes) {
        for (const auto& face : anode->faces()) {
            const WirePlaneId afid(kUnknownLayer, face->which(), face->anode());
            m_faces[afid.ident()] = face;
        }
    }
}

IAnodeFace::pointer ClusterLoader::face(WirePlaneId wpid) const
{
    WirePlaneId afid(kUnknownLayer, wpid.face(), wpid.apa());

    auto it = m_faces.find(afid.ident());
    if (it == m_faces.end()) {
        return nullptr;
    }
    return it->second;
}

cluster_graph_t ClusterLoader::load(const Json::Value& jgraph,
                                    const IFrame::vector& known_frames) const
{
    cluster_graph_t cgraph;

    std::unordered_map<int, cluster_vertex_t> vid2vtx;

    for (const auto& jvtx : jgraph["vertices"]) {
        const int vid = jvtx["ident"].asInt();
        auto vtx = boost::add_vertex(cgraph);
        vid2vtx[vid] = vtx;

        const char code = jvtx["type"].asString()[0];
        const auto& jobj = jvtx["data"];

        if (code == 'w') {      // no dependencies so pack early
            cgraph[vtx] = to_wire(jobj);
            continue;
        }

        if (code == 'b') {
            WirePlaneId afid(jobj["faceid"].asInt());
            if (! face(afid)) {
                THROW(ValueError() << errmsg{"unknown face ident"});
            }
        }
    }

    for (const auto& edge : jgraph["edges"]) {
        boost::add_edge(vid2vtx[edge[0].asInt()],
                        vid2vtx[edge[1].asInt()], cgraph);
    }

    chid2ichannel_t channels;
    for (const auto& jvtx : jgraph["vertices"]) {
        const int vid = jvtx["ident"].asInt();
        const auto vtx = vid2vtx[vid];

        const char code = jvtx["type"].asString()[0];
        const auto& jobj = jvtx["data"];

        if (code == 'c') {
            IWire::vector wires;
            for (auto nnvtx : mir(boost::adjacent_vertices(vtx, cgraph))) {
                const auto& nnobj = cgraph[nnvtx];
                if (nnobj.code() == 'w') {
                    auto iwire = std::get<cluster_node_t::wire_t>(nnobj.ptr);
                    wires.push_back(iwire);
                }
            }
            auto chnode = to_channel(jobj, wires);
            cgraph[vtx] = chnode;
            IChannel::pointer ich = std::get<cluster_node_t::channel_t>(chnode.ptr);
            channels[ich->ident()] = ich;
            continue;
        }
    }

    frid2iframe_t frames;
    for (const auto kf : known_frames) {
        frames[kf->ident()] = kf;
    }
    slid2islice_t slices;
    for (const auto& jvtx : jgraph["vertices"]) {
        const int vid = jvtx["ident"].asInt();
        auto vtx = vid2vtx[vid];

        const char code = jvtx["type"].asString()[0];
        const auto& jobj = jvtx["data"];

        if (code == 's') {
            int frameid = jobj["frameid"].asInt();
            auto iframe = frames[frameid];
            if (!iframe) {
                iframe = std::make_shared<Frame>(frameid);
                frames[frameid] = iframe;
            }
            auto snode = to_slice(jobj, iframe, channels);
            cgraph[vtx] = snode;
            auto islice = std::get<cluster_node_t::slice_t>(snode.ptr);
            slices[islice->ident()] = islice;
            continue;
        }
    }

    // measures
    for (const auto& jvtx : jgraph["vertices"]) {
        const int vid = jvtx["ident"].asInt();
        auto vtx = vid2vtx[vid];

        const char code = jvtx["type"].asString()[0];
        const auto& jobj = jvtx["data"];

        if (code == 'm') {
            IChannel::vector chans;
            for (auto nnvtx : mir(boost::adjacent_vertices(vtx, cgraph))) {
                const auto& nnobj = cgraph[nnvtx];
                if (nnobj.code() == 'c') {
                    auto ichan = std::get<cluster_node_t::channel_t>(nnobj.ptr);
                    chans.push_back(ichan);
                }
            }
            cgraph[vtx] = to_measure(jobj, chans);
            continue;
        }

        if (code == 'b') {
            IWire::vector wires;
            for (auto nnvtx : mir(boost::adjacent_vertices(vtx, cgraph))) {
                const auto& nnobj = cgraph[nnvtx];
                if (nnobj.code() == 'w') {
                    auto iwire = std::get<cluster_node_t::wire_t>(nnobj.ptr);
                    wires.push_back(iwire);
                }
            }
            const int sliceid = jobj["sliceid"].asInt();
            auto islice = slices[sliceid];
            WirePlaneId faceid(jobj["faceid"].asInt());
            cgraph[vtx] = to_blob(jobj, islice, face(faceid), wires);
            continue;
        }
    }


    return cgraph;
}
