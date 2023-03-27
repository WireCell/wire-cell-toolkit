#include "WireCellAux/ClusterHelpers.h"

#include "WireCellIface/IChannel.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWire.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/IFrame.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

using namespace WireCell;

using channel_t = cluster_node_t::channel_t;
using wire_t = cluster_node_t::wire_t;
using blob_t = cluster_node_t::blob_t;
using slice_t = cluster_node_t::slice_t;
using meas_t = cluster_node_t::meas_t;


static Json::Value size_stringer(const cluster_node_t& n) { return Json::nullValue; }

static Json::Value jpoint(const Point& p)
{
    Json::Value ret;
    for (int ind = 0; ind < 3; ++ind) {
        ret[ind] = p[ind];
    }
    return ret;
}

static Json::Value channel_jsoner(const cluster_node_t& n)
{
    auto ich = std::get<channel_t>(n.ptr);
    Json::Value ret = Json::objectValue;
    ret["ident"] = ich->ident();
    ret["index"] = ich->index();
    ret["wpid"] = ich->planeid().ident();
    return ret;
}

static Json::Value wire_jsoner(const cluster_node_t& n)
{
    auto iwire = std::get<wire_t>(n.ptr);
    Json::Value ret = Json::objectValue;
    ret["ident"] = iwire->ident();
    ret["index"] = iwire->index();
    ret["wpid"] = iwire->planeid().ident();
    ret["chid"] = iwire->channel();
    ret["seg"] = iwire->segment();
    auto r = iwire->ray();
    ret["tail"] = jpoint(r.first);
    ret["head"] = jpoint(r.second);
    return ret;
}

// this one needs extra info
Json::Value blob_jsoner(const cluster_node_t& n)
{
    auto iblob = std::get<blob_t>(n.ptr);
    IAnodeFace::pointer iface = iblob->face();
    IWirePlane::vector iplanes = iface->planes();
    auto iwire = iplanes[2]->wires()[0];

    auto islice = iblob->slice();

    // set "X" as the slice start time.  Note, this line represents a
    // breaking change but one which I think is for the best.  We used
    // to bake-in a drift speed in order to convert t0 to x0.  This
    // can now fall onto the responsibility of any consumers of the
    // data, if needed.  Better to keep the data clear from such
    // arbitrary transformation.
    double t0 = islice->start();

    const auto& coords = iface->raygrid();

    Json::Value ret = Json::objectValue;
    ret["span"] = islice->span();
    ret["ident"] = iblob->ident();
    ret["value"] = iblob->value();
    ret["error"] = iblob->uncertainty();
    {
        WirePlaneId afid(kUnknownLayer, iface->which(), iface->anode());
        ret["faceid"] = afid.ident();
    }
    ret["sliceid"] = iblob->slice()->ident();
    Json::Value jcorners = Json::arrayValue;
    const auto& blob = iblob->shape();
    for (const auto& c : blob.corners()) {
        Json::Value j = jpoint(coords.ray_crossing(c.first, c.second));
        j[0] = t0;
        jcorners.append(j);
    }
    ret["corners"] = jcorners;

    //-- strip info is redunant with edges connected to wire vertices.
    // Json::Value jstrips = Json::arrayValue;
    // for (const auto& strip : blob.strips()) {
    //     int plane_ind = strip.layer - 2; // fixme
    //     if (plane_ind <0) { continue; }
    //     Json ::Value j = Json::objectValue;
    //     j["wpid"] = iplanes[plane_ind]->planeid().ident();
    //     j["planeind"] = plane_ind;
    //     j["wip1"] = strip.bounds.first;
    //     j["wip2"] = strip.bounds.second;
    //     Json::Value jwires = Json::arrayValue;
    //     const auto& wires = iplanes[plane_ind]->wires();
    //     for (auto wip = strip.bounds.first; wip < strip.bounds.second; ++wip) {
    //         jwires.append(wires[wip]->ident());
    //     }
    //     j["wids"] = jwires;
    //     jstrips.append(j);
    // }
    // ret["strips"] = jstrips;

    return ret;
}

static
Json::Value slice_jsoner(const cluster_node_t& n)
{
    auto islice = std::get<slice_t>(n.ptr);
    Json::Value ret = Json::objectValue;

    ret["ident"] = islice->ident();
    ret["frameid"] = islice->frame()->ident();
    ret["start"] = islice->start();
    ret["span"] = islice->span();

    // note: used to be SOA, now AOS.
    Json::Value jsig = Json::objectValue;
    for (const auto& it : islice->activity()) {
        std::string ind = String::format("%d", it.first->ident());
        jsig[ind]["val"] = it.second.value();
        jsig[ind]["unc"] = it.second.uncertainty();
    }
    ret["signal"] = jsig;

    return ret;
}

static
Json::Value measurement_jsoner(const cluster_node_t& n)
{
    auto imeas = std::get<meas_t>(n.ptr);
    auto sig = imeas->signal();
    Json::Value ret = Json::objectValue;
    ret["ident"] = imeas->ident();
    ret["val"] = sig.value();
    ret["unc"] = sig.uncertainty();
    ret["wpid"] = imeas->planeid().ident();
    // note: used to have chids, which is redundant with the graph
    return ret;
}

Json::Value WireCell::Aux::jsonify(const cluster_graph_t& gr)
{
    std::vector<std::function<Json::Value(const cluster_node_t& ptr)> > jsoners{
        size_stringer, channel_jsoner, wire_jsoner, blob_jsoner, slice_jsoner, measurement_jsoner};
    auto asjson = [&](const cluster_node_t& n) { return jsoners[n.ptr.index()](n); };

    /*
     * - vertices :: a list of graph vertices
     * - edges :: a list of graph edges
     *
     * A vertex is represented as a JSON object with the following attributes
     * - ident :: an indexable ID number for the node, and referred to in "edges"
     * - type :: the letter "code" used in ICluster: one in "sbcwm"
     * - data :: an object holding information about the corresponding vertex object
     *
     * An edge is a pair of vertex ident numbers.
     *
     */

    Json::Value jvertices = Json::arrayValue;
    for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
        const auto& vobj = gr[vtx];
        if (!vobj.ptr.index()) {
            // warn?
            continue;
        }
        Json::Value jvtx = Json::objectValue;
        jvtx["ident"] = (int) vtx;
        jvtx["type"] = String::format("%c", vobj.code());
        jvtx["data"] = asjson(vobj);

        jvertices.append(jvtx);
    }

    Json::Value jedges = Json::arrayValue;
    for (auto eit : boost::make_iterator_range(boost::edges(gr))) {
        Json::Value jedge = Json::arrayValue;
        jedge[0] = (int) boost::source(eit, gr);
        jedge[1] = (int) boost::target(eit, gr);
        jedges.append(jedge);
    }

    Json::Value top = Json::objectValue;
    top["vertices"] = jvertices;
    top["edges"] = jedges;

    return top;
}
