#include "WireCellGen/AnodePlane.h"
#include "WireCellGen/AnodeFace.h"
#include "WireCellGen/WirePlane.h"

#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWireSchema.h"
#include "WireCellIface/SimpleWire.h"
#include "WireCellIface/SimpleChannel.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include <string>

WIRECELL_FACTORY(AnodePlane, WireCell::Gen::AnodePlane,
                 WireCell::INamed,
                 WireCell::IAnodePlane, WireCell::IConfigurable)

using namespace WireCell;
using namespace std;
using WireCell::String::format;

Gen::AnodePlane::AnodePlane()
  : Aux::Logger("AnodePlane", "geom")
{
}

const int default_nimpacts = 10;

struct channel_wire_collector_t {
    // To deal with wrapped wires we need to temporarily hold on to
    // channels in their concrete form, indexec by chanel ident.  This
    // struct explicitly does NOT free its new channels as it is
    // assumed they will all become wrapped in IChannel::pointers.
    std::unordered_map<int, SimpleChannel*> chid2sch;

    // Only call ONCE per iwire.
    void operator()(IWire::pointer iwire)
    {
        const int chid = iwire->channel();
        SimpleChannel* sc = nullptr;
        auto it = chid2sch.find(chid);
        if (it == chid2sch.end()) {
            sc = new SimpleChannel(chid);
            chid2sch[chid] = sc;
        }
        else {
            sc = it->second;
        }
        sc->add(iwire);
    }

    SimpleChannel* operator()(int chid) { return chid2sch[chid]; }
};


template<typename Face>
bool null_face(const Face& face)
{
    // we do not easily support anyOf to allow either a face structure
    // or a null so we assume if all three have default zero value
    // then the user didn't set. 
    return face.anode == 0 and face.response == 0 and face.cathode==0;
}
template<typename Faces>
bool good_faces(const Faces& faces)
{
    if (faces.size() != 2) return false;
    return not null_face(faces[0]) and null_face(faces[1]);
}

void Gen::AnodePlane::configure(const WireCell::Configuration& ccfg)
{
    // Note, normally we would NOT implment configure() and instead
    // leave the work all to the Aux::Configurable<Config> base class.
    // However, here we do some special "live" rewriting of the
    // configuration object in order to provide some backward
    // compatibility to older configuration schema.

    WireCell::Configuration cfg = ccfg;

    // A kludge to continue to support prior configuration promise
    { 
        auto& jfaces = cfg["faces"];
        if (jfaces.isNull() or jfaces.empty() or (jfaces[0].isNull() and jfaces[1].isNull())) {
            log->critical("at least two faces need to be defined, got:\n{}", jfaces);
            THROW(ValueError() << errmsg{"AnodePlane: error in configuration"});
        }
        for (int ind=0; ind<2; ++ind) {
            // allow a null instead of a default face
            if (jfaces[ind].isNull()) {
                jfaces[ind]["response"] = jfaces[ind]["cathode"] = 0;
            }
            // allow anode to be unspecified, taking response
            if (jfaces[ind]["anode"].isNull()) {
                jfaces[ind]["anode"] = jfaces[ind]["response"];
            }
        }
    }
    Aux::Configurable<Config>::configure(cfg);
}

void Gen::AnodePlane::configured()
{
    if (not good_faces(m_cfg.faces)) {
        log->critical("at least one or two faces need to be defined");
        THROW(ValueError() << errmsg{"AnodePlane: error in faces configuration"});
    }

    // get wire schema
    if (m_cfg.wire_schema.empty()) {
        log->critical("\"wire_schema\" parameter must specify an IWireSchema component");
        THROW(ValueError() << errmsg{"\"wire_schema\" parameter must specify an IWireSchema component"});
    }
    auto iws = Factory::find_tn<IWireSchema>(m_cfg.wire_schema);  // throws if not found
    const WireSchema::Store& ws_store = iws->wire_schema_store();

    // keep track which channels we know about in this anode.
    m_channels.clear();

    const WireSchema::Anode& ws_anode = ws_store.anode(ident());

    std::vector<WireSchema::Face> ws_faces = ws_store.faces(ws_anode);
    const size_t nfaces = ws_faces.size();

    channel_wire_collector_t chwcollector;

    m_faces.clear();
    m_faces.resize(nfaces);
    // note, WireSchema requires front/back face ordering in an anode
    for (size_t iface = 0; iface < nfaces; ++iface) {
        const auto& ws_face = ws_faces[iface];
        std::vector<WireSchema::Plane> ws_planes = ws_store.planes(ws_face);
        const size_t nplanes = ws_planes.size();

        // location of imaginary boundary planes
        bool sensitive_face = true;
        const auto& face = m_cfg.faces[iface];
        if (null_face(face)) {
            sensitive_face = false;
            log->debug("anode {} face {} is not sensitive", ident(), iface);
        }
        log->debug("X planes: \"cathode\"@ {}m, \"response\"@{}m, \"anode\"@{}m",
                   face.cathode / units::m,
                   face.response / units::m,
                   face.anode / units::m);

        IWirePlane::vector planes(nplanes);
        // note, WireSchema requires U/V/W plane ordering in a face.
        std::vector<BoundingBox> bbvols;  // sensitive volume for all planes
        for (size_t iplane = 0; iplane < nplanes; ++iplane) {
            const auto& ws_plane = ws_planes[iplane];

            WirePlaneId wire_plane_id(iplane2layer[iplane], iface, ident());
            if (wire_plane_id.index() < 0) {
                log->critical("Bad wire plane id: {}", wire_plane_id.ident());
                THROW(ValueError() << errmsg{format("bad wire plane id: %d", wire_plane_id.ident())});
            }

            // (wire,pitch) directions
            auto wire_pitch_dirs = ws_store.wire_pitch(ws_plane);
            auto ecks_dir = wire_pitch_dirs.first.cross(wire_pitch_dirs.second);

            std::vector<int> plane_chans = ws_store.channels(ws_plane);
            m_channels.insert(m_channels.end(), plane_chans.begin(), plane_chans.end());

            std::vector<WireSchema::Wire> ws_wires = ws_store.wires(ws_plane);
            const size_t nwires = ws_wires.size();
            IWire::vector wires(nwires);

            // note, WireSchema requires wire-in-plane ordering
            for (size_t iwire = 0; iwire < nwires; ++iwire) {
                const auto& ws_wire = ws_wires[iwire];
                const int chanid = plane_chans[iwire];

                Ray ray(ws_wire.tail, ws_wire.head);
                auto iwireptr =
                    make_shared<SimpleWire>(wire_plane_id, ws_wire.ident, iwire, chanid, ray, ws_wire.segment);
                wires[iwire] = iwireptr;
                m_c2wires[chanid].push_back(iwireptr);
                m_c2wpid[chanid] = wire_plane_id.ident();
                chwcollector(iwireptr);
            }  // wire

            const BoundingBox bb = ws_store.bounding_box(ws_plane);
            const Ray bb_ray = bb.bounds();
            const Vector plane_center = 0.5 * (bb_ray.first + bb_ray.second);

            const double pitchmin = wire_pitch_dirs.second.dot(wires[0]->center() - plane_center);
            const double pitchmax = wire_pitch_dirs.second.dot(wires[nwires - 1]->center() - plane_center);
            const Vector pimpos_origin(face.response, plane_center.y(), plane_center.z());

            log->debug("face:{}, plane:{}, origin:{} mm", iface, iplane, pimpos_origin / units::mm);

            Pimpos* pimpos = new Pimpos(nwires, pitchmin, pitchmax, wire_pitch_dirs.first, wire_pitch_dirs.second,
                                        pimpos_origin, m_cfg.nimpacts);

            IChannel::vector plane_channels;
            {
                for (auto w : wires) {
                    if (w->segment() > 0) {
                        continue;
                    }

                    const int chanid = w->channel();
                    SimpleChannel* sch = chwcollector(chanid);
                    sch->set_index(plane_channels.size());
                    IChannel::pointer ich(sch);
                    m_ichannels[chanid] = ich;
                    plane_channels.push_back(ich);
                }
            }
            sort(wires.begin(), wires.end(), IWireCompareIndex());  // redundant?
            planes[iplane] = make_shared<WirePlane>(ws_plane.ident, pimpos, wires, plane_channels);

            // compute sensitive volume for this plane
            const double mean_pitch = (pitchmax - pitchmin) / (nwires - 1);
            if (sensitive_face) {
                auto v1 = bb_ray.first;
                auto v2 = bb_ray.second;
                BoundingBox bb;
                Point pext;
                // Enlarge to anode/cathode planes in X and by 1/2 pitch in Y or Z.
                auto pitch_dir = wire_pitch_dirs.second;
                if (std::fabs(pitch_dir.z()) > 0.999)      // parallel to z. FIXME: probably need a better comparison
                    pext = Point(0, 0, 0.5 * mean_pitch);  // 1/2 pitch extension in Z
                else if (std::fabs(pitch_dir.y()) > 0.999)
                    pext = Point(0, 0.5 * mean_pitch, 0);

                Point p1(face.anode, std::min(v1.y(), v2.y()), std::min(v1.z(), v2.z()));
                p1 += pext * (-1.0);
                bb(p1);
                Point p2(face.cathode, std::max(v1.y(), v2.y()), std::max(v1.z(), v2.z()));
                p2 += pext;
                bb(p2);
                bbvols.push_back(bb);
            }
        }  // plane

        // intersection of sensitive volumes from all planes
        BoundingBox sensvol;
        if (not bbvols.empty()) {
            Ray r1 = bbvols.at(0).bounds();
            for (auto& bb : bbvols) {
                r1 = WireCell::box_intersect(r1, bb.bounds());
            }
            sensvol(r1);
        }

        log->debug("face:{} with {} planes and sensvol: {}",
                   ws_face.ident, planes.size(), sensvol.bounds());

        m_faces[iface] = make_shared<AnodeFace>(ws_face.ident, planes, sensvol);
    }  // face

    // remove any duplicate channels
    std::sort(m_channels.begin(), m_channels.end());
    auto chend = std::unique(m_channels.begin(), m_channels.end());
    m_channels.resize(std::distance(m_channels.begin(), chend));
}

IAnodeFace::pointer Gen::AnodePlane::face(int ident) const
{
    for (auto ptr : m_faces) {
        if (ptr->ident() == ident) {
            return ptr;
        }
    }
    return nullptr;
}

WirePlaneId Gen::AnodePlane::resolve(int channel) const
{
    const WirePlaneId bogus(0xFFFFFFFF);

    auto got = m_c2wpid.find(channel);
    if (got == m_c2wpid.end()) {
        return bogus;
    }
    return WirePlaneId(got->second);
}

IChannel::pointer Gen::AnodePlane::channel(int chident) const
{
    auto it = m_ichannels.find(chident);
    if (it == m_ichannels.end()) {
        return nullptr;
    }
    return it->second;
}

std::vector<int> Gen::AnodePlane::channels() const { return m_channels; }

IWire::vector Gen::AnodePlane::wires(int channel) const
{
    auto it = m_c2wires.find(channel);
    if (it == m_c2wires.end()) {
        return IWire::vector();
    }
    return it->second;
}
