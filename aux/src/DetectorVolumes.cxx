#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include <vector>
// #include <iostream>             // only debug

// Implementation is totally local to this comp. unit so no need for namespacing.
class DetectorVolumes;

WIRECELL_FACTORY(DetectorVolumes, DetectorVolumes,
                 WireCell::IDetectorVolumes, WireCell::IConfigurable)


using namespace WireCell;

class DetectorVolumes : public IDetectorVolumes, public IFiducial, public IConfigurable {

public:

    DetectorVolumes() {}
    virtual ~DetectorVolumes() {}

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["anodes"] = Json::arrayValue;
        return cfg;
    }

    // Return a new wpid with layer masked to 0 so it represents only
    // anode+face.  This is the value used to index m_faces.
    WirePlaneId wfid(WirePlaneId wpid) const {
        return WirePlaneId(WirePlaneLayer_t::kAllLayers, wpid.face(), wpid.apa());
    }

    virtual void configure(const Configuration& cfg) {
        m_faces.clear();
        for (const auto& janode : cfg["anodes"]) {
            const std::string anode_tn = janode.asString();
            auto anode = Factory::find_tn<IAnodePlane>(anode_tn);

            for (auto iface : anode->faces()) {
                if (! iface) continue;
                
                auto planes = iface->planes();
                auto wpid = wfid(planes[0]->planeid());
                if (!wpid) {
                    raise<ValueError>("got bogus wpid from anode %s", anode_tn);
                }

                // std::cerr << "DetectorVolumes face for: " << wpid << "\n";
                m_faces[wpid.ident()] = iface;
            }
        }
    }


    // IFiducial
    virtual bool contained(const Point& point) const {
        return contained_by(point);
    }


    // Rest is IDetectorVolumes
    virtual WirePlaneId contained_by(const Point& point) const {
        // This initial imp is perhaps too slow.  There are two options I can
        // think of immediately:
        // 
        // 1) Try to divine a way to represent the BBs on a regular grid of
        // boxes.  Calculate the 3D grid coordinates of a point directly, eg
        // i=floor((x-o)/w), etc for j/y and k/z.  Use ijk to look up iface to
        // do BB.inside() test.
        //
        // 2) Perhaps simpler, construct a k-d tree with BB corners.  Query to
        // find closest corner to point.  Associate corner back to iface and to
        // BB.inside() test.

        for (const auto& [wpident, iface] : m_faces) {
            auto bb = iface->sensitive();
            if (bb.inside(point)) {
                return WirePlaneId(wpident);
            }
        }
        return WirePlaneId(WirePlaneLayer_t::kUnknownLayer, -1, -1);
    }
    
    IAnodeFace::pointer get_face(WirePlaneId wpid) const {
        wpid = wfid(wpid);
        if (!wpid) {
            // std::cerr << "get_face false wpid: " << wpid << std::endl;
            return nullptr; 
        }
        auto it = m_faces.find(wpid.ident());
        if (it == m_faces.end()) {
            // std::cerr << "get_face no face for wpid: " << wpid << std::endl;
            return nullptr;
        }
        return it->second;
    }        

    IWirePlane::pointer get_plane(WirePlaneId wpid) const {
        if (! wpid.valid()) {
            // std::cerr << "get_plane invalid wpid: " << wpid << std::endl;
            return nullptr;
        }
        auto iface = get_face(wpid);
        if (!iface) {
            return nullptr;
        }
        return iface->planes()[wpid.index()];
    }

    virtual int face_dirx(WirePlaneId wpid) const {
        auto iface = get_face(wpid);
        if (!iface) return 0;
        return iface->dirx();
    }

    virtual Vector wire_direction(WirePlaneId wpid) const {
        auto iplane = get_plane(wpid);
        if (! iplane) { return Vector(0,0,0); }
        return iplane->pimpos()->axis(1);
    }
    
    virtual Vector pitch_vector(WirePlaneId wpid) const {
        auto iplane = get_plane(wpid);
        if (! iplane) { return Vector(0,0,0); }
        const auto* pimpos = iplane->pimpos();
        auto pdir = pimpos->axis(2);
        auto pmag = pimpos->region_binning().binsize();
        return pmag*pdir;
    }

private:

    // Map wpid with layer=0 to its face.
    std::map<int, IAnodeFace::pointer> m_faces;


};

