#include "WireCellClus/PRSegmentFunctions.h"

namespace WireCell::Clus::PR {


    bool break_segment(Graph& graph, SegmentPtr seg, Point point, double max_dist/*=1e9*/)
    {
        /// sanity checks
        if (! seg->descriptor_valid()) {
            raise<RuntimeError>("break_segment: segment has invalid descriptor\n");
        }
        auto ed = seg->get_descriptor();
        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);
        auto [_, ingraph] = boost::edge(vd1, vd2, graph);
        if (! ingraph) {
            raise<RuntimeError>("break_segment: segment not in graph\n");
        }

        const auto& fits = seg->fits();
        auto itfits = closest_point(fits, point, owp_to_point<Fit>);

        // reject if test point is at begin or end of fits.
        if (itfits == fits.begin() || itfits+1 == fits.end()) {
            return false;
        }

        const auto& wcpts = seg->wcpts();        
        auto itwcpts = closest_point(wcpts, point, owp_to_point<WCPoint>);

        // clamp the wcpts to not be first/last
        if (itwcpts == wcpts.begin()) {
            ++itwcpts;
        }
        else if (itwcpts+1 == wcpts.end()) {
            --itwcpts;
        }

        
        // update graph
        remove_segment(graph, seg);

        auto vtx1 = graph[vd1].vertex;
        auto vtx2 = graph[vd2].vertex;
        auto vtx = make_vertex(graph);

        // WARNING there is no "direction" in the graph.  You can not assume the
        // "source()" of a segment is closest to the segments first point.  As
        // of now, at least...
        auto seg1 = make_segment(graph, vtx, vtx1);
        auto seg2 = make_segment(graph, vtx, vtx2);


        // fill in the new objects.  All three get the middle thing

        seg1->wcpts(std::vector<WCPoint>(wcpts.begin(), itwcpts+1));
        seg2->wcpts(std::vector<WCPoint>(itwcpts, wcpts.end()));
        vtx->wcpt(*itwcpts);

        seg1->fits(std::vector<Fit>(fits.begin(), itfits+1));
        seg2->fits(std::vector<Fit>(itfits, fits.end()));
        vtx->fit(*itfits);

        //.... more for segment
        // dir_weak
        // flags (dir, shower traj, shower topo)
        // particle type and mass and score
        // points clouds
            
        return true;
    }

}
