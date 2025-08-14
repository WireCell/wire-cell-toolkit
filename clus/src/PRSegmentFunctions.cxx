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


    double segment_track_length(SegmentPtr seg, int flag)
    {
        double length = 0;
        
        if (flag == 1) {
            // Sum dx values from fits (equivalent to original flag==1 case)
            auto& fits = seg->fits();
            for (auto& fit : fits) {
                if (fit.valid() && fit.dx > 0) {
                    length += fit.dx;
                }
            }
        } else {
            // Calculate geometric length from fits (equivalent to original flag==0 case)
            
            const auto& fits = seg->fits();
            if (fits.size() < 2) {
                return 0.0;
            }
            
            for (size_t i = 0; i + 1 < fits.size(); i++) {
                const Point& p1 = fits[i].point;
                const Point& p2 = fits[i + 1].point;
                length += std::sqrt(
                    std::pow(p2.x() - p1.x(), 2) +
                    std::pow(p2.y() - p1.y(), 2) +
                    std::pow(p2.z() - p1.z(), 2)
                );
            }
        }
        
        return length;
    }

    double segment_median_dQ_dx(SegmentPtr seg)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        std::vector<double> vec_dQ_dx;
        vec_dQ_dx.reserve(fits.size());
        
        for (auto& fit : fits) {
            if (fit.valid() && fit.dx > 0 && fit.dQ >= 0) {
                // Add small epsilon to avoid division by zero (same as original)
                vec_dQ_dx.push_back(fit.dQ / (fit.dx + 1e-9));
            }
        }
        
        if (vec_dQ_dx.empty()) {
            return 0.0;
        }
        
        // Use nth_element to find median (same algorithm as original)
        size_t median_index = vec_dQ_dx.size() / 2;
        std::nth_element(vec_dQ_dx.begin(), 
                        vec_dQ_dx.begin() + median_index, 
                        vec_dQ_dx.end());
        
        return vec_dQ_dx[median_index];
    }

    double segment_track_length_threshold(SegmentPtr seg, double threshold)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        double length = 0;
        for (auto& fit : fits) {
            if (fit.valid() && fit.dx > 0 ) {
                double dQ_dx = fit.dQ / (fit.dx + 1e-9); // Add epsilon to avoid division by zero
                if (dQ_dx > threshold || threshold == 0) {
                    length += fit.dx;
                }
            }
        }
        
        return length;
    }


    // double segment_median_dQ_dx(const std::vector<double>& dQ, const std::vector<double>& dx)
    // {
    //     if (dQ.empty() || dx.empty() || dQ.size() != dx.size()) {
    //         return 0.0;
    //     }
        
    //     std::vector<double> vec_dQ_dx;
    //     vec_dQ_dx.reserve(dQ.size());
        
    //     for (size_t i = 0; i < dQ.size(); i++) {
    //         // Add small epsilon to avoid division by zero (same as original)
    //         vec_dQ_dx.push_back(dQ[i] / (dx[i] + 1e-9));
    //     }
        
    //     if (vec_dQ_dx.empty()) {
    //         return 0.0;
    //     }
        
    //     // Use nth_element to find median (same algorithm as original)
    //     size_t median_index = vec_dQ_dx.size() / 2;
    //     std::nth_element(vec_dQ_dx.begin(), 
    //                     vec_dQ_dx.begin() + median_index, 
    //                     vec_dQ_dx.end());
        
    //     return vec_dQ_dx[median_index];
    // }

    // double segment_track_length_threshold(const std::vector<double>& dQ,
    //                                     const std::vector<double>& dx,
    //                                     double threshold)
    // {
    //     if (dQ.empty() || dx.empty() || dQ.size() != dx.size()) {
    //         return 0.0;
    //     }
        
    //     double length = 0;
    //     for (size_t i = 0; i < dQ.size(); i++) {
    //         double dQ_dx = dQ[i] / (dx[i] + 1e-9); // Add epsilon to avoid division by zero
    //         if (dQ_dx > threshold || threshold == 0) {
    //             length += dx[i];
    //         }
    //     }
        
    //     return length;
    // }

    double segment_geometric_length(SegmentPtr seg)
    {
        return segment_track_length(seg, 0); // Always use geometric calculation
    }



}
