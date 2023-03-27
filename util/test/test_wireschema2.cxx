#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/svg.hpp"

#include <iostream>

using namespace WireCell;
using namespace WireCell::String;

// https://github.com/morhetz/gruvbox
std::vector<std::vector<svg::Color>> colors = {
    { // darker, good for light fg
        svg::Color(0x28,0x28,0x28), // bg
        svg::Color(0xcc,0x24,0x1d), // red
        svg::Color(0x98,0x97,0x1a), // green
        svg::Color(0xd7,0x99,0x21), // yellow
        svg::Color(0x45,0x85,0x88), // blue
        svg::Color(0xb1,0x62,0x86), // purple
        svg::Color(0x68,0x9d,0x6a), // aqua
        svg::Color(0xd6,0x5d,0x0e), // orange
    },
    { // lighter, good for dark fg
        svg::Color(0xeb,0xdb,0xb2), // fg
        svg::Color(0xfb,0x49,0x3f), // red
        svg::Color(0xb8,0xbb,0x26), // green
        svg::Color(0xfa,0xbd,0x2f), // yellow
        svg::Color(0x83,0xa5,0x98), // blue
        svg::Color(0xd3,0x86,0x9b), // purple
        svg::Color(0x8e,0xc0,0x7c), // auqua
        svg::Color(0xfe,0x80,0x19), // orange
    }
};

svg::Point point(const Point& p)
{
    return svg::Point(p.z(), p.y());
}

svg::Line line(const Ray& r,
               svg::Stroke const & stroke = svg::Stroke(0.15, svg::Color::Black))
{
    return svg::Line(point(r.first), point(r.second), stroke);
}

svg::Line line(const WireSchema::Wire& w,
               svg::Stroke const & stroke = svg::Stroke(0.15, svg::Color::Black))
{
    return svg::Line(point(w.tail), point(w.head), stroke);
}


int main(int argc, char* argv[])
{
    std::string fname = "microboone-celltree-wires-v2.1.json.bz2";
    if (argc > 1) {
        fname = argv[1];
    }
    auto slash = fname.rfind("/");
    if (slash == std::string::npos) {
        slash = 0;
    }
    else {
        ++slash;
    }
    const auto ext = fname.rfind(".json");
    const std::string oname = fname.substr(slash, ext-slash);

    auto store = WireSchema::load(fname.c_str());
    const auto anode = store.anodes()[0];

    const size_t wire_skip = 100;
    const double nominal_pitch = 5*units::mm;
    const double border = wire_skip*nominal_pitch;

    const std::vector<svg::Color> plane_colors = {
        svg::Color::Red, svg::Color::Green, svg::Color::Blue
    };
    const std::vector<svg::Color> layer_colors = {
        svg::Color::Orange, svg::Color::Purple,
        svg::Color::Red, svg::Color::Green, svg::Color::Blue
    };

    for (const auto& face : store.faces(anode)) {
        std::string svgname = format("test-wireschema-%s-f%d.svg", oname.c_str(), face.ident);
        std::cerr << svgname << "\n";

        auto bb = store.bounding_box(face).bounds();
        const double bb_width = std::abs(bb.first.z()-bb.second.z());
        const double bb_height = std::abs(bb.first.y()-bb.second.y());
        const double img_width = bb_width + 2*border;
        const double img_height = bb_height + 2*border;
        const svg::Dimensions dims(img_width, img_height);
        const double bb_z0 = std::min(bb.first.z(),bb.second.z());
        const double bb_y0 = std::min(bb.first.y(),bb.second.y());
        const double img_z0 = bb_z0 - border;
        const double img_y0 = bb_y0 - border;
        const svg::Point origin_offset(-img_z0, -img_y0);
        svg::Layout layout(dims, svg::Layout::TopLeft, 1.0, origin_offset);

        std::cerr << "img: z0=" << img_z0 << " y0=" << img_y0 << " w=" << img_width << " h=" << img_height << ")\n";
        std::cerr << " bb: z0=" << bb_z0 << " y0=" << bb_y0 << " w=" << bb_width << " h=" << bb_height << ")\n";

        svg::Document doc(svgname, layout);
        doc << svg::Rectangle(svg::Point(img_z0, img_y0), img_width, img_height, svg::Color::White);
        doc << svg::Rectangle(svg::Point(bb_z0, bb_y0), bb_width, bb_height,
                              svg::Fill(), svg::Stroke(nominal_pitch, svg::Color::Yellow));

        size_t layer_index=0;
        const auto raypairs = ray_pairs(store, face);
        for (const auto& [r1,r2] : raypairs) {
            const auto color = layer_colors[layer_index];
            doc << line(r1, svg::Stroke(0.5*nominal_pitch, color));
            doc << line(r2, svg::Stroke(0.5*nominal_pitch, color));
            ++layer_index;
        }
        
        size_t plane_index=0;
        for (const auto& plane : store.planes(face)) {

            Assert(! plane.wires.empty());
            auto ws = store.wires(plane);
            Assert(ws.size() == plane.wires.size());

            Assert(! ws.empty());
            {
                const size_t nhalf = ws.size()/2;
                Assert(nhalf + 1 < ws.size());
                const auto& w1 = ws[nhalf];
                const auto& w2 = ws[nhalf+1];
                Assert(w1.ident != w2.ident);
            }

            auto pmean = store.mean_pitch(plane);
            auto wmean = store.mean_wire(plane);
            std::cerr << "pmean="<< pmean << " wmean="<<wmean<< "\n";
            Assert(std::abs(pmean.magnitude()) > 1*units::mm); // future detector may be this fine?
            Assert(std::abs(wmean.magnitude()) > 1*units::mm);
            

            auto [W,P] = store.wire_pitch(plane);
            std::cerr << "W=" << W << " P=" << P << "\n";
            Assert(std::abs(W.magnitude() - 1) < 1e-6);
            Assert(std::abs(P.magnitude() - 1) < 1e-6);
            Assert(std::abs((pmean.norm()-P).magnitude()) < 1e-6);

            
            auto w1 = ws.front();
            auto w2 = ws.back();
            // w1.tail.x(0); w1.head.x(0); w2.tail.x(0); w2.head.x(0);

            auto half1 = 0.5*(w1.tail + w1.head);

            const double zh = half1.z();
            const double yh = half1.y();
            // make bullseye around center of first wire
            doc << svg::Circle(svg::Point(zh, yh), wire_skip*pmean,
                               svg::Fill(), svg::Stroke(2*pmean, plane_colors[plane_index]));
            doc << svg::Circle(svg::Point(zh, yh), 2*pmean,
                               plane_colors[plane_index]);

            std::cerr << "face=" << face.ident << " plane=" << plane.ident << "\n";
            std::cerr << "half="<< half1 << "\n";
            Assert (std::abs(pmean.magnitude() - 4*units::mm) < 2*units::mm);

            const size_t nwires = ws.size();
            for (size_t ind=1; ind<nwires; ++ind) {
                if (ind != 1 and ind != nwires-1 and ind % wire_skip) continue;
                auto one = ws[ind-1];
                auto two = ws[ind];

                Vector c = 0.5*(one.tail + one.head);

                Ray r1(one.tail, one.head), r2(two.tail, two.head);
                Ray rp = ray_pitch(r1, r2);
                Vector rpv = ray_vector(rp);
                // double rl = ray_length(rp);

                // Nudge text based on current plane to avoid overlaps
                const auto tnudge = Vector(0,2*pmean*(plane_index+1), 0);
                std::string lab;
                if (ind==1) {
                    lab = format("%d: PID=%d WIP=%d WID=%d pmean=%.1fmm",
                                 plane_index, plane.ident, ind-1, one.ident, pmean);
                }
                else {
                    lab = format("%d: WIP=%d WID=%d y=%.1fmm z=%.1fmm",
                                 plane_index, ind-1, one.ident, c.y(), c.z());
                }
                doc << svg::Text(point(c + tnudge), lab, svg::Color::Black, svg::Font(10, "Verdana"));
                doc << svg::Circle(point(c), pmean, plane_colors[plane_index]);

                doc << line(r1) << line(r2) << line(Ray(c, c+rpv));
            }

            // check central wire against extrapolated ray pairs
            const std::vector<size_t> nchecks = {
                (size_t)(0.25*nwires), (size_t)(0.5*nwires), (size_t)(0.75*nwires), nwires-1
            };
            for (const size_t icheck : nchecks) {
                auto wcen = ws[icheck];
                wcen.tail.x(0); wcen.head.x(0);
                const Point wcpt = 0.5*(wcen.tail + wcen.head); // central wire center

                auto zero = ws[0];
                zero.tail.x(0); zero.head.x(0);
                const Point zcpt  = 0.5*(zero.tail + zero.head); // zero wire center
                const Point mcpt = zcpt + pmean * icheck; // predicted mid center

                const Point tgpt = mcpt + W.dot(wcpt - mcpt) * W; // move along W closer to actual wire center
                const Vector whalf = (nominal_pitch * wire_skip) * W; // pick some longish extent


                // smaller, thinner bullseye 
                doc << svg::Circle(point(tgpt), 0.5*wire_skip*pmean,
                                   svg::Fill(), svg::Stroke(pmean, plane_colors[plane_index]));
                doc << line(wcen, svg::Stroke(pmean, plane_colors[plane_index])); // make easy to find
                doc << line(wcen); // draw precise 
                doc << line(Ray(tgpt + whalf, tgpt - whalf), // check if overlaps
                            svg::Stroke(0.15*units::mm, svg::Color::Magenta));
                std::cerr << "check: " << icheck <<  "\n";
                std::cerr << "\tmcpt=" << mcpt << " whalf=" << whalf << " tgpt=" << tgpt << "\n";
                std::cerr << "\tzcpt=" << zcpt
                          << " =?= 0.5 * (tail=" << zero.tail << " + head=" << zero.head << ")\n";
                std::cerr << "\tmissed by " << (tgpt-wcpt).magnitude()/units::mm << " mm\n";
            }

            std::cerr << "f" << face.ident << " p" << plane.ident
                      << "\n";

            // for (const auto& w : ws) {
            //     std::cerr << "w ["<<w.ident<<","<<w.channel<<","<<w.segment<<"]:"
            //               << w.tail << " --> " << w.head << "\n";
            // }

            ++plane_index;  
                               
        }

                    
        doc.save();
    }
    return 0;
}
