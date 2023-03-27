#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Math.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "TH1F.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "printer.h"

using namespace WireCell;
using namespace WireCell::String;
using namespace WireCell::WireSchema;


template<typename F>
void scale_one_axis(const std::vector<double>& x, F func)
{
    auto xmm = std::minmax_element(x.begin(), x.end());
    double xmin = *xmm.first, xmax = *xmm.second;
    if (xmin == xmax) {
        xmin -= 0.0001;
        xmin += 0.0001;
    }
    func(xmin, xmax);
}

void scale_both_axes(TGraph* g)
{
    const int n = g->GetN();
    std::vector<double> x(n), y(n);
    for (int i=0; i<n; ++i) {
        x[i] = g->GetPointX(i);
        y[i] = g->GetPointY(i);
    }
    scale_one_axis(x, [&](double min, double max) {
        g->GetXaxis()->SetRangeUser(min, max);
    });
    scale_one_axis(y, [&](double min, double max) {
        g->GetYaxis()->SetRangeUser(min, max);});
}

struct MG {
    Printer& print;
    std::string name;
    std::string cname;
    TGraph* operator()(const std::string& ytit, const std::string& xtit = "WIP index") {
        auto g = new TGraph;
        g->SetTitle((cname+", "+name).c_str());
        g->GetXaxis()->SetTitle(xtit.c_str());
        g->GetYaxis()->SetTitle(ytit.c_str());
        return g;
    }

    void operator()(TGraph* g)
    {
        
        g->Draw("A*");
        // g->GetYaxis()->SetMaxDigits(4);
        // g->GetYaxis()->SetLabelOffset(-0.8);
        scale_both_axes(g);
        // g->GetYaxis()->SetRangeUser(vmin, vmax);
        print.canvas.Modified();
        print.canvas.Update();
        print();
    }
};

void single(Printer& print, const Store& store, const std::string& cname)
{
    const auto& detectors = store.detectors();
    for (const auto& detector : detectors) {
        const auto anodes = store.anodes(detector);
        for (const auto& anode : anodes) {
            const auto faces = store.faces(anode);
            for (const auto& face : faces) {
                const auto planes = store.planes(face);
                for (const auto& plane : planes) {


                    const auto wires = store.wires(plane);
                    std::string name = format("d:%d a:%d f:%d p:%d nw:%d",
                                              detector.ident, anode.ident,
                                              face.ident, plane.ident,
                                              wires.size());
                    std::cerr << name << "\n";
                    MG mg{print, name, cname};

                    auto TX = mg("wire tail X");
                    auto TY = mg("wire tail Y");
                    auto TZ = mg("wire tail Z");
                    auto TYZ = mg("wire tail Z", "wire tail Y");

                    auto HX = mg("wire head X");
                    auto HY = mg("wire head Y");
                    auto HZ = mg("wire head Z");
                    auto HYZ = mg("wire head Z", "wire head Y");

                    auto DX = mg("wire angle with X-axis");
                    auto DY = mg("wire angle with Y-axis");
                    auto DZ = mg("wire angle with Z-axis");

                    auto P = mg("pitch");
                    auto PDX = mg("pitch angle with X-axis");
                    auto PDY = mg("pitch angle with Y-axis");
                    auto PDZ = mg("pitch angle with Z-axis");
                    auto PDYZ = mg("pitch Z direction", "pitch Y direction");

                    Ray rprev;
                    Wire wprev;
                    int pnum=0;

                    for (const auto& wire : wires) {
                        const Ray rnext(wire.tail, wire.head);

                        TX->SetPoint(pnum, pnum, wire.tail.x());
                        TY->SetPoint(pnum, pnum, wire.tail.y());
                        TZ->SetPoint(pnum, pnum, wire.tail.z());
                        TYZ->SetPoint(pnum, wire.tail.y(), wire.tail.z());

                        HX->SetPoint(pnum, pnum, wire.head.x());
                        HY->SetPoint(pnum, pnum, wire.head.y());
                        HZ->SetPoint(pnum, pnum, wire.head.z());
                        HYZ->SetPoint(pnum, wire.head.y(), wire.head.z());

                        const auto wdir = (wire.head - wire.tail).norm();
                        DX->SetPoint(pnum, pnum, 180/pi*acos(wdir.x()));
                        DY->SetPoint(pnum, pnum, 180/pi*acos(wdir.y()));
                        DZ->SetPoint(pnum, pnum, 180/pi*acos(wdir.z()));

                        if (pnum) {
                            const auto pray = ray_pitch(rprev, rnext);
                            const auto pvec = ray_vector(pray);
                            const double pmag = pvec.magnitude()/units::mm;

                            const auto pdir = pvec.norm();
                            const double pdx = 180/pi*acos(pdir.x());
                            const double pdy = 180/pi*acos(pdir.y());
                            const double pdz = 180/pi*acos(pdir.z());

                            const int gnum = pnum-1;
                            P->SetPoint(gnum, gnum, pmag);
                            PDX->SetPoint(gnum, gnum, pdx);
                            PDY->SetPoint(gnum, gnum, pdy);
                            PDZ->SetPoint(gnum, gnum, pdz);
                            PDYZ->SetPoint(gnum, pdir.y(), pdir.z());
                        }
                        rprev = rnext;
                        wprev = wire;
                        ++pnum;
                    } // wires

                    mg(TX);
                    mg(TY);
                    mg(TZ);
                    mg(TYZ);

                    mg(HX);
                    mg(HY);
                    mg(HZ);
                    mg(HYZ);

                    mg(DX);
                    mg(DY);
                    mg(DZ);

                    mg(P);
                    mg(PDX);
                    mg(PDY);
                    mg(PDZ);
                    mg(PDYZ);
                    
                }
            }
        }
    }
}


Vector wires_ray_pitch(const Wire& wa, const Wire& wb)
{
    const Ray ra(wa.tail, wa.head);
    const Ray rb(wb.tail, wb.head);
    return ray_vector(ray_pitch(ra, rb));
}
double wires_ray_angle(const Wire& wa, const Wire& wb)
{
    const Ray ra(wa.tail, wa.head);
    const Ray rb(wb.tail, wb.head);
    return acos(ray_unit(ra).dot(ray_unit(rb)));
}


void couple(Printer& print,
            const Store& a, const Store& b, 
            const std::string& cname)
{
    const auto& detectors = a.detectors();
    for (const auto& detector : detectors) {
        const auto anodes = a.anodes(detector);
        for (const auto& anode : anodes) {
            const auto faces = a.faces(anode);
            for (const auto& face : faces) {
                for (int iplane : face.planes) {
                    const auto aplane = a.planes()[iplane];
                    const auto bplane = b.planes()[iplane];
                    assert(aplane.ident == bplane.ident);
                    const auto plane = aplane;

                    const auto was = a.wires(aplane);
                    const auto wbs = b.wires(bplane);
                    assert(was.size() == wbs.size());

                    std::string name = format("d:%d a:%d f:%d p:%d nw:%d",
                                              detector.ident, anode.ident,
                                              face.ident, plane.ident,
                                              was.size());
                    std::cerr << name << "\n";
                    MG mg{print, name, cname};

                    auto TX = mg("wire tail X change (mm)");
                    auto TY = mg("wire tail Y change (mm)");
                    auto TZ = mg("wire tail Z change (mm)");

                    auto HX = mg("wire head X change (mm)");
                    auto HY = mg("wire head Y change (mm)");
                    auto HZ = mg("wire head Z change (mm)");

                    auto CD = mg("center distance btw old/new (mm)");
                    auto DP = mg("center distance along pitch btw old/new (mm)");
                    auto PD = mg("pitch distance btw old/new (mm)");
                    auto WD = mg("anglular distance btw old/new (deg)");
                    
                    const size_t nwires = was.size();
                    const size_t nhalf = nwires/2;
                    
                    const Vector pmean = wires_ray_pitch(wbs[nhalf], wbs[nhalf+1]).norm();

                    for (size_t wind=0; wind<nwires; ++wind) {
                        const auto& wa = was[wind];
                        const auto& wb = wbs[wind];

                        const auto dtail = wb.tail - wa.tail;
                        const auto dhead = wb.head - wa.head;
                        
                        TX->SetPoint(wind, wind, dtail.x()/units::mm);
                        TY->SetPoint(wind, wind, dtail.y()/units::mm);
                        TZ->SetPoint(wind, wind, dtail.z()/units::mm);
                        HX->SetPoint(wind, wind, dhead.x()/units::mm);
                        HY->SetPoint(wind, wind, dhead.y()/units::mm);
                        HZ->SetPoint(wind, wind, dhead.z()/units::mm);

                        const auto pvec = wires_ray_pitch(wa, wb);
                        const double pmag = pvec.magnitude();
                        const double ang = wires_ray_angle(wa, wb);
                        
                        const auto cenvec = 0.5*((wb.tail+wb.head) - (wa.tail+wa.head));
                        const double dcen = cenvec.magnitude();
                        CD->SetPoint(wind, wind, dcen/units::mm);
                        PD->SetPoint(wind, wind, pmag/units::mm);
                        WD->SetPoint(wind, wind, ang*180/pi);

                        const double dpitch = pmean.dot(cenvec);
                        DP->SetPoint(wind, wind, dpitch/units::mm);


                    } // wires

                    mg(TX);
                    mg(TY);
                    mg(TZ);

                    mg(HX);
                    mg(HY);
                    mg(HZ);

                    mg(CD);
                    mg(DP);
                    mg(PD);
                    mg(WD);
                }
            }
        }
    }
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

    std::vector<Store> stores;
    for (Correction cor=Correction::load; cor<Correction::ncorrections; ++cor) {
        auto store = WireSchema::load(fname.c_str(), (Correction)cor);
        std::cerr << fname << " correction: " << (int)cor << ":\n";
        stores.push_back(store);
        std::cerr << "\n";
    }

    for (size_t ind=0; ind<stores.size(); ++ind) {
        // if (ind !=3) continue;  // debug
        std::string cname = oname + "-single" + std::to_string(ind);
        std::string fname = argv[0];
        fname += "-" + cname;
        Printer print(fname);
        print.canvas.SetLeftMargin(0.3);
        single(print, stores[ind], cname);
    }

    const size_t ordered=1;
    for (size_t ind=ordered+1; ind<stores.size(); ++ind) {
        std::string cname = oname + "-couple" + std::to_string(ordered) + std::to_string(ind);
        std::string fname = argv[0];
        fname += "-" + cname;
        Printer print(fname);
        print.canvas.SetLeftMargin(0.3);
        couple(print, stores[ordered], stores[ind], oname);
    }        

    return 0;

}
