#include "WireCellIface/IWireSelectors.h"
#include "WireCellGen/WireGenerator.h"
#include "WireCellGen/WireParams.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Persist.h"

#include "TView.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "MultiPdf.h"

#include <boost/range.hpp>

#include <iostream>
#include <cmath>

using namespace WireCell;
using namespace WireCell::Test;
using namespace std;

void test1()
{
    WireParams* params = new WireParams;
    IWireParameters::pointer iwp(params);

    WireGenerator wg;
    IWire::shared_vector wires;
    bool ok = wg(iwp, wires);
    Assert(ok);
    Assert(wires);

    cerr << "Got " << wires->size() << " wires" << endl;
    Assert(wires->size());
    int last_plane = -1;
    int last_index = -1;
    for (auto wit = wires->begin(); wit != wires->end(); ++wit) {
        IWire::pointer wire = *wit;
        int iplane = wire->planeid().index();
        int ident = (1 + iplane) * 100000 + wire->index();
        Assert(ident == wire->ident());

        if (iplane == last_plane) {
            ++last_index;
        }
        else {
            last_plane = iplane;
            last_index = 0;
        }
        Assert(last_index == wire->index());
    }
}

void test2()
{
    double pitches[] = {10.0, 5.0, 3.0, -1};
//    int want[] = {371, 747, 1243, 0};
    int want[] = {375, 747, 1243, 0};
    for (int ind = 0; pitches[ind] > 0.0; ++ind) {
        WireParams* params = new WireParams;
        Configuration cfg = params->default_configuration();
        put(cfg, "pitch_mm.u", pitches[ind]);
        put(cfg, "pitch_mm.v", pitches[ind]);
        put(cfg, "pitch_mm.w", pitches[ind]);
        params->configure(cfg);
        IWireParameters::pointer iwp(params);
        WireGenerator wg;
        IWire::shared_vector wires;
        bool ok = wg(iwp, wires);
        Assert(ok);
        Assert(wires);

        int nwires = wires->size();
        cout << ind << ": pitch=" << pitches[ind] << " nwires=" << nwires << " (want=" << want[ind] << ")" << endl;
        AssertMsg(nwires == want[ind], "Wrong number of wires");
        AssertMsg(Persist::dumps(cfg).size(), "Failed to dump cfg");
    }
}

void test3D(MultiPdf& pdf, bool interactive)
{
    WireParams* params = new WireParams;
    IWireParameters::pointer iwp(params);
    WireGenerator wg;
    IWire::shared_vector wires;
    bool ok = wg(iwp, wires);
    Assert(ok);
    Assert(wires);
    AssertMsg(wires->size(), "Got no wires");

    const Ray& bbox = params->bounds();

    TView* view = TView::CreateView(1);
    view->SetRange(bbox.first.x(), bbox.first.y(), bbox.first.z(), bbox.second.x(), bbox.second.y(), bbox.second.z());
    view->ShowAxis();
    int colors[3] = {2, 4, 1};

    IWire::vector u_wires, v_wires, w_wires;
    copy_if(wires->begin(), wires->end(), back_inserter(u_wires), select_u_wires);
    copy_if(wires->begin(), wires->end(), back_inserter(v_wires), select_v_wires);
    copy_if(wires->begin(), wires->end(), back_inserter(w_wires), select_w_wires);
    size_t n_wires[3] = {u_wires.size(), v_wires.size(), w_wires.size()};

    double max_width = 5;
    for (auto wit = wires->begin(); wit != wires->end(); ++wit) {
        IWire::pointer wire = *wit;
        int iplane = wire->planeid().index();
        int index = wire->index();

        AssertMsg(n_wires[iplane], "Empty plane");
        double width = ((index + 1) * max_width) / n_wires[iplane];

        const Ray ray = wire->ray();
        TPolyLine3D* pl = new TPolyLine3D(2);
        pl->SetPoint(0, ray.first.x(), ray.first.y(), ray.first.z());
        pl->SetPoint(1, ray.second.x(), ray.second.y(), ray.second.z());
        pl->SetLineColor(colors[iplane]);
        pl->SetLineWidth(width);
        pl->Draw();
    }

    pdf();
}

int main(int argc, char** argv)
{
    MultiPdf pdf(argv[0]);
    test1();
    test2();
    test3D(pdf, argc > 1);

    return 0;
}
