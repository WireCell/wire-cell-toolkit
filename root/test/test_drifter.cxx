#include "WireCellGen/TrackDepos.h"

#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Testing.h"

#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IRandom.h"

#include "TCanvas.h"
#include "TView.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TColor.h"

using namespace WireCell;
using namespace std;

Gen::TrackDepos make_tracks()
{
    // warning, this bipases some logic in TrackDepos::configure().
    Gen::TrackDepos td;
    td.add_track(10 * units::us, Ray(Point(-10, 0, 0) * units::cm, Point(-100, 0, 10) * units::cm));
    td.add_track(120 * units::us, Ray(Point(-1, 0, 0) * units::cm, Point(-2, -100, 0) * units::cm));
    td.add_track(99 * units::us, Ray(Point(-130, 50, 50) * units::cm, Point(-11, -50, -30) * units::cm));
    return td;
}
IDepo::shared_vector get_depos()
{
    Gen::TrackDepos td = make_tracks();
    IDepo::vector* ret = new IDepo::vector;
    while (true) {
        IDepo::pointer depo;
        bool ok = td(depo);
        if (!ok or !depo) {
            break;
        }
        ret->push_back(depo);
    }
    return IDepo::shared_vector(ret);
}

std::vector<IDepo::pointer> track_depos(double t, const Ray& ray)
{
    Gen::TrackDepos td;
    td.add_track(t, ray);
    std::vector<IDepo::pointer> ret;
    while (true) {
        IDepo::pointer depo;
        bool ok = td(depo);
        ret.push_back(depo);
        if (!ok or !depo) {
            return ret;
        }
    }
}

std::deque<IDepo::pointer> test_depos(const std::vector<IDepo::pointer>& depos, std::string tn = "Drifter");
std::deque<IDepo::pointer> test_depos(const std::vector<IDepo::pointer>& depos, std::string tn)
{
    AssertMsg(depos.back() == nullptr, "Not given EOS terminated stream of depos");
    auto drifter = Factory::find_tn<IDrifter>(tn);
    WireCell::IDrifter::output_queue drifted;
    for (auto depo : depos) {
        WireCell::IDrifter::output_queue outq;
        if (!depo) {
            cerr << "test_depos: sending EOS to drifter\n";
        }
        bool ok = (*drifter)(depo, outq);
        Assert(ok);
        if (outq.size()) {
            cerr << "test_depos: got " << outq.size() << " drifted\n";
        }
        drifted.insert(drifted.end(), outq.begin(), outq.end());
    }
    AssertMsg(drifted.back() == nullptr, "Drifter did not pass on EOS");

    for (auto out : drifted) {
        if (!out) {
            break;
        }
        WireCell::IDepo::vector vec = depo_chain(out);
        AssertMsg(vec.size() > 1, "The history of the drifted deposition is truncated.");
    }

    cerr << "test_depos: start with: " << depos.size() << ", after drifting have: " << drifted.size() << endl;
    return drifted;
}

void test_tracks(std::string tn = "Drifter");
void test_tracks(std::string tn)
{
    std::vector<std::pair<double, double> > endpoints = {
        std::pair<double, double>{-1 * units::cm, 1 * units::cm},
        {20 * units::cm, 21 * units::cm},
        {-20 * units::cm, -21 * units::cm},
        {-10 * units::cm + 3 * units::mm, -10 * units::cm - 3 * units::mm},
        {+10 * units::cm + 3 * units::mm, +10 * units::cm - 3 * units::mm},
        {-2 * units::m + 3 * units::mm, -2 * units::m - 3 * units::mm},
        {+2 * units::m + 3 * units::mm, +2 * units::m - 3 * units::mm}};

    for (auto& ends : endpoints) {
        Ray ray(Point(ends.first, 0, 0), Point(ends.second, 0, 0));
        auto depos = track_depos(0, ray);
        cerr << "test_tracks: " << ray / units::cm << "cm produces " << depos.size() << "depos\n";
        test_depos(depos, tn);
    }
}

void test_time(std::string tn = "Drifter");
void test_time(std::string tn)
{
    Ray ray(Point(1 * units::m, 0, 0), Point(1 * units::m + 1 * units::cm, 0, 0));
    std::vector<IDepo::pointer> alldepos;
    for (double t : {0.0, 10.0 * units::us, 3 * units::ms}) {
        auto depos = track_depos(t, ray);
        Assert(depos.back() == nullptr);
        depos.pop_back();
        cerr << "test_time: " << t / units::us << "us produces " << depos.size() << "depos\n";
        alldepos.insert(alldepos.end(), depos.begin(), depos.end());
    }
    alldepos.push_back(nullptr);
    test_depos(alldepos, tn);
}

void test_order(std::string tn = "Drifter");
void test_order(std::string tn)
{
    const double x1 = 1 * units::m;
    const double x2 = 12 * units::cm;
    const double t1 = 1 * units::us;
    const double t2 = (x1 - x2) / (1.0 * units::mm / units::us) - 3 * units::us;
    const double t3 = 3 * units::ms;

    Ray ray1(Point(x1, 0, 0), Point(x1 + 1 * units::cm, 0, 0));
    Ray ray2(Point(x2, 0, 0), Point(x2 + 1 * units::cm, 0, 0));

    auto depos1 = track_depos(t1, ray1);
    auto depos2 = track_depos(t2, ray2);
    auto depos3 = track_depos(t3, ray1);
    std::vector<IDepo::pointer> alldepos;
    for (auto depos : {depos1, depos2, depos3}) {
        Assert(depos.back() == nullptr);
        depos.pop_back();
        alldepos.insert(alldepos.end(), depos.begin(), depos.end());
    }
    alldepos.push_back(nullptr);
    auto drifted = test_depos(alldepos, tn);
    for (auto& depo : drifted) {
        if (!depo) {
            break;
        }
        cerr << "t=" << depo->time() / units::us << "(" << depo->prior()->time() / units::us << ")"
             << "us, x=" << depo->pos().x() / units::cm << "(" << depo->prior()->pos().x() / units::cm << ")cm\n";
    }
}

IDepo::vector test_drifted(std::string tn)
{
    IDepo::vector result, activity(*get_depos());
    activity.push_back(nullptr);  // EOS

    auto drifter = Factory::find_tn<IDrifter>(tn);

    WireCell::IDrifter::output_queue outq;
    for (auto in : activity) {
        outq.clear();
        if (!in) {
            cerr << "test_drifter: sending EOS to drifter\n";
        }
        bool ok = (*drifter)(in, outq);
        Assert(ok);
        for (auto d : outq) {
            result.push_back(d);
        }
    }
    Assert(!result.back());

    for (auto out : result) {
        if (!out) {
            break;
        }
        WireCell::IDepo::vector vec = depo_chain(out);
        AssertMsg(vec.size() > 1, "The history of the drifted deposition is truncated.");
    }

    cerr << "test_drifter: start with: " << activity.size() << ", after drifting have: " << result.size() << endl;

    return result;
}

Ray make_bbox()
{
    BoundingBox bbox(Ray(Point(-1, -1, -1), Point(1, 1, 1)));
    IDepo::vector activity(*get_depos());
    for (auto depo : activity) {
        bbox(depo->pos());
    }
    Ray bb = bbox.bounds();
    cout << "Bounds: " << bb.first << " --> " << bb.second << endl;
    return bb;
}

int main(int argc, char* argv[])
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");
    {
        auto icfg = Factory::lookup<IConfigurable>("Random");
        auto cfg = icfg->default_configuration();
        icfg->configure(cfg);
    }
    {
        auto icfg = Factory::lookup<IConfigurable>("Drifter");
        auto cfg = icfg->default_configuration();
        // cerr << cfg << endl;
        cfg["drift_speed"] = 1.0 * units::mm / units::us;
        cfg["xregions"][0]["cathode"] = 2 * units::m;
        cfg["xregions"][0]["response"] = 20 * units::cm;
        cfg["xregions"][0]["anode"] = 10 * units::cm;
        cfg["xregions"][1]["anode"] = -10 * units::cm;
        cfg["xregions"][1]["response"] = -20 * units::cm;
        cfg["xregions"][1]["cathode"] = -2 * units::m;
        icfg->configure(cfg);
    }

    test_tracks("Drifter");
    test_time("Drifter");
    test_order("Drifter");

    IDepo::vector drifted = test_drifted("Drifter");

    Ray bb = make_bbox();

    TCanvas c("c", "c", 800, 800);

    TView* view = TView::CreateView(1);
    view->SetRange(bb.first.x(), bb.first.y(), bb.first.z(), bb.second.x(), bb.second.y(), bb.second.z());
    view->ShowAxis();

    // draw raw activity
    IDepo::vector activity(*get_depos());
    TPolyMarker3D orig(activity.size(), 6);
    orig.SetMarkerColor(2);
    int indx = 0;
    for (auto depo : activity) {
        const Point& p = depo->pos();
        orig.SetPoint(indx++, p.x(), p.y(), p.z());
    }
    orig.Draw();

    // draw drifted
    double tmin = -1, tmax = -1;
    for (auto depo : drifted) {
        if (!depo) {
            cerr << "Reached EOI" << endl;
            break;
        }
        auto history = depo_chain(depo);
        Assert(history.size() > 1);

        if (tmin < 0 && tmax < 0) {
            tmin = tmax = depo->time();
            continue;
        }
        tmin = min(tmin, depo->time());
        tmax = max(tmax, depo->time());
    }
    cerr << "Time bounds: " << tmin << " < " << tmax << endl;

    for (auto depo : drifted) {
        if (!depo) {
            cerr << "Reached EOI" << endl;
            break;
        }

        TPolyMarker3D* pm = new TPolyMarker3D(1, 8);
        const Point& p = depo->pos();
        pm->SetPoint(0, p.x(), p.y(), p.z());

        double rel = depo->time() / (tmax - tmin);
        int col = TColor::GetColorPalette(int(rel * TColor::GetNumberOfColors()));
        pm->SetMarkerColor(col);

        pm->Draw();
    }

    c.Print(String::format("%s.pdf", argv[0]).c_str());

    return 0;
}
