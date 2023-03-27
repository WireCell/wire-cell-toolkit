#include "WireCellUtil/Point.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Math.h"

#include <iomanip>
#include <iostream>
using namespace std;
using namespace WireCell;

int main()
{
    {
        Point ppz;
        Assert(ppz);
        Assert(ppz.x() == 0);
        Assert(ppz.y() == 0);
        Assert(ppz.z() == 0);
    }
    {
        Point pp1;
        pp1.invalidate();
        Assert(!pp1);
        Point pp2 = pp1;
        Assert(!pp2);
    }
    {
        std::vector<Point> pts(3);
        Assert(pts[2]);
        Assert(pts[2].x() == 0);
        Assert(pts[2].y() == 0);
        Assert(pts[2].z() == 0);
        
        std::vector<Point> pts2;
        pts2 = pts;
        Assert(pts2[2]);
        Assert(pts2[2].x() == 0);
        Assert(pts2[2].y() == 0);
        Assert(pts2[2].z() == 0);
    }

    {
        Point pp1(1, -500, -495);
        Point pp2(1, 500, -495);
        PointSet results;
        results.insert(pp1);
        results.insert(pp2);
        AssertMsg(2 == results.size(), "failed to insert");
    }

    {
        Point origin(0, 0, 0);
        Vector zdir(0, 0, 1);
        Point pt(0 * units::mm, 3.92772 * units::mm, 5.34001 * units::mm);
        double dot = zdir.dot(pt - origin);
        cerr << "origin=" << origin / units::mm << ", zdir=" << zdir << ", pt=" << pt / units::mm
             << " dot=" << dot / units::mm << endl;
    }

    Point p1(1, 2, 3);
    Assert(p1.x() == 1 && p1.y() == 2 && p1.z() == 3);

    Point p2 = p1;

    Assert(p1 == p2);

    Point p3(p1);
    Assert(p1 == p3);

    PointF pf(p1);
    Assert(Point(pf) == p1);
    Assert(pf == PointF(p1));

    Point ps = p1 + p2;
    Assert(ps.x() == 2);

    Assert(p1.norm().magnitude() == 1.0);

    double eps = (1 - 1e-11);
    Point peps = p1 * eps;
    cerr << "Epsilon=" << std::setprecision(12) << eps << " peps=" << peps << endl;
    PointSet pset;
    pset.insert(p1);
    pset.insert(p2);
    pset.insert(p3);
    pset.insert(ps);
    pset.insert(peps);
    for (auto pit = pset.begin(); pit != pset.end(); ++pit) {
        cerr << *pit << endl;
    }
    AssertMsg(pset.size() == 2, "tolerance set broken");

    Point pdiff = p1;
    pdiff.set(3, 2, 1);
    AssertMsg(p1 != pdiff, "Copy on write failed.");

    Point foo;
    /// temporarily make this a really big loop to test no memory leakage
    for (int ind = 0; ind < 1000; ++ind) {
        foo = p1;
        foo = p2;
        foo = p3;
    }

    const Ray r1(Point(3.75, -323.316, -500), Point(3.75, 254.034, 500));
    const Ray r2(Point(2.5, -254.034, 500), Point(2.5, 323.316, -500));
    const Ray c12 = ray_pitch(r1, r2);
    cerr << "r1=" << r1 << "\n"
         << "r2=" << r2 << "\n"
         << "rp=" << c12 << "\n";

    BoundingBox bb(pset.begin(), pset.end());
    bb.inside(Point());

    {
        Point p(1, 2, 3);
        Assert(p);
        cerr << "  valid: " << p << endl;
        p.invalidate();
        Assert(!p);
        // cerr << "invalid: " << p << endl;
        p.set();
        Assert(p);
        cerr << "revalid: " << p << endl;
    }
    {
        Ray rprev(Point(6.33552e-13, -1155.1, 10358.5), Point(6.34915e-13, -1148.67, 10369.6));
        Ray rnext(Point(6.34287e-13, -1155.1, 10364.5), Point(6.34915e-13, -1152.13, 10369.6));
        // Ray rprev(Point(0, -1155.1, 10358.5), Point(0, -1148.67, 10369.6));
        // Ray rnext(Point(0, -1155.1, 10364.5), Point(0, -1152.13, 10369.6));
        const auto wpdir = (rprev.second - rprev.first).norm();
        const auto wndir = (rnext.second - rnext.first).norm();
        const auto pray = ray_pitch(rprev, rnext);
        const auto pvec = ray_vector(pray);

        const double pmag = pvec.magnitude();
        const double ang = 180/pi*acos(ray_unit(rprev).dot(ray_unit(rnext)));
        std::cerr << "prev: " << rprev << " " << ray_length(rprev) << " " << ray_unit(rprev) << "\n"
                  << "next: " << rnext << " " << ray_length(rnext) << " " << ray_unit(rnext) << "\n"
                  << "pray: " << pray << " " << ray_length(pray) << "\n"
                  << "pvec: " << pvec << "\n"
                  << "pmag: " << pmag << "\n"
                  << "ang:  " << ang << " deg\n"
                  << "dYp ang: " << 180/pi*acos(wpdir.y()) << "\n"
                  << "dYn ang: " << 180/pi*acos(wndir.y()) << "\n"
            ;
    }

}
