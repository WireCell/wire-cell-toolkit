#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/Testing.h"

#include <iostream>
#include <random>

using namespace WireCell;
using namespace std;

int main()
{
    {
        Ray bounds(Point(-5, -500, -500), Point(5, 500, 500));
        Ray line(Point(1, 0, -495), Point(1, 1, -495));
        Vector dir = ray_unit(line);
        Ray hits;
        int hitmask = box_intersection(bounds, line.first, dir, hits);
        cerr << "line: " << line << " bounds:" << bounds << endl;
        cerr << "got: " << hitmask << " " << hits << endl;
        Assert(hitmask == 3);
    }

    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<> dist(-2, 2);

    Ray bounds(Point(-1, -1, -1), Point(1, 1, 1));
    for (int ind = 0; ind < 100; ++ind) {
        Point p1(dist(re), dist(re), dist(re));
        Point p2(dist(re), dist(re), dist(re));
        Vector dir = (p2 - p1).norm();

        cout << ind << " point=" << p1 << " dir=" << dir << endl;

        for (int axis = 0; axis < 3; ++axis) {
            Ray hits;
            int hitmask = box_intersection(axis, bounds, p1, dir, hits);
            cout << "\t axis=" << axis << " [" << hitmask << "]hits=" << hits << endl;
            Assert(hitmask >= 0);
        }

        {
            Ray ray(p1, p2), hits;
            int hitmask = box_intersection(bounds, p1, dir, hits);
            cout << "box: hitmask=" << hitmask << " ray=" << ray << " hits=" << hits << endl;
            if (point_contained(p1, bounds) && point_contained(p2, bounds)) {
                AssertMsg(hitmask == 3, "Inside box, but not enough hits");
            }
        }
    }

    return 0;
}
