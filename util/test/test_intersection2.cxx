#include "WireCellUtil/Point.h"
#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/Testing.h"

#include <iostream>
#include <cmath>
#include <map>

using namespace WireCell;
using namespace std;

int main()
{
    D3Vector<int> a(3, 4, 5), b(4, 3, 5), c(-5, -12, -13);
    cout << "a . b : " << a.dot(b) << endl;
    cout << "a x b : " << a.cross(b) << endl;
    cout << "a . b x c : " << a.triplescal(b, c) << endl;
    cout << "a x b x c : " << a.triplevec(b, c) << endl;

    Ray bounds(Point(0, 0, 0), Point(1, 1, 1));
    Vector direction = Point(1, 1, 1).norm();

    for (double x = -1.1; x <= 1; x += 0.5) {
        for (double y = -1.1; y <= 1; y += 0.5) {
            for (double z = -1.0; z <= 1; z += 0.5) {
                Vector point(x, y, z);
                Ray hits(Point(-111, -111, -111), Point(-222, -222, -222));

                int hitmask = box_intersection(bounds, point, direction, hits);
                cerr << "RESULT: " << hitmask << " p=" << point << " hits=" << hits << endl;
                Assert(hitmask >= 0);
            }
        }
    }

    return 0;
}
