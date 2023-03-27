#include "WireCellUtil/Point.h"
#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/Testing.h"

#include <iostream>
#include <vector>

using namespace WireCell;
using namespace std;

int main()
{
    // note, make this "upside-down" (from max to min) to assure
    // box_intersection() handles mismatch order.
    Ray bounds(Point(10,10,10), Point(0,0,0));

    {
        const Vector dir(0,0,1);
        std::vector<Point> points = {
            Point(5,5,-5),
            Point(5,5,5),
            Point(5,5,15)
        };
        for (const auto& point : points) {
            Ray hits;
            int hm = box_intersection(bounds, point, dir, hits);
            Assert(hm == 3);
            Assert(hits.first.x() == point.x());
            Assert(hits.second.x() == point.x());
            Assert(hits.first.y() == point.y());
            Assert(hits.second.y() == point.y());
            Assert(hits.first.z() == 0);
            Assert(hits.second.z() == 10);
            auto hdir = hits.second-hits.first;
            Assert(hdir.dot(dir) > 0);
        }
    }

    // Put "upside-down" min/max in conventional order
    std::swap(bounds.first, bounds.second);

    {                           // pierce corners
        Vector dir(1,1,1);
        dir = dir.norm();
        Point point(5,5,5);
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        Assert(hm == 3);
        Assert(hits.first == bounds.first);
        Assert(hits.second == bounds.second);
        auto hdir = hits.second-hits.first;
        Assert(hdir.dot(dir) > 0);
    }

    {
        Vector dir(1,1,-1);
        dir = dir.norm();
        {
            Point point(0,0,0);     // degenerate case!
            Ray hits;
            int hm = box_intersection(bounds, point, dir, hits);
            std::cerr << hm << std::endl;
            Assert(hm == 1 || hm == 2); // FP err could change order 
            Assert(hits.first == point);
            Assert(hits.first == hits.second);
            // auto hdir = hits.second-hits.first;
            // Assert(hdir.dot(dir) > 0);
        }
        {
            Point point = bounds.first + dir;
            Ray hits;
            int hm = box_intersection(bounds, point, dir, hits);
            std::cerr << hm << std::endl;
            Assert(hm == 1);
            Assert(hits.second == point);
            auto hdir = hits.second-hits.first;
            Assert(hdir.dot(dir) > 0);
        }
        {
            Point point = bounds.first - dir;
            Ray hits;
            int hm = box_intersection(bounds, point, dir, hits);
            std::cerr << hm << std::endl;
            Assert(hm == 2);
            Assert(hits.first == point);
            auto hdir = hits.second-hits.first;
            Assert(hdir.dot(dir) > 0);
        }
        {
            Point point = bounds.second + dir;
            Ray hits;
            int hm = box_intersection(bounds, point, dir, hits);
            std::cerr << hm << std::endl;
            Assert(hm == 1);
            Assert(hits.second == point);
            auto hdir = hits.second-hits.first;
            Assert(hdir.dot(dir) > 0);
        }
        {
            Point point = bounds.second - dir;
            Ray hits;
            int hm = box_intersection(bounds, point, dir, hits);
            std::cerr << hm << std::endl;
            Assert(hm == 2);
            Assert(hits.first == point);
            auto hdir = hits.second-hits.first;
            Assert(hdir.dot(dir) > 0);
        }
    }

    return 0;
}
