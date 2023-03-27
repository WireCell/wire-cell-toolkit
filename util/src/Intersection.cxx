#include "WireCellUtil/Intersection.h"

using namespace WireCell;

int WireCell::box_intersection(int axis0, const Ray& bounds, const Point& point, const Vector& dir, Ray& hits)
{
    hits = Ray(point,point);

    int hitmask = 0;
    if (0.0 == dir[axis0]) {
        return hitmask;         // parallel to side
    }

    Point bmin = bounds.first;
    Point bmax = bounds.second;
    for (size_t axis=0; axis<3; ++axis) {
        if (bmin[axis] > bmax[axis]) {
            std::swap(bmin[axis], bmax[axis]);
        }
    }
    
    int axis1 = (axis0 + 1) % 3;
    int axis2 = (axis1 + 1) % 3;

    {  // toward the min intercept
        double intercept = bmin[axis0];
        double scale = (intercept - point[axis0]) / dir[axis0];

        double one = point[axis1] + scale * dir[axis1];
        double two = point[axis2] + scale * dir[axis2];

        if (bmin[axis1] <= one && one <= bmax[axis1] && bmin[axis2] <= two && two <= bmax[axis2]) {
            hitmask |= 1;
            hits.first[axis0] = intercept;
            hits.first[axis1] = one;
            hits.first[axis2] = two;
        }
    }

    {  // toward the max intercept
        double intercept = bmax[axis0];
        double scale = (intercept - point[axis0]) / dir[axis0];

        double one = point[axis1] + scale * dir[axis1];
        double two = point[axis2] + scale * dir[axis2];

        if (bmin[axis1] <= one && one <= bmax[axis1] && bmin[axis2] <= two && two <= bmax[axis2]) {
            hitmask |= 2;
            hits.second[axis0] = intercept;
            hits.second[axis1] = one;
            hits.second[axis2] = two;
        }
    }

    auto hdir = hits.second-hits.first;
    if (hdir.dot(dir) < 0) {
        std::swap(hits.first, hits.second);
        hitmask = ((hitmask&0x1) << 1) | ((hitmask&0x2)>>1);
    }
    return hitmask;
}

int WireCell::box_intersection(const Ray& bounds, const Point& point, const Vector& dir, Ray& hits)
{
    PointSet results;

    // check each projection
    for (int axis = 0; axis < 3; ++axis) {
        Ray res;
        int got = box_intersection(axis, bounds, point, dir, res);

        if (got & 1) {
            // pair<PointSet::iterator, bool> what =
            results.insert(res.first);
        }
        if (got & 2) {
            // pair<PointSet::iterator, bool> what =
            results.insert(res.second);
        }
    }

    hits = Ray(point,point);

    int hitmask = 0;
    int count = 0;
    for(auto& hit : results) {
        if (count == 0) {
            hits.first = hit;
            hitmask |= 1;
        }
        if (count == 1) {
            hits.second = hit;
            hitmask |= 2;
        }
        ++count;
    }

    auto hdir = hits.second-hits.first;
    if (hdir.dot(dir) < 0) {
        std::swap(hits.first, hits.second);
        hitmask = ((hitmask&0x1) << 1) | ((hitmask&0x2)>>1);
    }
    return hitmask;
}
