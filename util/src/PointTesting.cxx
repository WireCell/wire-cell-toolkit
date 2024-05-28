#include "WireCellUtil/PointTesting.h"

using namespace WireCell;
using namespace WireCell::PointCloud;

PointCloud::Dataset PointTesting::make_janky_track(const Ray& box, double step_size)
{
    const double line_length = ray_length(box);
    const size_t nsteps = line_length / step_size;
    const auto step_dir = ray_unit(box);

    std::vector<double> x(nsteps,0), y(nsteps,0), z(nsteps,0), q(nsteps,0);
    std::vector<int> uwire_index(nsteps,0), vwire_index(nsteps,0), wwire_index(nsteps,0);

    // Walk the track and place points on grid of step_size to
    // construct the janky sampling
    for (size_t step=0; step<nsteps; ++step) {
        Point pt = box.first + step * step_size * step_dir;

        x[step] = step_size * ( (int) ( 0.5 + ( pt.x() / step_size ) ) );
        y[step] = step_size * ( (int) ( 0.5 + ( pt.y() / step_size ) ) );
        z[step] = step_size * ( (int) ( 0.5 + ( pt.z() / step_size ) ) );
        uwire_index[step] = step*1;
        vwire_index[step] = step*2;
        wwire_index[step] = step*3;

        // Accumulte relative distance away from line to the grid
        // point as some meaningless value to fill in for charge.
        q[step] = std::abs(
            pt.x()/step_size - ( (int) pt.x()/step_size ) +
            pt.y()/step_size - ( (int) pt.y()/step_size ) +
            pt.z()/step_size - ( (int) pt.z()/step_size ) );
    }

    Dataset ds;
    ds.add("x", Array(x));
    ds.add("y", Array(y));
    ds.add("z", Array(z));
    ds.add("q", Array(q));
    ds.add("uwire_index", Array(uwire_index));
    ds.add("vwire_index", Array(vwire_index));
    ds.add("wwire_index", Array(wwire_index));
    return ds;
}

