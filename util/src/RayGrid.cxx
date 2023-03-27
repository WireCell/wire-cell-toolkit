#include "WireCellUtil/RayGrid.h"
#include "WireCellUtil/Range.h"

using namespace WireCell;
using namespace WireCell::RayGrid;
using namespace WireCell::Range;

Coordinates::~Coordinates()
{
}

Coordinates::Coordinates()
    : m_nlayers(0)
{
}
Coordinates::Coordinates(const Coordinates& other)
{
    (*this) = other;
}
Coordinates& Coordinates::operator=(const Coordinates& other)
{
    m_nlayers = other.m_nlayers;
    m_pitch_mag = other.m_pitch_mag;
    m_pitch_dir = other.m_pitch_dir;
    m_center = other.m_center;
    m_zero_crossing = other.m_zero_crossing;
    m_ray_jump = other.m_ray_jump;
    m_a.resize(boost::extents[m_nlayers][m_nlayers][m_nlayers]);
    m_b.resize(boost::extents[m_nlayers][m_nlayers][m_nlayers]);
    m_a = other.m_a;
    m_b = other.m_b;
    return *this;
}

Coordinates::Coordinates(const ray_pair_vector_t& rays, int normal_axis, double normal_location)
  : m_nlayers(rays.size())
  , m_pitch_mag(m_nlayers, 0.0)
  , m_pitch_dir(m_nlayers)
  , m_center(m_nlayers)
  , m_zero_crossing(m_nlayers, m_nlayers)
  , m_ray_jump(m_nlayers, m_nlayers)
  , m_a(boost::extents[m_nlayers][m_nlayers][m_nlayers])
  , m_b(boost::extents[m_nlayers][m_nlayers][m_nlayers])
{
    // really we are working in 2D space, so project all vectors into the plane.
    auto project = [&](Vector v) {
        v[normal_axis] = normal_location;
        return v;
    };

    // must go through 1, 2 and 3 combinations

    // First, find the per-layer things
    for (layer_index_t ilayer = 0; ilayer < m_nlayers; ++ilayer) {
        const auto& rpair = rays[ilayer];
        const auto& r0 = rpair.first;
        const auto& r1 = rpair.second;

        // vector of closest approach between the two parallel rays
        const auto rpitch = ray_pitch(r0, r1);

        // relative pitch vector
        auto rpv = project(rpitch.second - rpitch.first);

        m_pitch_mag[ilayer] = rpv.magnitude();
        m_pitch_dir[ilayer] = rpv.norm();

        // center point of ray 0
        m_center[ilayer] = 0.5 * (project(r0.first + r0.second));
    }

    // Next find cross-layer things
    for (layer_index_t il = 0; il < m_nlayers; ++il) {
        for (layer_index_t im = 0; im < m_nlayers; ++im) {
            // ray pairs for layer l and m
            const auto& rpl = rays[il];
            const auto& rpm = rays[im];

            const auto& rl0 = rpl.first;
            const auto& rl1 = rpl.second;
            const auto& rm0 = rpm.first;
            const auto& rm1 = rpm.second;

            // Iterate only over triangle of indices to avoid extra work.
            if (il < im) {
                const auto r00 = ray_pitch(rl0, rm0);
                const auto& pl0 = r00.first;
                const auto& pm0 = r00.second;

                // These really should be the same after projection.
                m_zero_crossing(il, im) = project(pl0);
                m_zero_crossing(im, il) = project(pm0);

                // along l-layer ray 0, crossing of m-layer ray 1.
                {
                    const auto ray = ray_pitch(rl0, rm1);
                    const auto jump = project(ray.first - pl0);
                    m_ray_jump(il, im) = jump;
                }
                // along m-layer ray 0, crossing of l-layer ray 1.
                {
                    const auto ray = ray_pitch(rm0, rl1);
                    const auto jump = project(ray.first - pm0);
                    m_ray_jump(im, il) = jump;
                }
            }
            if (il == im) {
                // Initially we set diagonal elements to the "invalid"
                // point.  But then realized these special case values
                // are both useful and meaningful in context.  Setting
                // zero crossing to centers is apparently redundant
                // but m_center must be kept in order that centers()
                // may continue to return a const reference.
                m_zero_crossing(il, im) = m_center[il];
                m_ray_jump(il, im) = ray_unit(rpl.first);
            }
        }
    }

    // Finally, find triple-layer things (coefficients for
    // P^{lmn}_{ij}).  Needs some of the above completed.
    for (layer_index_t in = 0; in < m_nlayers; ++in) {
        const auto& pn = m_pitch_dir[in];
        const double cp = m_center[in].dot(pn);

        for (layer_index_t il = 0; il < m_nlayers; ++il) {
            if (il == in) {
                continue;
            }

            // triangle iteration
            for (layer_index_t im = 0; im < il; ++im) {
                if (im == in) {
                    continue;
                }

                const double rlmpn = m_zero_crossing(il, im).dot(pn);

                const double wlmpn = m_ray_jump(il, im).dot(pn);
                const double wmlpn = m_ray_jump(im, il).dot(pn);

                m_a[il][im][in] = wlmpn;
                m_a[im][il][in] = wmlpn;
                // b is symmetric.
                m_b[il][im][in] = rlmpn - cp;
                m_b[im][il][in] = rlmpn - cp;

                // j*a{lmn} + i*a{mln} + b{lmn}
            }
        }
    }
}

Vector Coordinates::zero_crossing(layer_index_t one, layer_index_t two) const
{
    return m_zero_crossing(one, two);
}

Vector Coordinates::ray_crossing(const coordinate_t& one, const coordinate_t& two) const
{
    const layer_index_t l = one.layer, m = two.layer;
    const auto r00 = m_zero_crossing(l, m);
    const auto& wlm = m_ray_jump(l, m);
    const auto& wml = m_ray_jump(m, l);
    const double i = one.grid, j = two.grid;
    Vector res = r00 + j * wlm + i * wml;
    return res;
}

double Coordinates::pitch_location(const coordinate_t& one, const coordinate_t& two, layer_index_t other) const
{
    const tensor_t::index il = one.layer, im = two.layer, in = other;
    const tensor_t::index i = one.grid, j = two.grid;
    return j * m_a[il][im][in] + i * m_a[im][il][in] + m_b[il][im][in];
}

vector_array1d_t Coordinates::ring_points(const crossings_t& corners) const
{
    const size_t npts = corners.size();
    vector_array1d_t points(npts);
    if (!npts) return points;
    Point avg;
    for (size_t ind : irange(npts)) {
        const auto& [one, two] = corners[ind];
        auto pt = ray_crossing(one, two);
        points[ind] = pt;
        avg += pt;
    }
    avg = avg / (double)npts;
    std::vector<double> ang(npts);
    std::vector<size_t> indices(npts);
    for (size_t ind : irange(npts)) {
        auto dir = points[ind] - avg;
        ang[ind] = atan2(dir[2], dir[1]);
        indices[ind] = ind;
    }
    std::sort(indices.begin(), indices.end(),
              [&](size_t a, size_t b) -> bool {
                  return ang[a] < ang[b];
              });
    
    std::vector<Point> ret(npts);
    for (size_t ind : irange(npts)) {
        ret[ind] = points[indices[ind]];
    }
    return ret;
}
