#include "WireCellGen/AnodeFace.h"

using namespace WireCell;
using namespace WireCell::Gen;

static ray_pair_vector_t get_raypairs(const BoundingBox& bb, const IWirePlane::vector& planes)
{
    ray_pair_vector_t raypairs;

    const Ray& bbray = bb.bounds();

    // We don't care which end point is which, but conceptually assume
    // point a is "lower" than "b" in Y and Z.  X gets zeroed.
    const Point& a = bbray.first;
    const Point& b = bbray.second;

    // Corners of a box in X=0 plane
    Point ll(0.0, a.y(), a.z());
    Point lr(0.0, a.y(), b.z());
    Point ul(0.0, b.y(), a.z());
    Point ur(0.0, b.y(), b.z());

    // boundary in the horizontal direction (rays point at Y+)
    Ray h1(ll, ul);
    Ray h2(lr, ur);
    raypairs.push_back(ray_pair_t(h1, h2));

    // boundary in the vertical direction (rays point at Z+)
    Ray v1(ll, lr);
    Ray v2(ul, ur);
    raypairs.push_back(ray_pair_t(v1, v2));

    // Now the wire planes.
    for (const auto& plane : planes) {
        const auto& wires = plane->wires();
        const auto wray0 = wires[0]->ray();
        const auto wray1 = wires[1]->ray();
        const auto pitray = ray_pitch(wray0, wray1);
        const auto pitvec = ray_vector(pitray);
        Ray r1(wray0.first - 0.5 * pitvec, wray0.second - 0.5 * pitvec);
        Ray r2(wray0.first + 0.5 * pitvec, wray0.second + 0.5 * pitvec);

        raypairs.push_back(ray_pair_t(r1, r2));
    }

    // for (size_t layer =0; layer<raypairs.size(); ++layer) {
    //     const auto& rp = raypairs[layer];
    //     std::cerr << "L" << layer << " " << planes[0]->planeid() << "\n"
    //               << "\twray0=" << rp.first << "\n"
    //               << "\twray1=" << rp.second << "\n";
    // }

    return raypairs;
}

AnodeFace::AnodeFace(int ident, IWirePlane::vector planes, const BoundingBox& bb, int which, int aid)
  : m_ident(ident)
  , m_planes(planes)
  , m_bb(bb)
  , m_coords(get_raypairs(bb, planes))
  , m_which(which)
  , m_aid(aid)
{
}
AnodeFace::~AnodeFace()
{
}

IWirePlane::pointer AnodeFace::plane(int ident) const
{
    for (auto ptr : m_planes) {
        if (ptr->ident() == ident) {
            return ptr;
        }
    }
    return nullptr;
}
