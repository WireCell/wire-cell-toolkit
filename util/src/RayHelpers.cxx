#include "WireCellUtil/RayHelpers.h"

using namespace WireCell;
using namespace WireCell::RayGrid;

ray_pair_vector_t RayGrid::symmetric_raypairs(double width, double height, double pitch_mag, double angle)
{
    ray_pair_vector_t raypairs;

    // corners
    const Vector ll(0, 0, 0), lr(0, 0, width), ul(0, height, 0), ur(0, height, width);

    // wire directions
    const Vector eckx(1, 0, 0);
    const Vector why(0, 1, 0);
    const Vector zee(0, 0, 1);

    // /-wires
    const Vector du(0, cos(angle), sin(angle));
    const Vector pu = eckx.cross(du).norm();
    Vector pjumpu = 0.5 * pitch_mag * pu;
    double mjumpu2 = pjumpu.dot(pjumpu);
    const Ray rayu0(ul + why * mjumpu2 / (why.dot(pjumpu)), ul + zee * mjumpu2 / (zee.dot(pjumpu)));
    pjumpu = 1.5 * pitch_mag * pu;
    mjumpu2 = pjumpu.dot(pjumpu);
    const Ray rayu1(ul + why * mjumpu2 / (why.dot(pjumpu)), ul + zee * mjumpu2 / (zee.dot(pjumpu)));

    // \-wires
    const Vector dv(0, cos(angle), -sin(angle));
    const Vector pv = eckx.cross(dv).norm();
    Vector pjumpv = 0.5 * pitch_mag * pv;
    double mjumpv2 = pjumpv.dot(pjumpv);
    const Ray rayv0(ll + why * mjumpv2 / (why.dot(pjumpv)), ll + zee * mjumpv2 / (zee.dot(pjumpv)));
    pjumpv = 1.5 * pitch_mag * pv;
    mjumpv2 = pjumpv.dot(pjumpv);
    const Ray rayv1(ll + why * mjumpv2 / (why.dot(pjumpv)), ll + zee * mjumpv2 / (zee.dot(pjumpv)));

    // |-wires
    // const Vector dw = why;
    const Vector pw = zee;
    const Vector pjumpw = pitch_mag * pw;
    const Ray rayw0(ll + 0.0 * pjumpw, ul + 0.0 * pjumpw);
    const Ray rayw1(ll + 1.0 * pjumpw, ul + 1.0 * pjumpw);

    // horizontal bounds
    raypairs.push_back(make_pair(Ray(ll, lr), Ray(ul, ur)));

    // vertical bounds
    raypairs.push_back(make_pair(Ray(ll, ul), Ray(lr, ur)));

    // pitch1
    raypairs.push_back(make_pair(rayu0, rayu1));

    // pitch2
    raypairs.push_back(make_pair(rayv0, rayv1));

    // pitch3
    raypairs.push_back(make_pair(rayw0, rayw1));
    return raypairs;
}

std::vector<measure_t> RayGrid::make_measures(const Coordinates& coords, const std::vector<Point>& points)
{
    int nlayers = coords.nlayers();
    std::vector<measure_t> measures(nlayers);
    const auto& pitches = coords.pitch_dirs();
    const auto& centers = coords.centers();
    const auto& pitch_mags = coords.pitch_mags();

    for (size_t ipt = 0; ipt < points.size(); ++ipt) {
        const auto& p = points[ipt];
        for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
            const auto& pit = pitches[ilayer];
            const auto& cen = centers[ilayer];
            const auto rel = p - cen;
            const int pit_ind = pit.dot(rel) / pitch_mags[ilayer];
            if (pit_ind < 0) {
                continue;
            }
            if (ilayer <= 1) {
                if (pit_ind >= 1 or pit_ind < 0) {
                    // debug("mm: pit_ind={} with ipt={}", pit_ind, ipt);
                    // if (pit_ind == 1) {
                    //     debug("\tpit={} cen={} rel={}", pit, cen, rel);
                    // }
                    continue;
                }
            }
            measure_t& m = measures[ilayer];
            if ((int) m.size() <= pit_ind) {
                m.resize(pit_ind + 1, 0.0);
            }

            m[pit_ind] += 1.0;
        }
    }

    return measures;
}

activities_t RayGrid::make_activities(const Coordinates& coords, std::vector<measure_t>& measures)
{
    int nlayers = coords.nlayers();
    activities_t activities;
    for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
        auto& m = measures[ilayer];
        Activity activity(ilayer, {m.begin(), m.end()});
        activities.push_back(activity);
    }
    return activities;
}

