// This is NOT a real header.  It's code shared between a few tests in
// an ugly way.  It should be #include'd before main().

#include "TCanvas.h"
#include "TMarker.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TArrow.h"
#include "TH1F.h"
#include "TStyle.h"

#include "printer.h"

void draw_text(const Point& pt, const std::string text, int color = 1, int align = 22);

void dump(const blobs_t& blobs)
{
    info("-----dumping {} blobs:", blobs.size());
    for (const auto& b : blobs) {
        info("\t{}", b.as_string());
    }
    info("------");
}

void draw_point(const Point& p, float size = 1, int style = 20, int color = 1);
void draw_point(const Point& p, float size, int style, int color)
{
    TMarker m;
    m.SetMarkerColor(color);
    m.SetMarkerSize(size);
    m.SetMarkerStyle(style);
    m.DrawMarker(p.z(), p.y());
}
void draw_ray(const Ray& ray, int color = 1, float width = 1.0)
{
    TArrow l;
    l.SetLineColor(color);
    l.SetLineWidth(width);
    l.DrawLine(ray.first.z(), ray.first.y(), ray.second.z(), ray.second.y());
}
void draw_arrow(const Ray& ray, int color = 1, float width = 1.0, float asize = 0, const char* opt = "")
{
    TArrow l;
    l.SetLineColor(color);
    l.SetLineWidth(width);
    l.DrawArrow(ray.first.z(), ray.first.y(), ray.second.z(), ray.second.y(), asize, opt);
}

TH1F* draw_frame(TCanvas& canvas, std::string title)
{
    auto* frame = canvas.DrawFrame(-1.0 * border, -1.0 * border, width + border, height + border);
    frame->SetTitle(title.c_str());
    return frame;
}

const std::vector<int> layer_colors{1, 1, 2, 3, 4};

void draw_strip(const Point& head, const Point& tail, const Vector& raydir, int color, bool outline = true);
void draw_strip(const Point& head, const Point& tail, const Vector& raydir, int color, bool outline)
{
    const double shoot = 2 * std::max(width, height);
    std::vector<Point> points{tail + raydir * shoot, tail - raydir * shoot, head - raydir * shoot,
                              head + raydir * shoot};

    TPolyLine* pl = new TPolyLine;  // like a sieve
    pl->SetLineColor(color);
    pl->SetFillColorAlpha(color, 0.1);
    for (const auto& p : points) {
        pl->SetNextPoint(p.z(), p.y());
    }
    pl->SetNextPoint(points.front().z(), points.front().y());
    pl->Draw("f");
    if (outline) {
        pl->Draw("");
    }
}

void draw_layer(Coordinates& coords, int ilayer, double pitch_mag, const Point& pitch, const Point& center,
                const std::vector<double>& measure)
{
    const Vector ecks(1, 0, 0);
    const auto raydir = ecks.cross(pitch);

    for (size_t ind = 0; ind < measure.size(); ++ind) {
        int color = layer_colors[ilayer];
        if (measure[ind] <= 0.0) {
            continue;
        }
        const auto tail = center + ind * pitch_mag * pitch;
        const auto head = center + (ind + 1) * pitch_mag * pitch;
        draw_strip(tail, head, raydir, color);
    }
}

void draw_strips(Coordinates& coords, const strips_t& strips, bool outline = true);
void draw_strips(Coordinates& coords, const strips_t& strips, bool outline)
{
    const Vector ecks(1, 0, 0);

    for (const auto& strip : strips) {
        int color = layer_colors[strip.layer];
        const auto& pitch = coords.pitch_dirs()[strip.layer];
        const auto raydir = ecks.cross(pitch);
        const auto& center = coords.centers()[strip.layer];

        const double pitch_mag = coords.pitch_mags()[strip.layer];

        const auto pind1 = strip.bounds.first;
        const auto pind2 = strip.bounds.second;
        const double pitch_dist = std::abs(pind1 - pind2) * pitch_mag;

        const auto tail = center + pind1 * pitch_mag * pitch;
        const auto head = tail + pitch_dist * pitch;
        draw_strip(tail, head, raydir, color, outline);
    }
}

Point draw_blob(Coordinates& coords, const Blob& blob, int color = 1)
{
    const auto& corners = blob.corners();
    if (corners.empty()) {
        return Point();
    }

    std::vector<Point> points;
    Point center;
    for (const auto& corn : corners) {
        const auto p = coords.ray_crossing(corn.first, corn.second);
        center += p;
        points.push_back(p);
    }
    center = center * (1.0 / points.size());
    sort(points.begin(), points.end(), [&](const Point& a, const Point& b) {
        const Point ac = a - center;
        const Point bc = b - center;
        const double anga = atan2(ac.y(), ac.z());
        const double angb = atan2(bc.y(), bc.z());
        return anga > angb;
    });

    TPolyLine* pl = new TPolyLine;  // like a sieve
    pl->SetLineColor(color);
    pl->SetLineWidth(2);

    for (const auto& p : points) {
        pl->SetNextPoint(p.z(), p.y());
    }
    pl->SetNextPoint(points.front().z(), points.front().y());
    pl->Draw();
    return center;
}

void draw_points_blobs(Coordinates& coords, Printer& print, const std::vector<Point>& points, const blobs_t& blobs)
{
    int nstrips = 0;
    for (const auto& b : blobs) {
        nstrips += b.strips().size();
    }

    draw_frame(print.canvas, Form("%d points, %d blobs, %d strips", (int) points.size(), (int) blobs.size(), nstrips));
    for (size_t ipt = 0; ipt < points.size(); ++ipt) {
        const auto& p = points[ipt];
        draw_point(p, 1, 24, ipt + 1);
    }
    for (size_t ib = 0; ib < blobs.size(); ++ib) {
        draw_blob(coords, blobs[ib], 1);
    }
}

void draw_points_blobs_solved(Coordinates& coords, Printer& print, const std::vector<Point>& points,
                              const blobs_t& blobs, const std::unordered_map<size_t, float>& blob_charge)
{
    int nstrips = 0;
    for (const auto& b : blobs) {
        nstrips += b.strips().size();
    }

    draw_frame(print.canvas, Form("%d points, %d blobs, %d strips", (int) points.size(), (int) blobs.size(), nstrips));

    for (size_t ipt = 0; ipt < points.size(); ++ipt) {
        const auto& p = points[ipt];
        draw_point(p, 1, 24, ipt + 1);
    }

    for (auto it : blob_charge) {
        size_t ind = it.first;
        const auto& blob = blobs[ind];
        float q = it.second;
        int color = 1;
        if (q < 1.0) {
            color = 2;
        }
        auto center = draw_blob(coords, blob, color);
        if (q < 1.0) {
            draw_text(center, Form("s%d", (int) ind), color, 22);
        }
        else {
            draw_text(center, Form("s%d:%.1f", (int) ind, q), color, 22);
        }
        info("center:{} ind:{} q:[]", center, ind, q);
    }
}

void draw_text(const Point& pt, const std::string text, int color, int align)
{
    TLatex l;
    l.SetTextAlign(align);
    l.SetTextFont(52);
    l.SetTextSize(0.03);
    l.SetTextColor(color);
    l.DrawLatex(pt.z(), pt.y(), text.c_str());
}
void draw_raygrid(Printer& print, const Coordinates& coords,
                  const ray_pair_vector_t& raypairs, TH1F* frame=nullptr);
void draw_raygrid(Printer& print, const Coordinates& coords,
                  const ray_pair_vector_t& raypairs, TH1F* frame)
{
    if (!frame) {
        frame = print.canvas.DrawFrame(-110, 50 - 110, 110, 50 + 110);
        frame->SetTitle("Ray Grid");
    }

    info("got {} ray pairs", raypairs.size());
    std::vector<Point> centers0;
    std::vector<Point> centers1;
    const int ui = 2, vi = 3, wi = 4;
    int index[] = {2, 3, 4};
    int colors[] = {2, 4, 1};
    for (int ind = 0; ind < 3; ++ind) {
        int ri = index[ind];

        auto r0 = raypairs[ri].first;
        auto r1 = raypairs[ri].second;

        draw_ray(r0, colors[ind], 2.0);
        draw_ray(r1, colors[ind], 2.0);
        // auto dir = ray_direction(raypairs[ind].first).norm();
        auto pit = ray_vector(ray_pitch(r0, r1));
        auto cen = 0.5 * (r0.first + r0.second);
        centers0.push_back(cen);
        centers1.push_back(0.5 * (r1.first + r1.second));
        draw_arrow(Ray(cen, cen + pit), colors[ind], 1, 0.01, ">");
    }

    const auto r00 = coords.zero_crossing(ui, vi);
    const auto r0n = coords.ray_crossing({ui, -1}, {vi, 1});
    const auto r1n = coords.ray_crossing({ui, 1}, {vi, -1});

    const auto r01 = coords.ray_crossing({ui, 0}, {vi, 1});
    const auto r10 = coords.ray_crossing({ui, 1}, {vi, 0});

    const int i = 10, j = 16;

    const auto rij = coords.ray_crossing({ui, i}, {vi, j});
    const auto r0j = coords.ray_crossing({ui, 0}, {vi, j});
    const auto ri0 = coords.ray_crossing({ui, i}, {vi, 0});

    const double wpij = coords.pitch_location({ui, i}, {vi, j}, wi);
    const auto wray0 = raypairs[wi].first;
    const auto wray1 = raypairs[wi].second;
    const auto wtail = wray0.first;
    const auto whead = wray0.second;
    const auto wpit = ray_vector(ray_pitch(wray0, wray1));
    const double wpind = wpij / wpit.magnitude();

    draw_point(r00, 1, 24, 1);
    draw_text(r00 + Point(0, 0, -10), "r^{01}_{00}");
    draw_text(r00 + Point(0,  10, 0), "w^{01}", 2);
    draw_text(r00 + Point(0, -10, 0), "w^{10}", 4);

    // draw_point(r01, 1, 24, 2);
    // draw_point(r10, 1, 24, 4);
    draw_point(rij, 1, 20, 1);
    draw_text(rij + Point(0, 0, 10), "r^{01}_{ij}", 9); // purple?

    draw_ray(Ray(r00, centers0[0]), colors[0], 1.0);
    draw_ray(Ray(r00, centers0[1]), colors[1], 1.0);

    draw_ray(Ray(r1n, centers1[0]), colors[0], 1.0);
    draw_ray(Ray(r0n, centers1[1]), colors[1], 1.0);

    draw_text(centers0[0] + Point(0,  10,   0), "p^{0}", colors[0]);
    draw_text(centers0[1] + Point(0, -10,   0), "p^{1}", colors[1]);
    draw_text(centers0[2] + Point(0,   0,  -5), "p^{2}", colors[2]);

    // draw p^lmn_ij
    Point ptail = rij, phead=rij;
    ptail[2] = centers0[2][2];
    ptail[1] += 5;
    phead[1] += 5;
    draw_arrow(Ray(ptail, phead), 9, 2.0, 0.01, ">");
    draw_text(0.5*(phead+ptail) + Point(0,7,0), "p^{012}_{ij}", 9);

    draw_arrow(Ray(r00, r01), colors[0], 2.0, 0.01, ">");
    draw_arrow(Ray(r00, r10), colors[1], 2.0, 0.01, ">");

    draw_arrow(Ray(r00, rij), 1, 1.0, 0.01, ">");
    // draw_arrow(Ray(r00, r0j), colors[0], 2.0, 0.01, ">");
    draw_arrow(Ray(r00, ri0), colors[1], 2.0, 0.01, ">");
    draw_arrow(Ray(ri0, rij), colors[0], 2.0, 0.01, ">");

    draw_text(ri0 + Point(0, 0, -15), "iw^{10}", 4);
    draw_text(0.5 * (ri0 + rij) + Point(0, 0, 20), "jw^{01}", 2);

    Point lower(0,-10,0);
    Point raise(0,+40,0);
    draw_ray(Ray(raise + wtail + wpind * wpit,
                 lower + whead + wpind * wpit), 9, 1.0);

    print();
}
