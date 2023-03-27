#include "WireCellUtil/Rectangles.h"
#include "WireCellUtil/svg.hpp"

#include <random>
#include <string>
#include <iostream>


using namespace svg;
using namespace WireCell;
using boost::icl::length;
// <Value,XKey,YKey>
using rects_t = Rectangles<char, double, int>;

// scale from [0,1] to cartesian domain axis.
const double scale = 1000;
const double pen = 1.0;

// https://github.com/morhetz/gruvbox
std::vector<std::vector<Color>> colors = {
    { // darker, good for light fg
        Color(0x28,0x28,0x28), // bg
        Color(0xcc,0x24,0x1d), // red
        Color(0x98,0x97,0x1a), // green
        Color(0xd7,0x99,0x21), // yellow
        Color(0x45,0x85,0x88), // blue
        Color(0xb1,0x62,0x86), // purple
        Color(0x68,0x9d,0x6a), // aqua
        Color(0xd6,0x5d,0x0e), // orange
    },
    { // lighter, good for dark fg
        Color(0xeb,0xdb,0xb2), // fg
        Color(0xfb,0x49,0x3f), // red
        Color(0xb8,0xbb,0x26), // green
        Color(0xfa,0xbd,0x2f), // yellow
        Color(0x83,0xa5,0x98), // blue
        Color(0xd3,0x86,0x9b), // purple
        Color(0x8e,0xc0,0x7c), // auqua
        Color(0xfe,0x80,0x19), // orange
    }
};


std::string setstr(const std::set<char>& cs)
{
    std::stringstream letters;
    std::string comma = "";
    letters << "{";
    for (const auto& one : cs) {
        letters << comma << one;
        comma=",";
    }
    letters << "}";
    return letters.str();
}
std::string charstr(char c)
{
    std::stringstream letter;
    letter << "[" << c << "]";
    return letter.str(); 
}

template<typename Key = double>
rects_t::interval_t<Key> line(std::function<double()>& r)
{
    Key a = scale*r(), b = scale*r();
    if (a > b) {
        std::swap(a,b);
    }
    return rects_t::interval_t<Key>(a,b);
}

Document rect2svg(Document& doc,
                  const rects_t::rectangle_t& rect,
                  const std::string& label,
                  const Color& bg,
                  const Color& fg)
{
    const auto& [xi, yi] = rect;

    Group grp;

    Point corner(Point(xi.lower(), yi.lower()));

    grp << Rectangle(corner, length(xi), length(yi), bg);

    grp << Text(corner, label, fg);

    doc << grp;
    return doc;
}


template<typename Region>
Document region2svg(Document& doc,
                    const Region& reg,
                    const Color& bg,
                    const Color& fg)
{
    const auto& [rect, cs] = reg;
    const auto& [xi,yi] = rect;

    Group grp;

    Point corner(Point(xi.lower(), yi.lower()));

    grp << Rectangle(corner, length(xi), length(yi), 
                     Fill(), Stroke(pen, bg));

    Point tpt(Point(xi.lower(), yi.lower() + 0.5*length(yi)));
    grp << Text(tpt, setstr(cs), fg);

    doc << grp;
    return doc;
}

Color randcolor(std::function<double()>& r)
{
    return Color((int)(256*r()),(int)(256*r()),(int)(256*r()));
}
Color randgray(std::function<double()>& r)
{
    int val = 256*r();
    return Color(val,val,val);
}

rects_t doit(int nrecs, std::function<double()> r, Document& doc)
{
    rects_t recs;
    for (int count=0; count < nrecs; ++count) {
        auto rect = rects_t::rectangle_t(line<double>(r), line<int>(r));
        char letter = 'a'+count;
        auto ele = rects_t::element_t(rect, letter);
        std::cerr << "Adding " << letter << ": " << rect.first << " x " << rect.second << std::endl;
        recs += ele;
        auto bg = colors[0][1+count%7];
        auto fg = colors[0][0];
        rect2svg(doc, rect, charstr(letter), bg, fg);
    }
    
    int count = 0;
    for (const auto& reg : recs.regions()) {
        auto bg = colors[1][1+count%7];
        auto fg = colors[1][0];
        region2svg(doc, reg, bg, fg);
        ++count;
    }
    return recs;
}

void doquery(const rects_t& rects,
             const rects_t::xinterval_t& xi,
             const rects_t::yinterval_t& yi,
             Document& doc)
{
    doc << Rectangle(Point(xi.lower(), yi.lower()),
                     length(xi), length(yi),
                     Fill(colors[0][0], 0.1));

    std::cerr << "X<- " << xi << " Y<- " << yi << std::endl;
    int count = 0;
    for (const auto& [rect, qs] : rects.intersection(xi,yi)) {
        const auto& [qxi, qyi] = rect;
        std::cerr << "X-> " << qxi
                  << " Y-> " << qyi
                  << " " << setstr(qs)
                  << std::endl;

        doc << Rectangle(Point(qxi.lower(), qyi.lower()),
                         length(qxi), length(qyi),
                         Fill(colors[1][1+count%7], 0.5));
        ++count;
    }
}


int main(int argc, char* argv[])
{
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<> dist(0,1);
    auto rand = [&]() { return dist(re); };

    Dimensions dimensions(scale, scale);
    std::string fname = argv[0];
    fname += ".svg";
    Document doc(fname, Layout(dimensions, Layout::TopLeft));

    doc << Rectangle(Point(0,0), scale, scale, Color::Silver);
    const double rad = 100;

    doc << Circle(Point(0,0), rad, colors[1][1]);
    doc << Circle(Point(0,scale), rad, colors[1][2]);
    doc << Circle(Point(scale,0), rad, colors[1][3]);
    doc << Circle(Point(scale,scale), rad, colors[1][4]);

    const int nrecs=5;
    auto recs = doit(nrecs, rand, doc);

    const double ks = 0.1*scale;
    auto kxi = rects_t::xinterval_t(0.5*scale-ks, 0.5*scale+ks);
    auto kyi = rects_t::yinterval_t(0.5*scale-ks, 0.5*scale+ks);
    doquery(recs, kxi, kyi, doc);
    doc.save();

    {
        const double scale2 = scale/10;
        const size_t nx = scale2;
        const size_t ny = scale2;
        auto irecs = pixelize(recs,
                              Binning(2*nx, -scale, scale),
                              Binning(2*ny, -scale, scale));
        Dimensions dimensions2(2*scale2, 2*scale2);
        std::string fname2 = argv[0];
        fname2 += "2.svg";
        Document doc2(fname2, Layout(dimensions2, Layout::TopLeft));

        int count = 0;
        for (const auto& ireg : irecs.regions()) {
            auto bg = colors[1][1+count%7];
            auto fg = colors[1][0];
            region2svg(doc2, ireg, bg, fg);
            ++count;
        }
        doc2.save();
    }

    return 0;
}
