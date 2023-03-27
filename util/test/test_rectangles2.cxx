#include "WireCellUtil/Rectangles.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace WireCell;

int main(int argc, char* argv[])
{
    using rect_t = Rectangles<char, int>;
    rect_t recs;

    // one way to add rectangle + value data.
    recs.add(0,7, 0, 4, 'a');

    // another way to do it.
    recs += rect_t::element_t(
        rect_t::rectangle_t(rect_t::xinterval_t(3,11),
                            rect_t::yinterval_t(2, 6)),
        'b');

    std::string fname = argv[0];
    fname += ".euk";
    std::ofstream fstr(fname);

    fstr << "frame -1, -1, 12, 7\n";

    std::vector<std::string> colors = {
        "red", "green", "blue", "cyan", "magenta"
    };

    int ind=0;
    for (const auto& [rec, cs] : recs.regions()) {
        const auto& [xi, yi] = rec;
        rect_t::xkey_t x1 = xi.lower();
        rect_t::xkey_t x2 = xi.upper();
        rect_t::ykey_t y1 = yi.lower();
        rect_t::ykey_t y2 = yi.upper();

        fstr << "draw ["
             << "point(" << x1 << "," << y1 << ")."
             << "point(" << x1 << "," << y2 << ")."
             << "point(" << x2 << "," << y2 << ")."
             << "point(" << x2 << "," << y1 << ")."
             << "point(" << x1 << "," << y1 << ")] "
             << colors[ind++] << "\n";
    }
    std::cerr << "eukleides --output=" << argv[0] << ".pdf " << fname << "\n";
    return 0;
}
