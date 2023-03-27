/*
This test does this:

     0  3   7   11
     |  |   |   |
                                              
0 ─  ┌──────┐                  ┌──────┐
     │ a    │                  └──────┘a
2 ─  │  ┌───┼───┐              ┌─┐┌───┐┌──┐
     │  │   │   │   disjoin    │a││a,b││b │
4 ─  └──┼───┘   │  ─────────>  └─┘└───┘└──┘
        │ b     │                 ┌───────┐
6 ─     └───────┘                 └───────┘b

With this:
*/
#include "boost/icl/interval_map.hpp"
/*
Except I got it transposed so it is Y-major instead of X-major.

Main thing is to confirm we may nest.

Thank you herbstluftwm(1) for the drawing.
*/

#include <iostream>
#include <fstream>
#include <string>

using key_t = int;                  // coordinates
typedef boost::icl::interval<key_t>::interval_type interval_t;

using set_t = std::set<char>; // payload

using ymap_t = boost::icl::interval_map<key_t, set_t>;
using xmap_t = boost::icl::interval_map<key_t, ymap_t>;

int main(int argc, char* argv[])
{
    xmap_t xmap;

    {
        ymap_t ymap;
        ymap += std::make_pair(interval_t::right_open(0,  4), set_t{'a'});
        xmap += std::make_pair(interval_t::right_open(0,  7), ymap);
    }
    {
        ymap_t ymap;        
        ymap += std::make_pair(interval_t::right_open(2,  6), set_t{'b'});
        xmap += std::make_pair(interval_t::right_open(3, 11), ymap);
    }
    

    std::string fname = argv[0];
    fname += ".euk";
    std::ofstream fstr(fname);

    fstr << "frame -1, -1, 12, 7\n";

    std::vector<std::string> colors = {
        "red", "green", "blue", "cyan", "magenta"
    };

    int ind=0;
    for (const auto& [xi, ymap] : xmap) {
        key_t x1 = xi.lower();
        key_t x2 = xi.upper();

        for (const auto& [yi, cs] : ymap) {
            key_t y1 = yi.lower();
            key_t y2 = yi.upper();

            fstr << "draw ["
                 << "point(" << x1 << "," << y1 << ")."
                 << "point(" << x1 << "," << y2 << ")."
                 << "point(" << x2 << "," << y2 << ")."
                 << "point(" << x2 << "," << y1 << ")."
                 << "point(" << x1 << "," << y1 << ")] "
                 << colors[ind++] << "\n";
        }
    }
    std::cerr << "eukleides --output=" << argv[0] << ".pdf " << fname << "\n";
    return 0;
}
