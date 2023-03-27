#include "WireCellUtil/svggpp.hpp"

#include <cassert>
#include <fstream>
#include <iostream>

using namespace svggpp;

int main(int argc, char* argv[])
{
    // fixme: for now, svggpp treats css as a text blob.  It may be
    // extended to represent css also with nlohmann/json.
    auto css = style(R"(
.small {
  font: italic 13px sans-serif;
}
.heavy {
  font: bold 30px sans-serif;
}
.target {
  stroke: green;
  fill: yellow;
}

/* Note that the color of the text is set with the    *
 * fill property, the color property is for HTML only */
.Rrrrr {
  font: italic 40px serif;
  fill: red;
}
)");
    auto c = circle(50,50,30, css_class("target"));
    // in firefox, load file with "#target" appended to URL
    auto v = view("target");
    viewbox(v, view_t{20,20,60,60});
    std::cerr << v.dump() << std::endl;
    auto pg = polygon({ {20,35}, {40,35}, {55,55}, {65,55} }, { {"fill","blue"} });

    auto top = svg(svg_header, {
            css,
            anchor(link(v), {}, c),
            v,
            pg,
            text("My", 20,35,css_class("small")),
            text("cat", 40,35,css_class("heavy")),
            text("is", 55,55,css_class("small")),
            text("Grumpy", 65,55,css_class("Rrrrr"))});
    viewbox(top, view_t{0,0,240,80});
    attrs(top).update(css_class("top"));

    std::cerr << top.dump(4) << std::endl;
    std::cerr << "-----\n";
    std::cout << dumps(top) << std::endl;

    std::string oname = argv[0];
    oname += ".svg";
    std::ofstream out(oname);
    out << dumps(top) << std::endl;
    std::cerr << oname << std::endl;
    return 0;
}
