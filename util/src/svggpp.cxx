#include "WireCellUtil/svggpp.hpp"
#include <stdexcept>

#include <regex>

// tag: required SVG element name
const std::string tag_key = "tag";

// attr: optional array of objects providing attributes as key/value, values are string.
const std::string attr_key = "attr";

// body: optional element body.
const std::string body_key = "body";

static
std::vector<std::string> split(const std::string str,
                               const std::string regex_str) {
    std::regex regexz(regex_str);
    return {std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
            std::sregex_token_iterator()};
}

static
std::string join(const std::vector<std::string>& parts, const std::string& delim=" ")
{
    std::stringstream ss;
    std::string comma="";
    for (const auto& one : parts) {
        ss << comma << one;
        comma = delim;
    }
    return ss.str();
}

using namespace svggpp;

xml_t& svggpp::attrs(xml_t& elem)
{
    if (elem.contains(attr_key)) {
        return elem[attr_key];
    }
    return elem;
}
const xml_t& svggpp::attrs(const xml_t& elem)
{
    if (elem.contains(attr_key)) {
        return elem[attr_key];
    }
    return elem;
}

xml_t& svggpp::body(xml_t& elem)
{
    if (elem.contains(body_key)) {
        return elem[body_key];
    }
    return elem;
}
const xml_t& body(const xml_t& elem)
{
    if (elem.contains(body_key)) {
        return elem[body_key];
    }
    return elem;
}

bool svggpp::has(const xml_t& xml, const std::string& key)
{
    if (xml.contains(tag_key)) {
        return has(attrs(xml), key);
    }
    return xml.contains(key);
}

std::string svggpp::get(const xml_t& xml, const std::string& key, const std::string& def)
{
    if (xml.contains(tag_key)) {
        return get(attrs(xml), key, def);
    }
    if (xml.contains(key)) {
        return xml[key];
    }
    return def;
}
std::string svggpp::get(const xml_t& xml, const std::string& key)
{
    if (xml.contains(tag_key)) {
        return get(attrs(xml), key);
    }
    if (xml.contains(key)) {
        return xml[key];
    }
    throw std::invalid_argument("no key in object " + key);
}

void svggpp::set(xml_t& xml, const std::string& key, const std::string& value)
{
    if (xml.contains(tag_key)) {
        set(attrs(xml), key, value);
        return;
    }
    xml[key] = value;
}

xml_t svggpp::subset(const xml_t& xml, const std::string& key)
{
    std::vector<std::string> keys = {key};
    return subset(xml, keys);
}

xml_t svggpp::subset(const xml_t& xml, const std::vector<std::string>& keys)
{
    xml_t ret;
    for (const auto& key : keys) {
        if (has(xml, key)) {
            ret[key] = get(xml,key);
        }
    }

    return ret;
}

std::string svggpp::id(const xml_t& xml)
{
    return get(xml, "id", "");
}
void svggpp::id(xml_t& xml, const std::string& name)
{
    set(xml, "id", name);
}

std::string svggpp::link(const xml_t& xml)
{
    return "#" + id(xml);
}


xml_t svggpp::point(double x, double y)
{
    return xml_t::object({
        {"x", std::to_string(x)},
        {"y", std::to_string(y)}});
}

xml_t svggpp::point(const point_t& pt) {
    return point(pt.first, pt.second);
}

xml_t svggpp::points(const std::vector<point_t>& pts)
{
    std::stringstream ss;
    std::string space="";
    for (const auto& p : pts) {
        ss << space << p.first << "," << p.second;
        space = " ";
    }
    return xml_t::object({ {"points", ss.str()} });
}

xml_t svggpp::center(double cx, double cy)
{
    return xml_t::object({
        {"cx", std::to_string(cx)},
        {"cy", std::to_string(cy)}});
}

xml_t svggpp::size(double width, double height)
{
    return xml_t::object({
        {"width", std::to_string(width)},
        {"height", std::to_string(height)}});
}


xml_t svggpp::box(double x, double y, double width, double height)
{
    xml_t xml = point(x,y);
    xml.update(size(width,height));
    return xml;
}

xml_t svggpp::endpoints(double x1, double y1, double x2, double y2)
{
    return xml_t::object({
        {"x1", std::to_string(x1)},
        {"y1", std::to_string(y1)},
        {"x2", std::to_string(x2)},
        {"y2", std::to_string(y2)}});
}

view_t svggpp::viewbox(const xml_t& xml)
{
    if (xml.contains(tag_key)) {
        return viewbox(attrs(xml));
    }

    std::string vb = xml["viewBox"];
    auto parts = split(vb, " ");

    view_t v;
    v.x = cast<double>(parts[0]);
    v.y = cast<double>(parts[1]);
    v.width = cast<double>(parts[2]);
    v.height = cast<double>(parts[3]);
    return v;
}

void svggpp::viewbox(xml_t& xml, const view_t& view)
{
    if (xml.contains(tag_key)) {
        viewbox(attrs(xml), view);
        return;
    }
    std::stringstream ss;
    ss << view.x << " " << view.y << " " << view.width << " " << view.height;
    xml["viewBox"] = ss.str();
}

view_t svggpp::viewport(const xml_t& xml)
{
    if (xml.contains(tag_key)) {
        return viewport(attrs(xml));
    }
    view_t v;
    if (xml.contains("x")) {
        v.x = cast<double>(xml["x"]);
    }
    if (xml.contains("y")) {
        v.y = cast<double>(xml["y"]);
    }
    if (xml.contains("width")) {
        v.width = cast<double>(xml["width"]);
    }
    if (xml.contains("height")) {
        v.height = cast<double>(xml["height"]);
    }
    return v;
}

void svggpp::viewport(xml_t& xml, const view_t& view)
{
    if (xml.contains(tag_key)) {
        viewport(attrs(xml), view);
        return;
    }

    xml["x"] = std::to_string(view.x);
    xml["y"] = std::to_string(view.y);
    xml["width"] = std::to_string(view.width);
    xml["height"] = std::to_string(view.height);
}


xml_t svggpp::css_style(const std::string& css)
{
    return xml_t::object({ {"style",css} });
}

xml_t svggpp::css_class(const std::string& cname)
{
    return xml_t::object({ {"class",cname} });
}

xml_t svggpp::css_class(const std::vector<std::string>& cnames)
{
    return css_class(join(cnames));
}



//
// Eleme
// 

static void assert_is_attr(const svggpp::xml_t& one)
{
    for (auto& [key, val] : one.items()) {
        if (! val.is_string()) {
            throw std::invalid_argument("illegal svg attribute value type");
        }
    }
}
static void assert_is_elem(const svggpp::xml_t& one)
{
    if (one.empty()) {
        throw std::invalid_argument("element is null");
    }
    if (one.is_array()) {
        throw std::invalid_argument("element is array");
    }
    if (one.is_object()) {
        if (! one.contains(tag_key)) {
            throw std::invalid_argument("element has no tag");
        }
        if (! one[tag_key].is_string()) {
            throw std::invalid_argument("element tag is not a string");
        }
        if (one.contains(attr_key)) {
            if (! one[attr_key].is_object()) {
                throw std::invalid_argument("element attributes not an object");            
            }
        }
        if (one.contains(body_key)) {
            if (! one[body_key].is_array()) {
                throw std::invalid_argument("element body not an array");            
            }
        }
        return;
    }
    // all other types okay (string, number, bool)
}

svggpp::xml_t svggpp::element(const std::string& tag,
                               const svggpp::xml_t& cattrs,
                               const svggpp::xml_t& cbody)
{
    xml_t attrs = cattrs;
    xml_t body = cbody;
    xml_t elem = {{tag_key,tag},
                   {attr_key, xml_t::object()},
                   {body_key, xml_t::array()}};

    // attrs
    if (attrs.empty()) { attrs = xml_t::object(); }
    if (attrs.is_object()) {
        attrs = xml_t::array({attrs});
    }
    if (attrs.is_array()) {
        for (auto& one : attrs) {
            assert_is_attr(one);
            elem[attr_key].update(one);
        }
    }

    // body
    if (body.empty()) { body = xml_t::array(); }
    if (! body.is_array()) {
        body = xml_t::array({body});
    }
    for (auto& one : body) {
        if (one.empty()) {
            continue;
        }
        if (one.is_object()) {
            assert_is_elem(one);
        }
        elem[body_key].push_back(one);
    }

    return elem;
}

// xml_t& svggpp::update(xml_t& elem, xml_t attrs)
// {
//     elem[attr_key].update(attrs);
//     return elem;
// }

xml_t& svggpp::append(xml_t& parent, const xml_t& cchild)
{
    xml_t child = cchild;
    if (! child.is_array()) {
        child = xml_t::array({child});
    }
    for (auto& c : child) {
        parent[body_key].push_back(c);
    }
    return parent;
}

xml_t svggpp::text(const std::string& t, double x, double y, xml_t attrs)
{
    attrs.update(svggpp::point(x,y));
    return element("text", attrs, t);
}

xml_t svggpp::circle(double cx, double cy, double r, xml_t attrs)
{
    attrs.update(svggpp::center(cx,cy));
    attrs.update({ {"r", std::to_string(r) } });
    return element("circle", attrs);
}

xml_t svggpp::circle(point_t cen, double r, xml_t attrs)
{
    return circle(cen.first, cen.second, r, attrs);
}

xml_t svggpp::line(double x1, double y1, double x2, double y2,
                    xml_t attrs)
{
    attrs.update(svggpp::endpoints(x1,y1,x2,y2));
    return element("line", attrs);
}

xml_t svggpp::rect(double x, double y, double width, double height,
                    xml_t attrs)
{
    attrs.update(svggpp::point(x,y));
    attrs.update(svggpp::size(width,height));
    return element("rect", attrs);
}

xml_t svggpp::polygon(const std::vector<point_t>& pts, xml_t attrs)
{
    attrs.update(svggpp::points(pts));
    return element("polygon", attrs);
}

xml_t svggpp::style(const std::string& css)
{
    return element("style", nullptr, css);
}

xml_t svggpp::title(const std::string& tit,
                    double x, double y, double width, double height,
                    xml_t attrs)
{
    attrs.update(svggpp::point(x,y));
    attrs.update(svggpp::size(width,height));
    xml_t xtit = tit;
    return element("title", attrs, xtit);

}
xml_t svggpp::title(const std::string& tit, xml_t attrs)
{
    xml_t xtit = tit;
    return element("title", attrs, xtit);
}



xml_t svggpp::view(const std::string& id, xml_t attr, xml_t body)
{
    if (attr.empty()) {
        attr = xml_t::object();
    }
    attr["id"] = id;
    return element("view", attr, body);
}        

xml_t svggpp::anchor(const std::string& url, xml_t attr, xml_t body)
{
    if (attr.empty()) {
        attr = xml_t::object();
    }
    attr["href"] = url;
    return element("a", attr, body);
}

xml_t svggpp::group(xml_t attrs, xml_t body)
{
    return element("g", attrs, body);
}

xml_t svggpp::svg(xml_t attrs, xml_t body)
{
    return element("svg", attrs, body);
}

xml_t svggpp::svg(double x, double y, double width, double height, 
                   xml_t attrs, xml_t body)
{
    attrs.update(svggpp::point(x,y));
    attrs.update(svggpp::size(width,height));
    return element("svg", attrs, body);
}
xml_t svggpp::use(const std::string& url)
{
    return element("use", { {"href",url} });
}

xml_t svggpp::defs(xml_t body)
{
    return element("defs", {}, body);
}

//
// I/O
//

std::string svggpp::dumps_level(const xml_t& elem, int indent, int level)
{
    if (elem.is_null()) return "";

    std::stringstream ss;
    
    const std::string tab(indent*level, ' ');
    const std::string tabp1(indent*(level+1), ' ');

    const size_t ninline = 1;

    if (elem.is_array()) {
        for (const auto& sub : elem) {
            ss << dumps_level(sub, indent, level); // don't indent one more 
        }
        return ss.str();
    }
    
    if (elem.is_object()) {
        if (! elem.contains(tag_key)) {
            throw std::invalid_argument("no key \"" + tag_key + "\" in " + elem.dump() );
        }
        std::string tag = elem[tag_key];
        ss << tab << "<" << tag;
        const auto& a = attrs(elem);
        if (a.empty()) {
        }
        else {
            if (a.size() <= ninline) {
                for (const auto& kv : a.items()) {
                    ss << " " << kv.key() << "=" << kv.value();
                }
            }
            else {
                for (const auto& kv : a.items()) {
                    ss << "\n" << tabp1 << kv.key() << "=" << kv.value();
                }
            }
        }
        const auto& b = elem[body_key];
        if (b.empty()) {
            ss << "/>\n";
        }
        else {
            ss << ">\n";
            ss << dumps_level(b, indent, level+1);
            ss << tab << "</" << tag << ">\n";
            return ss.str();
        }
        return ss.str();
    }

    if (elem.is_string()) {
        std::string s = elem;
        ss << s << "\n"; // flush left
        return ss.str();
    }

    // number, boolean
    ss << tab << elem.dump() << "\n";
    return ss.str();
}

std::string svggpp::dumps(const xml_t& elem, int indent)
{
    return dumps_level(elem, indent);
}
