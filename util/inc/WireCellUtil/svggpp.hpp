/*
  A simple API for generating SVG.

  It makes use of JSON for Modern C++ (nlohmann/json) as an
  intermediate representation of XML.

*/

#ifndef SVGPPP_HPP
#define SVGPPP_HPP

// Fixme: we ride on custard's coat tails to provide json.hpp,
// eventually this will become a first class citizen in WCT
#include "WireCellUtil/custard/nlohmann/json.hpp"

namespace svggpp {

    // Represent all XML types.  An SVG element type is "elem" is of
    // XML type object and which has exactly three attributes: "tag"
    // of XML type string, "attr" of XML type object and "body" of XML
    // type array of "elem" or of XML type string.  As xml_t is
    // dynamic, "elem", "attr", "body" variable names will be used to
    // indicate when those specific forms are expected.
    using xml_t = nlohmann::json;

    // Provide attribute fragments for top, file-level <svg>
    const xml_t svg_header = {{"version","1.1"}, {"xmlns","http://www.w3.org/2000/svg"}};

    // Return attribute set of element.
    xml_t& attrs(xml_t& elem);
    const xml_t& attrs(const xml_t& elem);

    // Return the body array of elements
    xml_t& body(xml_t& elem);
    const xml_t& body(const xml_t& elem);

    // Return value for key in the attribute set of an elem or
    // directly from the attr.
    std::string get(const xml_t& elem_or_attr, const std::string& key, const std::string& def);
    // As above but throws invalid_argument if key not defined
    std::string get(const xml_t& elem_or_attr, const std::string& key);

    // Return true if element or attribute set has key
    bool has(const xml_t& elem_or_attr, const std::string& key);

    // Return the subset of a with key.
    xml_t subset(const xml_t& elem_or_attr, const std::string& key);

    // return the subset of a with keys.
    xml_t subset(const xml_t& elem_or_attr, const std::vector<std::string>& keys);

    // SVG scalar values are always strings, even when they are
    // numbers.  This wraps up the double cast to get them as a number
    // type.
    template<typename T>
    T cast(const xml_t& scalar)
    {
        return xml_t::parse(scalar.get<std::string>()).get<T>();
    }

    // Set key's value.
    void set(xml_t& elem_or_attr, const std::string& key, const std::string& value);

    // Return the ID or empty string
    std::string id(const xml_t& elem_or_attr);
    // Set id on element (or attr)
    void id(xml_t& elem_or_attr, const std::string& name);

    std::string link(const xml_t& elem_or_attr);

    // represent a point in 2-space
    using point_t = std::pair<double,double>;

    // Function to return attribute object fragments to foist SVG
    // schema into C++.

    // Attribute generators return an attributes object
    xml_t point(double x, double y);
    xml_t point(const point_t& pt);
    xml_t points(const std::vector<point_t>& pts);
    xml_t center(double cx, double cy);
    xml_t size(double width, double height);
    xml_t box(double x, double y, double width, double height);
    xml_t endpoints(double x1, double y1, double x2, double y2);

    // A "view" is any min-x/min-y point with a size
    struct view_t { double x{0}, y{0}, width{0}, height{0}; };
        
    // A viewBox is an <svg> attr defining coordinate system in which
    // children are drawn.  A viewBox may also be used in a <view> to
    // override the viewBox of the <svg> parent of the <view>.

    // Retrieve a viewBox as a view
    view_t viewbox(const xml_t& elem_or_attr);

    // Set a viewBox with a view
    void viewbox(xml_t& elem_or_attr, const view_t& view);

    // A view port is the x,y,width,height attributes of an <svg> and
    // determine where the bounds of the <svg> are located in the
    // parent of the <svg>

    // Retrieve a view port as a view
    view_t viewport(const xml_t& attr);

    // Set a view port with a view
    void viewport(xml_t& attr, const view_t& view);



    /// Produce a CSS style attribute.  Note, there is also a <style>
    /// element.  For now, CSS is treated only as an opaque string.
    xml_t css_style(const std::string& css);
    /// Return a CSS "class" attribute
    xml_t css_class(const std::string& cname);
    xml_t css_class(const std::vector<std::string>& cnames);


    /// Create and return a new element.  User code should only call
    /// this function to construct new elements.  The attrs and body
    /// may evalute to false.  The attrs may be an object or array of
    /// object.  The value of the attributes of these objects must be
    /// of type string.  The body may be of type xml_t object or
    /// string or number or array of these.
    svggpp::xml_t element(const std::string& tag,
                          const svggpp::xml_t& attrs = xml_t::object(),
                          const svggpp::xml_t& body = xml_t::array());

    /// Append child to parent body, return parent.  Parent body will
    /// be an array of elements.  Child may be scalar or array.  See
    /// comments on xml_t above.
    xml_t& append(xml_t& parent, const xml_t& child);

    /// Update attributes of element, returning elem.
    // xml_t& update(xml_t& elem, xml_t attrs);
    // use attr(elem).update(attrs)

    /// Following functions foist specific SVG tags to a C++ function
    /// API.

    /// Shape things:
        
    xml_t text(const std::string& t, double x=0, double y=0,
               xml_t attrs = xml_t::object());
    xml_t line(double x1=0, double y1=0, double x2=1, double y2=1,
               xml_t attrs = { {"stroke","black"} });
    xml_t circle(double cx, double cy, double r=1,
                 xml_t attrs = xml_t::object());
    xml_t circle(point_t cen, double r=1,
                 xml_t attrs = xml_t::object());
    xml_t rect(double x=0, double y=0, double width=1, double height=1,
               xml_t attrs = xml_t::object());
    xml_t polygon(const std::vector<point_t>& pts,
                  xml_t attrs = xml_t::object());

    xml_t title(const std::string& tit,
                double x, double y, double width, double height,
                xml_t attrs = xml_t::object());
    xml_t title(const std::string& tit,
                xml_t attrs = xml_t::object());

    /// Abstract things and grouping

    /// Produces a <style> element.  Note there is also a style
    /// attribute, see css_style().  For now, CSS is treated only as
    /// an opaque string.
    xml_t style(const std::string& css);
    xml_t view(const std::string& id, xml_t attrs = {}, xml_t body = {});
    xml_t anchor(const std::string& url, xml_t attrs = {}, xml_t body = {});
    /// A group <g> element.
    xml_t group(xml_t attrs = {}, xml_t body = {});
    /// An <svg> element
    xml_t svg(xml_t attrs = svg_header, xml_t body=nullptr);
    xml_t svg(double x, double y, double width, double height, 
              xml_t attrs = svg_header, xml_t body=nullptr);
    /// Use some element by referencing it's URL, typically an
    /// "#name" ID.  See Attr::link()
    xml_t use(const std::string& url);
    /// A defs is a spot to park elements for later reference via <use>
    xml_t defs(xml_t body = {});


    
    // I/O
    
    std::string dumps_level(const xml_t& elem, int indent=4, int level=0);

    std::string dumps(const xml_t& elem, int indent=4);

}



#endif
