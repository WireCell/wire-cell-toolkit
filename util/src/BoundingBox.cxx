#include "WireCellUtil/BoundingBox.h"

WireCell::BoundingBox::BoundingBox(const Point& initial)
  : m_bounds(initial, initial)
  , m_initialized(true)
{
}
WireCell::BoundingBox::BoundingBox(const Ray& initial)
  : m_bounds(initial)
  , m_initialized(true)
{
}
void WireCell::BoundingBox::operator()(const Ray& r)
{
    (*this)(r.first);
    (*this)(r.second);
}
void WireCell::BoundingBox::operator()(const Point& p)
{
    if (empty()) {
        m_bounds.first = p;
        m_bounds.second = p;
        m_initialized = true;
        return;
    }

    for (int ind = 0; ind < 3; ++ind) {
        if (p[ind] < m_bounds.first[ind]) {
            m_bounds.first[ind] = p[ind];
        }
        if (p[ind] > m_bounds.second[ind]) {
            m_bounds.second[ind] = p[ind];
        }
    }
}

bool WireCell::BoundingBox::inside(const Point& point) const
{
    // if (empty()) {
    //     return false;
    // }
    // for (int ind = 0; ind < 3; ++ind) {
    //     const double p = point[ind];
    //     const double b1 = m_bounds.first[ind];
    //     const double b2 = m_bounds.second[ind];

    //     if (b1 < b2) {
    //         if (p < b1 or p > b2) return false;
    //         continue;
    //     }
    //     if (b2 < b1) {
    //         if (p < b2 or p > b1) return false;
    //         continue;
    //     }

    //     // if equal, then zero width dimension, don't test.
    //     continue;
    // }
    // return true;

    if (empty()) {
        return false;
    }
    
    for (int ind = 0; ind < 3; ++ind) {
        const double p = point[ind];
        const double min_val = std::min(m_bounds.first[ind], m_bounds.second[ind]);
        const double max_val = std::max(m_bounds.first[ind], m_bounds.second[ind]);
        
        // Skip zero-width dimensions
        if (min_val == max_val) {
            continue;
        }
        
        // Check if point coordinate is outside bounds
        if (p < min_val || p > max_val) {
            return false;
        }
    }
    
    return true;
}

void WireCell::BoundingBox::pad_rel(double relative)
{
    auto vec = ray_vector(m_bounds)*relative;
    m_bounds.first = m_bounds.first - vec;
    m_bounds.second = m_bounds.second + vec;
}
void WireCell::BoundingBox::pad_abs(double distance)
{
    auto vec = ray_unit(m_bounds)*distance;
    m_bounds.first = m_bounds.first - vec;
    m_bounds.second = m_bounds.second + vec;
}


// Add to BoundingBox.cxx
WireCell::BoundingBox WireCell::BoundingBox::unite(const BoundingBox& other) const
{
    if (other.empty()) {
        return *this;
    }
    if (this->empty()) {
        return other;
    }
    
    BoundingBox result;
    result.m_initialized = true;
    
    for (int ind = 0; ind < 3; ++ind) {
        result.m_bounds.first[ind] = std::min(m_bounds.first[ind], other.m_bounds.first[ind]);
        result.m_bounds.second[ind] = std::max(m_bounds.second[ind], other.m_bounds.second[ind]);
    }
    
    return result;
}

double WireCell::BoundingBox::distance(const Point& point) const
{
    if (empty()) {
        return std::numeric_limits<double>::max();
    }
    
    if (inside(point)) {
        return 0.0;
    }
    
    double dist_squared = 0.0;
    
    for (int ind = 0; ind < 3; ++ind) {
        const double p = point[ind];
        const double min_val = m_bounds.first[ind];
        const double max_val = m_bounds.second[ind];
        
        if (p < min_val) {
            dist_squared += (min_val - p) * (min_val - p);
        }
        else if (p > max_val) {
            dist_squared += (p - max_val) * (p - max_val);
        }
    }
    
    return std::sqrt(dist_squared);
}

WireCell::Point WireCell::BoundingBox::center() const
{
    if (empty()) {
        return Point();
    }
    
    return Point(
        (m_bounds.first.x() + m_bounds.second.x()) / 2.0,
        (m_bounds.first.y() + m_bounds.second.y()) / 2.0,
        (m_bounds.first.z() + m_bounds.second.z()) / 2.0
    );
}

WireCell::Vector WireCell::BoundingBox::dimensions() const
{
    if (empty()) {
        return Vector(0, 0, 0);
    }
    
    return Vector(
        std::abs(m_bounds.second.x() - m_bounds.first.x()),
        std::abs(m_bounds.second.y() - m_bounds.first.y()),
        std::abs(m_bounds.second.z() - m_bounds.first.z())
    );
}

double WireCell::BoundingBox::volume() const
{
    if (empty()) {
        return 0.0;
    }
    
    Vector dims = dimensions();
    return dims.x() * dims.y() * dims.z();
}

WireCell::BoundingBox WireCell::BoundingBox::combine(const std::vector<BoundingBox>& boxes)
{
    BoundingBox result;
    
    for (const auto& box : boxes) {
        if (!box.empty()) {
            result = result.empty() ? box : result.unite(box);
        }
    }
    
    return result;
}

// Implementation for intersect method
WireCell::BoundingBox WireCell::BoundingBox::intersect(const BoundingBox& other) const
{
    if (empty() || other.empty()) {
        // If either box is empty, intersection is empty
        return BoundingBox();
    }
    
    // Create a new bounding box for the intersection
    BoundingBox result;
    result.m_initialized = false;
    
    // For each dimension, find the intersection range
    for (int ind = 0; ind < 3; ++ind) {
        const double this_min = std::min(m_bounds.first[ind], m_bounds.second[ind]);
        const double this_max = std::max(m_bounds.first[ind], m_bounds.second[ind]);
        const double other_min = std::min(other.m_bounds.first[ind], other.m_bounds.second[ind]);
        const double other_max = std::max(other.m_bounds.first[ind], other.m_bounds.second[ind]);
        
        // Calculate intersection
        const double intersect_min = std::max(this_min, other_min);
        const double intersect_max = std::min(this_max, other_max);
        
        // If no intersection in any dimension, return empty box
        if (intersect_min > intersect_max) {
            return BoundingBox();
        }
        
        // Set this dimension's intersection bounds
        if (!result.m_initialized) {
            result.m_bounds.first[ind] = intersect_min;
            result.m_bounds.second[ind] = intersect_max;
        }
        else {
            result.m_bounds.first[ind] = intersect_min;
            result.m_bounds.second[ind] = intersect_max;
        }
    }
    
    // If we got here, we have a valid intersection
    result.m_initialized = true;
    return result;
}

// Implementation for overlaps method
bool WireCell::BoundingBox::overlaps(const BoundingBox& other) const
{
    if (empty() || other.empty()) {
        return false;
    }
    
    // Check for overlap in each dimension
    for (int ind = 0; ind < 3; ++ind) {
        const double this_min = std::min(m_bounds.first[ind], m_bounds.second[ind]);
        const double this_max = std::max(m_bounds.first[ind], m_bounds.second[ind]);
        const double other_min = std::min(other.m_bounds.first[ind], other.m_bounds.second[ind]);
        const double other_max = std::max(other.m_bounds.first[ind], other.m_bounds.second[ind]);
        
        // If there's no overlap in any dimension, boxes don't overlap
        if (this_max < other_min || this_min > other_max) {
            return false;
        }
    }
    
    // If we get here, boxes overlap in all dimensions
    return true;
}

// Implementation for closest_point method
WireCell::Point WireCell::BoundingBox::closest_point(const Point& point) const
{
    if (empty()) {
        return Point();
    }
    
    if (inside(point)) {
        return point;  // The point itself is the closest if it's inside
    }
    
    Point closest;
    
    for (int ind = 0; ind < 3; ++ind) {
        const double p = point[ind];
        const double min_val = std::min(m_bounds.first[ind], m_bounds.second[ind]);
        const double max_val = std::max(m_bounds.first[ind], m_bounds.second[ind]);
        
        // Clamp the coordinate to be within the box bounds
        closest[ind] = std::max(min_val, std::min(p, max_val));
    }
    
    return closest;
}

std::vector<double> WireCell::BoundingBox::axis_distances(const Point& point, int axis) const
{
    const int axis1 = (axis+1)%3;
    const int axis2 = (axis+2)%3;
    const auto vmin = m_bounds.first;
    const auto vmax = m_bounds.second;

    std::vector<double> adists;

    if (point[axis1] < vmin[axis1] || point[axis1] > vmax[axis1]) return adists;
    if (point[axis2] < vmin[axis2] || point[axis2] > vmax[axis2]) return adists;

    return {vmin[axis] - point[axis], vmax[axis] - point[axis]};
}
