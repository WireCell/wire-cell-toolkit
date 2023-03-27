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
    if (empty()) {
        return false;
    }
    for (int ind = 0; ind < 3; ++ind) {
        const double p = point[ind];
        const double b1 = m_bounds.first[ind];
        const double b2 = m_bounds.second[ind];

        if (b1 < b2) {
            if (p < b1 or p > b2) return false;
            continue;
        }
        if (b2 < b1) {
            if (p < b2 or p > b1) return false;
            continue;
        }

        // if equal, then zero width dimension, don't test.
        continue;
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
