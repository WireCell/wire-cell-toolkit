/** Vector a 3-vector of double.
 *
 * See also WireCell::Point.
 */

#ifndef WIRECELLDATA_VECTOR
#define WIRECELLDATA_VECTOR

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>  // for ostream

namespace WireCell {

    /** Dimension-3 vector class.
     *
     * Adapted in laziness from:
     * http://rosettacode.org/wiki/Vector_products#C.2B.2B
     *
     */
    template <class T>
    class D3Vector {
        template <class U>
        friend std::ostream& operator<<(std::ostream&, const D3Vector<U>&);

        typedef std::vector<T> D3VectorStore;
        D3VectorStore m_v{3};

      public:

        using coordinate_t = T;

        /// Construct from elements.
        D3Vector(const T& a = 0, const T& b = 0, const T& c = 0)
            : m_v(3,0)
        {
            this->set(a, b, c);
        }

        // Copy constructor.
        D3Vector(const D3Vector& o)
            : m_v(3,0)
        {
            if (o) {
                this->set(o.x(), o.y(), o.z());
            }
            else {
                this->invalidate();
            }
        }

        // Move constructor.
        D3Vector(D3Vector&& o)
            : m_v(std::move(o.m_v))
        { }

        D3Vector(const T d[3])
            : m_v(3,0)
        {
            this->set(d[0], d[1], d[2]);
        }

        // Assignment.
        D3Vector& operator=(const D3Vector& o)
        {
            if (o) {
                this->set(o.x(), o.y(), o.z());
            }
            else {
                this->invalidate();
            }
            return *this;
        }

        /// Set vector from elements;
        void set(const T& a = 0, const T& b = 0, const T& c = 0)
        {
            m_v.resize(3, 0);
            m_v[0] = a;
            m_v[1] = b;
            m_v[2] = c;
        }
        T x(const T& val) { return m_v[0] = val; }
        T y(const T& val) { return m_v[1] = val; }
        T z(const T& val) { return m_v[2] = val; }

        /// Convert from other typed vector.
        template <class TT>
        D3Vector(const D3Vector<TT>& o)
            : m_v(3,0)
        {
            this->set(o.x(), o.y(), o.z());
        }

        /// Access elements by name.
        T x() const { return m_v[0]; }
        T y() const { return m_v[1]; }
        T z() const { return m_v[2]; }

        /// Access elements by copy.
        T operator[](std::size_t index) const { return m_v.at(index); }

        /// Access elements by reference.
        T& operator[](std::size_t index)
        {
            return m_v.at(index);  // throw if out of bounds
        }

        /// Return the dot product of this vector and the other.
        T dot(const D3Vector& rhs) const
        {
            T scalar = x() * rhs.x() + y() * rhs.y() + z() * rhs.z();
            return scalar;
        }

        /// Return the magnitude of this vector.
        T magnitude() const { return std::sqrt(x() * x() + y() * y() + z() * z()); }

        /// Return a normalized vector in the direction of this vector.
        D3Vector norm() const
        {
            T m = this->magnitude();
            if (m <= 0) {
                return D3Vector();
            }
            return D3Vector(x() / m, y() / m, z() / m);
        }

        /// Return the cross product of this vector and the other.
        D3Vector cross(const D3Vector& rhs) const
        {
            T a = y() * rhs.z() - z() * rhs.y();
            T b = z() * rhs.x() - x() * rhs.z();
            T c = x() * rhs.y() - y() * rhs.x();
            D3Vector product(a, b, c);
            return product;
        }

        /// Return the triple cross product of this vector and the other two.
        D3Vector triplevec(D3Vector& a, D3Vector& b) const { return cross(a.cross(b)); }
        /// Return the dot-cross product of this vector and the other two.
        T triplescal(D3Vector& a, D3Vector& b) const { return dot(a.cross(b)); }

        bool operator<(const D3Vector& rhs) const
        {
            if (z() < rhs.z()) return true;
            if (y() < rhs.y()) return true;
            if (x() < rhs.x()) return true;
            return false;
        }

        D3Vector& operator+=(const D3Vector& other)
        {
            this->set(x() + other.x(), y() + other.y(), z() + other.z());
            return *this;
        }

        D3Vector& operator-=(const D3Vector& other)
        {
            this->set(x() - other.x(), y() - other.y(), z() - other.z());
            return *this;
        }

        bool operator!() const { return m_v.size() != 3; }
        operator bool() const { return m_v.size() == 3; }
        // can call set(x,y,z) to revalidate.
        void invalidate()
        {
            m_v.clear();
            m_v.shrink_to_fit();
        }
    };

    template <class T>
    std::ostream& operator<<(std::ostream& os, const D3Vector<T>& vec)
    {
        os << "(" << vec.x() << " " << vec.y() << " " << vec.z() << ")";
        return os;
    }

    template <class T>
    D3Vector<T> operator-(const D3Vector<T>& a, const D3Vector<T>& b)
    {
        return D3Vector<T>(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
    }

    template <class T>
    D3Vector<T> operator+(const D3Vector<T>& a, const D3Vector<T>& b)
    {
        return D3Vector<T>(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
    }

    template <class T, typename N>
    D3Vector<T> operator*(const D3Vector<T>& a, const N& s)
    {
        return D3Vector<T>(a.x() * s, a.y() * s, a.z() * s);
    }

    template <class T, typename N>
    D3Vector<T> operator*(const N& s, const D3Vector<T>& a)
    {
        return D3Vector<T>(a.x() * s, a.y() * s, a.z() * s);
    }

    template <class T, typename N>
    D3Vector<T> operator/(const D3Vector<T>& a, const N& s)
    {
        return D3Vector<T>(a.x() / s, a.y() / s, a.z() / s);
    }

    template <class T>
    bool operator==(const D3Vector<T>& a, const D3Vector<T>& b)
    {
        return a.x() == b.x() && a.y() == b.y() && a.z() == b.z();
    }

    template <class T>
    bool operator!=(const D3Vector<T>& a, const D3Vector<T>& b)
    {
        return !(a == b);
    }

}  // namespace WireCell

#endif
