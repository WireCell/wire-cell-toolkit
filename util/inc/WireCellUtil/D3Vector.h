/** Implement 3-vector arithmetic with std::vector store.
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

        T m_v[3];

      public:

        using value_type = T;   // mimic std::vector
        using coordinate_t = T;

        /// Construct from elements.
        D3Vector(const T& a = 0, const T& b = 0, const T& c = 0)
        {
            this->set(a, b, c);
        }

        // Copy constructor.
        D3Vector(const D3Vector& o)
        {
            if (o.size() == 3) {
                this->set(o.x(), o.y(), o.z());
            }
            else {
                this->invalidate();
            }
        }

        // arrays do not have move constructors or move assignment operators.
        // Move constructor.
        // D3Vector(D3Vector&& o)
        //     : m_v{std::move(o.m_v)}
        // { }

        D3Vector(const T d[3])
        {
            this->set(d[0], d[1], d[2]);
        }

        // Assignment.
        D3Vector& operator=(const D3Vector& o)
        {
            if (o.size()) {
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
            m_v[0] = a;
            m_v[1] = b;
            m_v[2] = c;
        }
        T x(const T& val) { return m_v[0] = val; }
        T y(const T& val) { return m_v[1] = val; }
        T z(const T& val) { return m_v[2] = val; }

        // make this look like std::vector
        const T& at(size_t index) const {
            return m_v[index];
        }
        T& at(size_t index) {
            return m_v[index];
        }
        const T* data() const { return m_v; }
        T* data() { return m_v; }
        const size_t size() const { return 3; }
        void clear() { m_v[0] = m_v[1] = m_v[2] = 0; }
        void resize(size_t /*s*/) { /* no-op */ }

        /// Convert from other typed vector.
        template <class TT>
        D3Vector(const D3Vector<TT>& o)
        {
            this->set(o.x(), o.y(), o.z());
        }

        /// Access elements by name.
        T x() const { return m_v[0]; }
        T y() const { return m_v[1]; }
        T z() const { return m_v[2]; }

        /// Access elements by copy.
        T operator[](std::size_t index) const { return m_v[index]; }

        /// Access elements by reference.
        T& operator[](std::size_t index)
        {
            return m_v[index];  // throw if out of bounds
        }

        /// Return the dot product of this vector and the other.
        T dot(const D3Vector& rhs) const
        {
            T scalar = x() * rhs.x() + y() * rhs.y() + z() * rhs.z();
            return scalar;
        }

        /// Return angle between this vector and the other.
        T angle(const D3Vector& rhs) const
        {
            T m1 = this->magnitude();
            T m2 = rhs.magnitude();
            if (m1 <= 0 || m2 <= 0) {
                return 0;
            }
            T cosine = this->dot(rhs) / (m1 * m2);
            return std::acos(std::min(std::max(cosine, T(-1)), T(1)));
        }

        /// Return the magnitude of this vector.
        T magnitude() const { return std::sqrt(magnitude2()); }
        /// Return the magnitude-squared of this vector.
        T magnitude2() const { return x() * x() + y() * y() + z() * z(); }

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

        template <typename N>
        D3Vector& operator*=(const N& a)
        {
            this->set(x()*a, y()*a, z()*a);
            return *this;
        }

        template <typename N>
        D3Vector& operator/=(const N& a)
        {
            this->set(x()/a, y()/a, z()/a);
            return *this;
        }

        /// defining these opens a fairly nightmarish door.
        /// https://www.artima.com/articles/the-safe-bool-idiom
        // bool operator!() const { return m_v.size() != 3; }
        // operator bool() const { return m_v.size() == 3; }

        // can call set(x,y,z) to revalidate.
        void invalidate()
        {
            /// TODO: no op?
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
