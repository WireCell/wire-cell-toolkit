/** Implement 2-vector arithmetic with std::vector store.
 *
 */

#ifndef WIRECELLDATA_VECTOR2D
#define WIRECELLDATA_VECTOR2D

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>  // for ostream

namespace WireCell {

    /** Dimension-2 vector class.
     * modified from D2Vector.h
     *
     */
    template <class T>
    class D2Vector {
        template <class U>
        friend std::ostream& operator<<(std::ostream&, const D2Vector<U>&);

        T m_v[2];

      public:

        using value_type = T;   // mimic std::vector
        using coordinate_t = T;

        /// Construct from elements.
        D2Vector(const T& a = 0, const T& b = 0)
        {
            this->set(a, b);
        }

        // Copy constructor.
        D2Vector(const D2Vector& o)
        {
            if (o.size() == 2) {
                this->set(o.x(), o.y());
            }
            else {
                this->invalidate();
            }
        }

        // arrays do not have move constructors or move assignment operators.
        // Move constructor.
        // D2Vector(D2Vector&& o)
        //     : m_v{std::move(o.m_v)}
        // { }

        D2Vector(const T d[2])
        {
            this->set(d[0], d[1]);
        }

        // Assignment.
        D2Vector& operator=(const D2Vector& o)
        {
            if (o.size()) {
                this->set(o.x(), o.y());
            }
            else {
                this->invalidate();
            }
            return *this;
        }

        /// Set vector from elements;
        void set(const T& a = 0, const T& b = 0)
        {
            m_v[0] = a;
            m_v[1] = b;
        }
        T x(const T& val) { return m_v[0] = val; }
        T y(const T& val) { return m_v[1] = val; }

        // make this look like std::vector
        const T& at(size_t index) const {
            return m_v[index];
        }
        T& at(size_t index) {
            return m_v[index];
        }
        const T* data() const { return m_v; }
        T* data() { return m_v; }
        const size_t size() const { return 2; }
        void clear() { m_v[0] = m_v[1] = 0; }
        void resize(size_t /*s*/) { /* no-op */ }

        /// Convert from other typed vector.
        template <class TT>
        D2Vector(const D2Vector<TT>& o)
        {
            this->set(o.x(), o.y());
        }

        /// Access elements by name.
        T x() const { return m_v[0]; }
        T y() const { return m_v[1]; }

        /// Access elements by copy.
        T operator[](std::size_t index) const { return m_v[index]; }

        /// Access elements by reference.
        T& operator[](std::size_t index)
        {
            return m_v[index];  // throw if out of bounds
        }

        /// Return the dot product of this vector and the other.
        T dot(const D2Vector& rhs) const
        {
            T scalar = x() * rhs.x() + y() * rhs.y();
            return scalar;
        }

        /// Return angle between this vector and the other.
        T angle(const D2Vector& rhs) const
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
        T magnitude() const { return std::sqrt(x() * x() + y() * y()); }

        /// Return a normalized vector in the direction of this vector.
        D2Vector norm() const
        {
            T m = this->magnitude();
            if (m <= 0) {
                return D2Vector();
            }
            return D2Vector(x() / m, y() / m);
        }

        bool operator<(const D2Vector& rhs) const
        {
            if (y() < rhs.y()) return true;
            if (x() < rhs.x()) return true;
            return false;
        }

        D2Vector& operator+=(const D2Vector& other)
        {
            this->set(x() + other.x(), y() + other.y());
            return *this;
        }

        D2Vector& operator-=(const D2Vector& other)
        {
            this->set(x() - other.x(), y() - other.y());
            return *this;
        }

        /// defining these opens a fairly nightmarish door.
        /// https://www.artima.com/articles/the-safe-bool-idiom
        // bool operator!() const { return m_v.size() != 2; }
        // operator bool() const { return m_v.size() == 2; }

        // can call set(x,y) to revalidate.
        void invalidate()
        {
            /// TODO: no op?
        }
    };

    template <class T>
    std::ostream& operator<<(std::ostream& os, const D2Vector<T>& vec)
    {
        os << "(" << vec.x() << " " << vec.y() << ")";
        return os;
    }

    template <class T>
    D2Vector<T> operator-(const D2Vector<T>& a, const D2Vector<T>& b)
    {
        return D2Vector<T>(a.x() - b.x(), a.y() - b.y());
    }

    template <class T>
    D2Vector<T> operator+(const D2Vector<T>& a, const D2Vector<T>& b)
    {
        return D2Vector<T>(a.x() + b.x(), a.y() + b.y());
    }

    template <class T, typename N>
    D2Vector<T> operator*(const D2Vector<T>& a, const N& s)
    {
        return D2Vector<T>(a.x() * s, a.y() * s);
    }

    template <class T, typename N>
    D2Vector<T> operator*(const N& s, const D2Vector<T>& a)
    {
        return D2Vector<T>(a.x() * s, a.y() * s);
    }

    template <class T, typename N>
    D2Vector<T> operator/(const D2Vector<T>& a, const N& s)
    {
        return D2Vector<T>(a.x() / s, a.y() / s);
    }

    template <class T>
    bool operator==(const D2Vector<T>& a, const D2Vector<T>& b)
    {
        return a.x() == b.x() && a.y() == b.y();
    }

    template <class T>
    bool operator!=(const D2Vector<T>& a, const D2Vector<T>& b)
    {
        return !(a == b);
    }

}  // namespace WireCell

#endif
