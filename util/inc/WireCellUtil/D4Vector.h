/** Implement 4-vector arithmetic for energy-momentum vectors with std::vector store.
 *
 * This implements a 4-vector primarily for representing particle 4-momentum (E, px, py, pz).
 * See also WireCell::D3Vector.
 */

#ifndef WIRECELLUTIL_D4VECTOR
#define WIRECELLUTIL_D4VECTOR

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>  // for ostream

namespace WireCell {

    /** Dimension-4 vector class for energy-momentum representation.
     *
     * Components: (E, px, py, pz) where E is energy and px,py,pz are momentum components.
     * Adapted from D3Vector.
     *
     */
    template <class T>
    class D4Vector {
        template <class U>
        friend std::ostream& operator<<(std::ostream&, const D4Vector<U>&);

        T m_v[4];

      public:

        using value_type = T;   // mimic std::vector
        using coordinate_t = T;

        /// Construct from elements (E, px, py, pz).
        D4Vector(const T& e = 0, const T& px = 0, const T& py = 0, const T& pz = 0)
        {
            this->set(e, px, py, pz);
        }

        // Copy constructor.
        D4Vector(const D4Vector& o)
        {
            if (o.size() == 4) {
                this->set(o.e(), o.px(), o.py(), o.pz());
            }
            else {
                this->invalidate();
            }
        }

        D4Vector(const T d[4])
        {
            this->set(d[0], d[1], d[2], d[3]);
        }

        // Assignment.
        D4Vector& operator=(const D4Vector& o)
        {
            if (o.size()) {
                this->set(o.e(), o.px(), o.py(), o.pz());
            }
            else {
                this->invalidate();
            }
            return *this;
        }

        /// Set vector from elements (E, px, py, pz);
        void set(const T& e = 0, const T& px = 0, const T& py = 0, const T& pz = 0)
        {
            m_v[0] = e;
            m_v[1] = px;
            m_v[2] = py;
            m_v[3] = pz;
        }
        T e(const T& val) { return m_v[0] = val; }
        T px(const T& val) { return m_v[1] = val; }
        T py(const T& val) { return m_v[2] = val; }
        T pz(const T& val) { return m_v[3] = val; }

        // make this look like std::vector
        const T& at(size_t index) const {
            return m_v[index];
        }
        T& at(size_t index) {
            return m_v[index];
        }
        const T* data() const { return m_v; }
        T* data() { return m_v; }
        const size_t size() const { return 4; }
        void clear() { m_v[0] = m_v[1] = m_v[2] = m_v[3] = 0; }
        void resize(size_t /*s*/) { /* no-op */ }

        /// Convert from other typed vector.
        template <class TT>
        D4Vector(const D4Vector<TT>& o)
        {
            this->set(o.e(), o.px(), o.py(), o.pz());
        }

        /// Access elements by name.
        T e() const { return m_v[0]; }   // energy
        T px() const { return m_v[1]; }  // momentum x
        T py() const { return m_v[2]; }  // momentum y
        T pz() const { return m_v[3]; }  // momentum z

        /// Access elements by copy.
        T operator[](std::size_t index) const { return m_v[index]; }

        /// Access elements by reference.
        T& operator[](std::size_t index)
        {
            return m_v[index];  // throw if out of bounds
        }

        /// Return the 4-vector dot product (Minkowski inner product: E1*E2 - p1·p2).
        T dot(const D4Vector& rhs) const
        {
            T scalar = e() * rhs.e() - px() * rhs.px() - py() * rhs.py() - pz() * rhs.pz();
            return scalar;
        }

        /// Return the invariant mass squared (m² = E² - |p|²).
        T mass2() const 
        { 
            return e() * e() - px() * px() - py() * py() - pz() * pz(); 
        }
        
        /// Return the invariant mass (√(E² - |p|²)).
        T mass() const 
        { 
            T m2 = mass2();
            return m2 >= 0 ? std::sqrt(m2) : 0;
        }

        /// Return the 3-momentum magnitude.
        T p() const { return std::sqrt(px() * px() + py() * py() + pz() * pz()); }
        /// Return the 3-momentum magnitude squared.
        T p2() const { return px() * px() + py() * py() + pz() * pz(); }

        /// Return the transverse momentum.
        T pt() const { return std::sqrt(px() * px() + py() * py()); }
        /// Return the transverse momentum squared.
        T pt2() const { return px() * px() + py() * py(); }

        /// Return rapidity (0.5 * ln((E + pz)/(E - pz))).
        T rapidity() const 
        {
            T eplus = e() + pz();
            T eminus = e() - pz();
            if (eplus <= 0 || eminus <= 0) return 0;
            return 0.5 * std::log(eplus / eminus);
        }

        /// Return pseudorapidity (-ln(tan(θ/2)) where θ is polar angle).
        T eta() const 
        {
            T p_mag = p();
            if (p_mag == 0) return 0;
            T theta = std::acos(pz() / p_mag);
            if (theta <= 0 || theta >= M_PI) return 0;
            return -std::log(std::tan(theta / 2.0));
        }

        /// Return azimuthal angle φ.
        T phi() const 
        {
            return std::atan2(py(), px());
        }

        /// Return the beta factor (|p|/E).
        T beta() const 
        {
            return e() > 0 ? p() / e() : 0;
        }

        /// Return the gamma factor (E/m).
        T gamma() const 
        {
            T m = mass();
            return m > 0 ? e() / m : 0;
        }

        bool operator<(const D4Vector& rhs) const
        {
            if (pz() < rhs.pz()) return true;
            if (py() < rhs.py()) return true;
            if (px() < rhs.px()) return true;
            if (e() < rhs.e()) return true;
            return false;
        }

        D4Vector& operator+=(const D4Vector& other)
        {
            this->set(e() + other.e(), px() + other.px(), py() + other.py(), pz() + other.pz());
            return *this;
        }

        D4Vector& operator-=(const D4Vector& other)
        {
            this->set(e() - other.e(), px() - other.px(), py() - other.py(), pz() - other.pz());
            return *this;
        }

        template <typename N>
        D4Vector& operator*=(const N& a)
        {
            this->set(e()*a, px()*a, py()*a, pz()*a);
            return *this;
        }

        template <typename N>
        D4Vector& operator/=(const N& a)
        {
            this->set(e()/a, px()/a, py()/a, pz()/a);
            return *this;
        }

        // can call set(e,px,py,pz) to revalidate.
        void invalidate()
        {
            /// TODO: no op?
        }
    };

    template <class T>
    std::ostream& operator<<(std::ostream& os, const D4Vector<T>& vec)
    {
        os << "(" << vec.e() << "; " << vec.px() << ", " << vec.py() << ", " << vec.pz() << ")";
        return os;
    }

    template <class T>
    D4Vector<T> operator-(const D4Vector<T>& a, const D4Vector<T>& b)
    {
        return D4Vector<T>(a.e() - b.e(), a.px() - b.px(), a.py() - b.py(), a.pz() - b.pz());
    }

    template <class T>
    D4Vector<T> operator+(const D4Vector<T>& a, const D4Vector<T>& b)
    {
        return D4Vector<T>(a.e() + b.e(), a.px() + b.px(), a.py() + b.py(), a.pz() + b.pz());
    }

    template <class T, typename N>
    D4Vector<T> operator*(const D4Vector<T>& a, const N& s)
    {
        return D4Vector<T>(a.e() * s, a.px() * s, a.py() * s, a.pz() * s);
    }

    template <class T, typename N>
    D4Vector<T> operator*(const N& s, const D4Vector<T>& a)
    {
        return D4Vector<T>(a.e() * s, a.px() * s, a.py() * s, a.pz() * s);
    }

    template <class T, typename N>
    D4Vector<T> operator/(const D4Vector<T>& a, const N& s)
    {
        return D4Vector<T>(a.e() / s, a.px() / s, a.py() / s, a.pz() / s);
    }

    template <class T>
    bool operator==(const D4Vector<T>& a, const D4Vector<T>& b)
    {
        return a.e() == b.e() && a.px() == b.px() && a.py() == b.py() && a.pz() == b.pz();
    }

    template <class T>
    bool operator!=(const D4Vector<T>& a, const D4Vector<T>& b)
    {
        return !(a == b);
    }

}  // namespace WireCell

#endif
