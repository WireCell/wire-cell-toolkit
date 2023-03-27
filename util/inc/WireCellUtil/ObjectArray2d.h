/** This class provides a 2D array-like structure which holds things
 * of templated type.  It's not fancy.
 */

#ifndef WIRECELL_OBJECTARRAY2D
#define WIRECELL_OBJECTARRAY2D

#include <vector>

namespace WireCell {

    // little helper to give something that looks like a 2D array of objects.
    template <typename Thing>
    class ObjectArray2d {
       public:
        typedef std::vector<Thing> store_t;
        typedef typename store_t::iterator iterator;
        typedef typename store_t::const_iterator const_iterator;

        ObjectArray2d(const ObjectArray2d& other)
        {
            (*this) = other;
        }

        ObjectArray2d& operator=(const ObjectArray2d& other)
        {
            m_nrows = other.m_nrows;
            m_ncols = other.m_ncols;
            m_things.clear();
            m_things.reserve(m_nrows*m_ncols);
            for (const auto& other_thing : other) {
                m_things.emplace_back(other_thing);
            }
            return *this;
        }

        ObjectArray2d(size_t nrows = 0, size_t ncols = 0)
        {
            resize(nrows, ncols);

        }

        void resize(size_t nrows, size_t ncols)
        {
            m_nrows = nrows;
            m_ncols = ncols;
            reset();
        }

        void reset()
        {
            m_things.clear();
            m_things.reserve(m_nrows*m_ncols);
            for (size_t ind=0; ind<m_nrows*m_ncols; ++ind) {
                m_things.emplace_back();
            }
        }

        size_t nrows() const { return m_nrows; }
        size_t ncols() const { return m_ncols; }
        size_t size() const { return m_things.size(); }

        const Thing& operator()(size_t irow, size_t icol) const { return m_things.at(icol + m_ncols * irow); }
        Thing& operator()(size_t irow, size_t icol) { return m_things.at(icol + m_ncols * irow); }

        // range based access to the underlying things
        iterator begin() { return m_things.begin(); }
        iterator end() { return m_things.end(); }
        const_iterator begin() const { return m_things.begin(); }
        const_iterator end() const { return m_things.end(); }

       private:
        // column major
        size_t m_nrows{0}, m_ncols{0};
        store_t m_things;
    };
}  // namespace WireCell

#endif
