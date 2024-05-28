#ifndef WIRECELLUTIL_DISJOINTRANGE
#define WIRECELLUTIL_DISJOINTRANGE

#include "WireCellUtil/Logging.h" // debugging
#include "WireCellUtil/Type.h"

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/range.hpp>
#include <map>
#include <vector>
#include <stdexcept>

namespace WireCell {

    /** A collection of iterator ranges presented as a flat collection of iterators.

        The disjoint_range defines the following iterator types:

        - minor iterator :: iterates within each iterator range.

        - major iterator :: iterates over the iterator ranges themselves.

        - iterator :: aka disjoint cursor.  It yields a pair holding
          major and minor iterators.
     */

    template<typename DisjointRange, typename ElementType>
    class disjoint_cursor;

    using disjoint_index = std::pair<size_t, size_t>;

    template<typename UserMinorRangeType>
    class disjoint_range {
      public:
        using self_type = disjoint_range<UserMinorRangeType>;

        // using minor_range = UserMinorRangeType;
        using minor_range = boost::sub_range<UserMinorRangeType>;
        using value_type = typename minor_range::value_type;

        // We store minor ranges in vector form to allow fast major/minor indexing.
        // using accum_minor = std::pair<size_t, minor_range>;
        // using major_vec = std::vector<accum_minor>;
        using major_map = std::map<size_t, minor_range>;
        using major_vec = std::vector<size_t>; // random-access map-keys.

        // The flattened iterator
        using iterator = disjoint_cursor<self_type, value_type>;
        using const_iterator = disjoint_cursor<self_type const, value_type const>;

        disjoint_range()
            : last(this, {0,0})
        {
            SPDLOG_TRACE("disjoint_range: ctor empty");
        }
        ~disjoint_range()
        {
            SPDLOG_TRACE("disjoint_range: dtor size {}", size());
        }

        // construct with container-of-ranges type.
        template<typename Container>
        explicit disjoint_range(Container& con)
            : last(this, {0,0})
        {
            fill(con);
            SPDLOG_TRACE("disjoint_range: container ctor size {}", size());
        }

        // construct with a single range
        template<typename RangeIter>
        disjoint_range(RangeIter beg, RangeIter end)
            : last(this, {0,0})
        {
            append(beg, end);
            SPDLOG_TRACE("disjoint_range: iterator ctor size {}", size());
        }


        // Return the flat index coresponding to major/minor indices.
        size_t flat_index(const disjoint_index& djind) const
        {
            const size_t acc = ind2acc[djind.first];
            return acc + djind.second;
            // return ranges[djind.first].first + djind.second;
        }

        iterator begin()
        {
            return iterator(this, {0,0});
        }

        iterator end()
        {
            return iterator(this, {ind2acc.size(), 0});
        }

        const_iterator begin() const
        {
            return const_iterator(this, {0,0});
        }

        const_iterator end() const
        {
            return const_iterator(this, {ind2acc.size(),0});
        }

        // Append each range in container of ranges.
        template<typename Container>
        void fill(Container& con) {
            for (auto& r : con) {
                append(r);
            }
        }

        // Append a range-like
        template<typename Range>
        void append(Range& r) {
            if (r.empty()) return;
            const size_t olds = size();
            acc2mai.emplace(olds, minor_range(r));
            ind2acc.push_back(olds);
            // ranges.emplace_back(size(), minor_range(r));
        }

        // Append an iteration over a range
        template<typename RangeIter>
        void append(RangeIter ribeg, RangeIter riend) {
            if (ribeg == riend) {
                return;
            }
            const size_t olds = size();
            acc2mai.emplace(olds, minor_range(ribeg, riend));
            ind2acc.push_back(olds);
            // ranges.emplace_back(size(), minor_range(ribeg, riend));
        }

        bool empty() const {
            return ind2acc.empty();
        }
        size_t size() const {
            if (ind2acc.empty()) return 0;
            size_t ret = ind2acc.back();
            ret += mai(ret).size();
            return ret;
        }
        size_t nminors() const {
            return ind2acc.size();
        }

        void clear() {
            //ranges.clear();
            acc2mai.clear();
            ind2acc.clear();
        }

        const value_type& operator[](size_t ind) const {
            const size_t rel = ind - last.flat_index();
            last = last + rel;
            return *last;
        }

        value_type& operator[](size_t ind) {
            const size_t flat = last.flat_index();
            if (ind == flat) {
                value_type& ref = *last;
                return ref;
            }
            last += ind - flat;
            return *last;
        }

        const value_type& at(size_t ind) const {
            if (ind < size()) {
                return (*this)[ind];
            }
            throw std::out_of_range("index out of range");
        }
        value_type& at(size_t ind) {
            if (ind < size()) {
                return (*this)[ind];
            }
            throw std::out_of_range("index out of range");
        }

        
        const value_type& operator[](disjoint_index const& djind) const {
            return mai(djind)[djind.second];
        }
        value_type& operator[](disjoint_index const& djind) {
            return mai(djind)[djind.second];
        }
        const value_type& at(disjoint_index const& djind) const {
            const auto& mr = mai(djind);
            if (djind.second < mr.size()) {
                return mr[djind.second];
            }
            throw std::out_of_range("index out of range");
        }
        value_type& at(disjoint_index const& djind) {
            auto& mr = mai(djind);
            if (djind.second < mr.size()) {
                return mr[djind.second];
            }
            throw std::out_of_range("index out of range");
        }

      private:
            
        friend iterator;
        friend const_iterator;

        // fixme: once this all settles down, probably best to try to change
        // major_vec to hold major_map::iterator to avoid the cost of many map
        // lookups.
        //
        // Note, we put the minor range in the map to assure stable iterators
        // over the history of filling.
        major_map acc2mai{};
        major_vec ind2acc{};

        minor_range& mai(size_t acc) {
            auto it = acc2mai.find(acc);
            if (it == acc2mai.end()) {
                throw std::out_of_range("accumulant out of range");
            }
            return it->second;        
        }
        const minor_range& mai(size_t acc) const {
            auto it = acc2mai.find(acc);
            if (it == acc2mai.end()) {
                throw std::out_of_range("accumulant out of range");
            }
            return it->second;        
        }
        minor_range& mai(disjoint_index const& djind) {
            return mai(ind2acc.at(djind.first));
        }
        const minor_range& mai(disjoint_index const& djind) const {
            return mai(ind2acc.at(djind.first));
        }

        // Cache last random access by flat index to reduce number of major
        // jumps to next random access by flat index.  
        mutable iterator last;  
    };

        
    template <typename DisjointRange, typename ElementValue>
    class disjoint_cursor
        : public boost::iterator_facade<disjoint_cursor<DisjointRange, ElementValue>
                                        // value
                                        , ElementValue
                                        // cagegory
                                        , boost::random_access_traversal_tag
                                        // reference
                                        , ElementValue&
                                        >
    {
      public:

        // iterator facade types
        using base_type =
            boost::iterator_facade<disjoint_cursor<DisjointRange, ElementValue>
                                   , ElementValue
                                   , boost::random_access_traversal_tag
                                   , ElementValue&
                                   >;
        using difference_type = typename base_type::difference_type;
        using value_type = typename base_type::value_type;
        using pointer = typename base_type::pointer;
        using reference = typename base_type::reference;

      private:
        DisjointRange* dr;
        // Major/minor index pair
        disjoint_index djind{0,0};

      public:

        // An "end" iterator has major index equal to the major size of the
        // disjoing range.  Otherwise the maj/mai indices should be valid.
        disjoint_cursor(DisjointRange* dr, const disjoint_index& djind)
            : dr(dr), djind(djind) {}

        // copy ctor
        template<typename OtherDisjointCursor>
        disjoint_cursor(OtherDisjointCursor& o)
            : dr(o.dr), djind(o.djind) { }

        // assignment
        template<typename OtherDisjointCursor>
        disjoint_cursor& operator=(OtherDisjointCursor& o) {
            dr = o.dr; djind = o.djind;
        }

        disjoint_index index() const { return djind; }

        // Number of points before the current minor range.
        size_t accumulated() const {
            //return dr->ranges[djind.first].first;
            return dr->ind2acc[djind.first];
        }

        // The number of points in the current minor range.
        size_t minor_size() const {
            //return dr->ranges[djind.first].second.size();
            return dr->mai(djind).size();
        }

        // The flat index of the element at which this iterator resides.
        size_t flat_index() const
        {
            if (dr->empty()) return 0;
            if (at_end()) {
                // at end, ignore minor index
                return dr->size();
            }
            return dr->flat_index(djind);
        }

      private:
        friend class boost::iterator_core_access;

        bool equal(disjoint_cursor const& o) const {
            return dr == o.dr && flat_index() == o.flat_index();
        }

        void increment() {
            if (at_end()) {
                throw std::out_of_range("increment beyond end");
            }
            ++djind.second;
            if (djind.second == minor_size()) {
                ++djind.first;
                djind.second = 0;
            }
        }
        void decrement() {
            if (djind.first == 0 && djind.second == 0) {
                throw std::out_of_range("decrement beyond begin");
            }

            if (at_end()) {
                djind.second=0;                     // set up for next line
            }

            // if sitting one into "the next" min range, back up to previous.
            if (djind.second == 0) {
                --djind.first;
                djind.second = minor_size() - 1;
                return;
            }

            // we are in the middle of a minor range
            --djind.second;
        }

        difference_type distance_to(disjoint_cursor const& other) const
        {
            return other.flat_index() - this->flat_index();
        }
        
        reference dereference() const
        {
            return dr->at(djind);
        }


        bool at_end() const { return djind.first == dr->nminors(); }
        bool at_begin() const { return djind.first==0 && djind.second == 0; }

        void advance(difference_type n)
        {
            if (n && dr->empty()) {
                throw std::out_of_range("advance on empty range");
            }

            while (true) {      // state machine on n


                if (n < 0) { // decrementing.

                    if (at_begin()) {
                        throw std::out_of_range("decrement before begin");
                    }
                    
                    // We are at end or at the start of a minor range.  Jump to
                    // start of previous range.
                    if ( at_end() || djind.second == 0) {
                        --djind.first;
                        djind.second = 0;
                        n += minor_size(); // n becomes less negative
                        continue;
                    }

                    // We are inside a minor range.

                    // Our jump keeps us in our minor range.  Make last jump.
                    if (-n <= (int)djind.second) {
                        djind.second += n;
                        n = 0;
                        continue;
                    }

                    // Need to jump before the start of current minor range.
                    --djind.first;
                    n += djind.second + minor_size(); // n becomes less negative
                    djind.second = 0;
                    continue;
                }

                if (n > 0) {    // incrementing

                    if (at_end()) {
                        throw std::out_of_range("advace beyond end");
                    }

                    const int n_left = minor_size() - djind.second;

                    // Our jump keeps us in the minor range.  Make last jump.
                    if (n < n_left) {
                        djind.second += n;
                        n = 0;
                        continue;
                    }

                    // Need to jump beyond current minor range
                    ++djind.first;
                    djind.second = 0;
                    n -= n_left;
                    continue;
                }

                return;         // n == 0
            }
        }
    };
}

#endif
