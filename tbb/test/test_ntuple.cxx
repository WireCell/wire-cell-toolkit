#include <tuple>

template<typename T0, typename ... TN>
struct MyTuple
{
    typedef std::tuple<std::size_t, T0, TN...> tagged_tuple;
};


template<typename T, unsigned N, typename... REST>
struct generate_tuple_type
{
    typedef typename generate_tuple_type<T, N-1, T, REST...>::type type;
};

template<typename T, typename... REST>
struct generate_tuple_type<T, 0, REST...>
{
    typedef MyTuple<REST...> type;
};


int main()
{
    generate_tuple_type<float, 3>::type::tagged_tuple mt(0, 1.1, 2.2, 3.3);
    return std::get<0>(mt);
}

// also https://en.cppreference.com/w/cpp/utility/integer_sequence
