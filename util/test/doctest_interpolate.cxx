#include "WireCellUtil/Interpolate.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#pragma GCC diagnostic pop

#include <iomanip>
#include <iostream>

#include "WireCellUtil/doctest.h"

// https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/cardinal_cubic_b.html

using boost::math::interpolators::cardinal_cubic_b_spline;

using namespace WireCell;

template<typename Terp>
void do_terp(Terp& terp)
{
    REQUIRE(terp(0) == 0);
    REQUIRE(terp(9) == 9);
    REQUIRE(terp(-1) == 0);
    REQUIRE(terp(10) == 9);

    std::vector<typename Terp::ytype> terped;
    const size_t num = 100 + 1;
    terp(std::back_inserter(terped), num, 0, 0.1);
    // for (size_t ind=0; ind<num; ++ind) {
    //     std::cerr << "[" << ind << "] = " << terped[ind] << "\n";
    // }
    REQUIRE(terped.size() == num);
    REQUIRE(terped[0] == 0);
    REQUIRE(terped.back() == 9);
    std::cerr << "terped[90]-9=" << terped[90]-9 << "\n";
    REQUIRE(std::abs(terped[90]-9) < 2e-6);
}

template<typename X, typename Y=X>
void test_linterp()
{
    std::vector<Y> ydata{0,1,2,3,4,5,6,7,8,9};
    linterp<X,Y> terp(ydata.begin(), ydata.end(), 0, 1);
    do_terp(terp);
}
template<typename X, typename Y=X>
void test_irrterp()
{
    std::map<X,Y> data{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9}};
    irrterp<X,Y> terp(data.begin(), data.end());
    do_terp(terp);
}

void test_boost()
{
    std::vector<double> f{0.01, -0.02, 0.3, 0.8, 1.9, -8.78, -22.6};
    const double xstep = 0.01;
    const double x0 = xstep;
    cardinal_cubic_b_spline<double> spline(f.begin(), f.end(), x0, xstep);

    linterp<double> lin(f.begin(), f.end(), x0, xstep);

    for (double x = 0; x < x0 + xstep * 10; x += 0.1 * xstep) {
        std::cout << std::setprecision(3) << std::fixed << "x=" << x << "\tlin(x)=" << lin(x)
                  << "\tspline(x)=" << spline(x) << "\n";
    }
}

TEST_CASE("util interpolation")
{
    std::cerr << "testing irrterp:\n";
    test_irrterp<float>();
    test_irrterp<double>();
    test_irrterp<float,  float>();
    test_irrterp<double, float>();
    test_irrterp<float,  double>();
    test_irrterp<double, double>();
    std::cerr << "testing linterp:\n";
    test_linterp<float>();
    test_linterp<double>();
    test_linterp<float,  float>();
    test_linterp<double, float>();
    test_linterp<float,  double>();
    test_linterp<double, double>();
}

TEST_CASE("util regular linear interpolation")
{
   std::cout << "--- Testing linterp<double, double> ---" << std::endl;

    // 1. Default Constructor
    linterp<double, double> interp1;
    REQUIRE(interp1.get_data().empty());
    REQUIRE(interp1.get_le() == 0.0);
    REQUIRE(interp1.get_re() == 0.0);
    REQUIRE(interp1.get_step() == 0.0);
    std::cout << "Default constructor: OK" << std::endl;

    // 2. Prepare an original object
    linterp<double, double> original_interp;
    original_interp.set_data({1.1, 2.2, 3.3});
    original_interp.set_params(0.0, 10.0, 1.0);
    REQUIRE(original_interp.get_data().size() == 3);
    REQUIRE(original_interp.get_le() == 0.0);

    // 3. Copy Constructor
    linterp<double, double> copy_constructed_interp = original_interp;
    REQUIRE(copy_constructed_interp.get_data().size() == 3);
    REQUIRE(copy_constructed_interp.get_data()[0] == 1.1);
    REQUIRE(copy_constructed_interp.get_le() == 0.0);
    REQUIRE(copy_constructed_interp.get_re() == 10.0);
    REQUIRE(copy_constructed_interp.get_step() == 1.0);
    // Ensure original is unchanged
    REQUIRE(original_interp.get_data().size() == 3);
    std::cout << "Copy constructor: OK" << std::endl;

    // 4. Move Constructor
    linterp<double, double> another_original_interp;
    another_original_interp.set_data({4.4, 5.5});
    another_original_interp.set_params(20.0, 30.0, 2.0);

    linterp<double, double> move_constructed_interp = std::move(another_original_interp);
    REQUIRE(move_constructed_interp.get_data().size() == 2);
    REQUIRE(move_constructed_interp.get_data()[0] == 4.4);
    REQUIRE(move_constructed_interp.get_le() == 20.0);
    REQUIRE(move_constructed_interp.get_re() == 30.0);
    REQUIRE(move_constructed_interp.get_step() == 2.0);
    // Ensure original is in a valid, moved-from state (vector empty, primitives default)
    REQUIRE(another_original_interp.get_data().empty());
    REQUIRE(another_original_interp.get_le() == 0.0);
    REQUIRE(another_original_interp.get_re() == 0.0);
    REQUIRE(another_original_interp.get_step() == 0.0);
    std::cout << "Move constructor: OK" << std::endl;

    // 5. Copy Assignment Operator
    linterp<double, double> assigned_copy_interp;
    assigned_copy_interp.set_data({9.9});
    assigned_copy_interp.set_params(100.0, 200.0, 10.0);

    linterp<double, double> source_for_copy_assign;
    source_for_copy_assign.set_data({6.6, 7.7});
    source_for_copy_assign.set_params(40.0, 50.0, 3.0);

    assigned_copy_interp = source_for_copy_assign;
    REQUIRE(assigned_copy_interp.get_data().size() == 2);
    REQUIRE(assigned_copy_interp.get_data()[0] == 6.6);
    REQUIRE(assigned_copy_interp.get_le() == 40.0);
    REQUIRE(assigned_copy_interp.get_re() == 50.0);
    REQUIRE(assigned_copy_interp.get_step() == 3.0);
    // Ensure source is unchanged
    REQUIRE(source_for_copy_assign.get_data().size() == 2);
    std::cout << "Copy assignment operator: OK" << std::endl;

    // 6. Move Assignment Operator
    linterp<double, double> assigned_move_interp;
    assigned_move_interp.set_data({12.12});
    assigned_move_interp.set_params(300.0, 400.0, 20.0);

    linterp<double, double> source_for_move_assign;
    source_for_move_assign.set_data({8.8, 9.9, 10.10});
    source_for_move_assign.set_params(60.0, 70.0, 4.0);

    assigned_move_interp = std::move(source_for_move_assign);
    REQUIRE(assigned_move_interp.get_data().size() == 3);
    REQUIRE(assigned_move_interp.get_data()[0] == 8.8);
    REQUIRE(assigned_move_interp.get_le() == 60.0);
    REQUIRE(assigned_move_interp.get_re() == 70.0);
    REQUIRE(assigned_move_interp.get_step() == 4.0);
    // Ensure source is in a valid, moved-from state
    REQUIRE(source_for_move_assign.get_data().empty());
    REQUIRE(source_for_move_assign.get_le() == 0.0);
    REQUIRE(source_for_move_assign.get_re() == 0.0);
    REQUIRE(source_for_move_assign.get_step() == 0.0);
    std::cout << "Move assignment operator: OK" << std::endl;

    std::cout << "All linterp<double, double> tests passed!" << std::endl << std::endl;
}

TEST_CASE("util irregular linear interpolation")
{
    std::cout << "--- Testing irrterp<double, double> ---" << std::endl;

    // 1. Default Constructor
    irrterp<double, double> interp1;
    REQUIRE(interp1.get_points().empty());
    std::cout << "Default constructor: OK" << std::endl;

    // 2. Prepare an original object
    irrterp<double, double> original_interp;
    original_interp.add_point(1.0, 10.0);
    original_interp.add_point(2.0, 20.0);
    original_interp.add_point(3.0, 30.0);
    REQUIRE(original_interp.get_points().size() == 3);
    REQUIRE(original_interp.get_points().at(1.0) == 10.0);

    // 3. Copy Constructor
    irrterp<double, double> copy_constructed_interp = original_interp;
    REQUIRE(copy_constructed_interp.get_points().size() == 3);
    REQUIRE(copy_constructed_interp.get_points().at(1.0) == 10.0);
    REQUIRE(copy_constructed_interp.get_points().at(2.0) == 20.0);
    REQUIRE(copy_constructed_interp.get_points().at(3.0) == 30.0);
    // Ensure original is unchanged
    REQUIRE(original_interp.get_points().size() == 3);
    std::cout << "Copy constructor: OK" << std::endl;

    // 4. Move Constructor
    irrterp<double, double> another_original_interp;
    another_original_interp.add_point(4.0, 40.0);
    another_original_interp.add_point(5.0, 50.0);

    irrterp<double, double> move_constructed_interp = std::move(another_original_interp);
    REQUIRE(move_constructed_interp.get_points().size() == 2);
    REQUIRE(move_constructed_interp.get_points().at(4.0) == 40.0);
    REQUIRE(move_constructed_interp.get_points().at(5.0) == 50.0);
    // Ensure original is empty after move
    REQUIRE(another_original_interp.get_points().empty());
    std::cout << "Move constructor: OK" << std::endl;

    // 5. Copy Assignment Operator
    irrterp<double, double> assigned_copy_interp;
    assigned_copy_interp.add_point(100.0, 1000.0); // Some initial data

    irrterp<double, double> source_for_copy_assign;
    source_for_copy_assign.add_point(6.0, 60.0);
    source_for_copy_assign.add_point(7.0, 70.0);

    assigned_copy_interp = source_for_copy_assign;
    REQUIRE(assigned_copy_interp.get_points().size() == 2);
    REQUIRE(assigned_copy_interp.get_points().at(6.0) == 60.0);
    REQUIRE(assigned_copy_interp.get_points().at(7.0) == 70.0);
    // Ensure source is unchanged
    REQUIRE(source_for_copy_assign.get_points().size() == 2);
    std::cout << "Copy assignment operator: OK" << std::endl;

    // 6. Move Assignment Operator
    irrterp<double, double> assigned_move_interp;
    assigned_move_interp.add_point(200.0, 2000.0); // Some initial data

    irrterp<double, double> source_for_move_assign;
    source_for_move_assign.add_point(8.0, 80.0);
    source_for_move_assign.add_point(9.0, 90.0);

    assigned_move_interp = std::move(source_for_move_assign);
    REQUIRE(assigned_move_interp.get_points().size() == 2);
    REQUIRE(assigned_move_interp.get_points().at(8.0) == 80.0);
    REQUIRE(assigned_move_interp.get_points().at(9.0) == 90.0);
    // Ensure source is empty after move
    REQUIRE(source_for_move_assign.get_points().empty());
    std::cout << "Move assignment operator: OK" << std::endl;

    std::cout << "All irrterp<double, double> tests passed!" << std::endl << std::endl;
}


