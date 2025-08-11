// test_raygrid.cpp
#include "WireCellSpng/RayGrid.h"
#include "WireCellSpng/RayTest.h"
#include <iostream>

// Name collission for "CHECK" between torch and doctest.
#undef CHECK
#include "WireCellUtil/doctest.h"

using namespace torch::indexing;
using namespace WireCell::Spng::RayGrid;

TEST_CASE("spng ray grid three view") {
    std::cerr << "Running RayGrid tests..." << std::endl;

    // Test data (example from Python, adjust as needed)
    // This example uses 3 views for simplicity, the original Python code used 5.
    // (N-views, 2 endpoints, 2 coordinates)
    torch::Tensor views = torch::tensor({
        {{0.0, 0.0}, {1.0, 0.0}}, // View 0: horizontal, pitch along x-axis
        {{0.0, 0.0}, {0.0, 1.0}}, // View 1: vertical, pitch along y-axis
        {{0.0, 0.0}, {0.70710678, 0.70710678}} // View 2: diagonal, pitch along y=x
    }, torch::kDouble);

    // Test constructor and basic attributes
    Coordinates coords(views);
    std::cerr << "nviews: " << coords.nviews() << std::endl;
    REQUIRE(coords.nviews() == 3);

    std::cerr << "pitch_mag: " << coords.pitch_mag << std::endl;
    std::cerr << "pitch_dir: " << coords.pitch_dir << std::endl;
    std::cerr << "center: " << coords.center << std::endl;
    std::cerr << "ray_dir: " << coords.ray_dir << std::endl;
    std::cerr << "zero_crossings: " << coords.zero_crossings << std::endl;
    std::cerr << "ray_jump: " << coords.ray_jump << std::endl;
    std::cerr << "a: " << coords.a << std::endl;
    std::cerr << "b: " << coords.b << std::endl;

    // Test bounding_box
    torch::Tensor expected_bbox = torch::tensor({{0.0, 1.0}, {0.0, 1.0}}, torch::kDouble); // For the first two views
    // The bounding box logic in Python is simplified and assumes specific view orientations.
    // For general views, this might need more robust calculation.
    // For the given example, pitch_dir[0,0] is 1.0 (not 0), so it goes to the else branch.
    // x0 = center[0,1] = 0.0
    // y0 = center[1,0] = 0.0
    // x1 = x0 + pitch_mag[0] = 0.0 + 1.0 = 1.0
    // y1 = y0 + pitch_mag[1] = 0.0 + 1.0 = 1.0
    torch::Tensor actual_bbox = coords.bounding_box();
    std::cerr << "bounding_box: " << actual_bbox << std::endl;
    REQUIRE(are_tensors_close(actual_bbox, expected_bbox));


    // Test point_pitches
    torch::Tensor test_points = torch::tensor({
        {0.5, 0.5},
        {0.1, 0.9},
        {1.0, 0.0}
    }, torch::kDouble);
    torch::Tensor pitches = coords.point_pitches(test_points);
    std::cerr << "point_pitches:\n" << pitches << std::endl;
    // Expected pitches for the given views:
    // Point (0.5, 0.5):
    // View 0 (x-axis pitch): 0.5
    // View 1 (y-axis pitch): 0.5
    // View 2 (diagonal, pitch_dir = (0.707, 0.707)): dot((0.5,0.5), (0.707,0.707)) = 0.707
    // Point (0.1, 0.9):
    // View 0: 0.1
    // View 1: 0.9
    // View 2: dot((0.1,0.9), (0.707,0.707)) = 0.0707 + 0.6363 = 0.707
    // Point (1.0, 0.0):
    // View 0: 1.0
    // View 1: 0.0
    // View 2: dot((1.0,0.0), (0.707,0.707)) = 0.707
    torch::Tensor expected_pitches = torch::tensor({
        {0.5, 0.5, 0.70710678},
        {0.1, 0.9, 0.70710678},
        {1.0, 0.0, 0.70710678}
    }, torch::kDouble);
    REQUIRE(are_tensors_close(pitches, expected_pitches));


    // Test point_indices
    torch::Tensor indices = coords.point_indices(test_points);
    std::cerr << "point_indices:\n" << indices << std::endl;
    // Expected indices (floor of pitches, since pitch_mag is 1 for views 0 and 1, and 1 for view 2)
    torch::Tensor expected_indices = torch::tensor({
        {0, 0, 0},
        {0, 0, 0},
        {1, 0, 0}
    }, torch::kLong);
    REQUIRE(are_tensors_close(indices.to(torch::kDouble), expected_indices.to(torch::kDouble))); // Compare as double for allclose


    // Test ray_crossing (scalar input)
    {
        torch::Tensor one_scalar = torch::tensor({0, 1}, torch::kLong); // view 0, ray 1
        torch::Tensor two_scalar = torch::tensor({1, 1}, torch::kLong); // view 1, ray 1
        torch::Tensor crossing_point_scalar = coords.ray_crossing(
            one_scalar[0], one_scalar[1],
            two_scalar[0], two_scalar[1]);
        std::cerr << "ray_crossing (scalar): " << crossing_point_scalar << std::endl;
        REQUIRE(are_tensors_close(crossing_point_scalar, torch::tensor({1.0, 1.0}, torch::kDouble)));
    }
    {
        torch::Tensor one_scalar = torch::tensor({0, 1}, torch::kLong); // view 0, ray 1
        torch::Tensor two_scalar = torch::tensor({1, 0}, torch::kLong); // view 1, ray 0
        torch::Tensor crossing_point_scalar = coords.ray_crossing(
            one_scalar[0], one_scalar[1],
            two_scalar[0], two_scalar[1]);
        std::cerr << "ray_crossing (scalar): " << crossing_point_scalar << std::endl;
        REQUIRE(are_tensors_close(crossing_point_scalar, torch::tensor({1.0, 0.0}, torch::kDouble)));
    }

    // Expose this one to use below
    torch::Tensor one_scalar = torch::tensor({0, 0}, torch::kLong); // view 0, ray 0
    torch::Tensor two_scalar = torch::tensor({1, 0}, torch::kLong); // view 1, ray 0
    torch::Tensor crossing_point_scalar = coords.ray_crossing(
            one_scalar[0], one_scalar[1],
            two_scalar[0], two_scalar[1]);
    std::cerr << "ray_crossing (scalar): " << crossing_point_scalar << std::endl;
    REQUIRE(are_tensors_close(crossing_point_scalar, torch::tensor({0.0, 0.0}, torch::kDouble)));


    // Test ray_crossing (batched input)
    torch::Tensor one_batch = torch::tensor({
            {0, 0}, {0, 1}, {0, 1}
    }, torch::kLong);
    torch::Tensor two_batch = torch::tensor({
            {1, 0}, {1, 0}, {1, 1}
    }, torch::kLong);
    torch::Tensor crossing_point_batch = coords.ray_crossing(
            one_batch.index({Slice(), 0}), one_batch.index({Slice(), 1}),
            two_batch.index({Slice(), 0}), two_batch.index({Slice(), 1}));

    std::cerr << "ray_crossing (batch):\n" << crossing_point_batch << std::endl;
    // (view 0, ray 0) and (view 1, ray 0) -> (0,0)
    // (view 0, ray 1) and (view 1, ray 0) -> (1,0) (ray 1 of view 0 is x=1, ray 0 of view 1 is y=0)
    torch::Tensor expected_crossing_batch = torch::tensor({
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0}
    }, torch::kDouble);
    REQUIRE(are_tensors_close(crossing_point_batch, expected_crossing_batch));


    // Test pitch_location (scalar input)
    torch::Tensor view_idx_scalar = torch::tensor(2, torch::kLong); // View 2
    torch::Tensor pitch_loc_scalar = coords.pitch_location(
        one_scalar[0], one_scalar[1],
        two_scalar[0], two_scalar[1], view_idx_scalar);
    std::cerr << "pitch_location (scalar): " << pitch_loc_scalar << std::endl;
    // For (0,0) and (1,0) crossing (which is (0,0)), its pitch in view 2 (diagonal) is 0.
    REQUIRE(are_tensors_close(pitch_loc_scalar, torch::tensor(0.0, torch::kDouble)));

    // Test pitch_location (batched input)
    torch::Tensor view_idx_batch = torch::tensor({2, 2, 2}, torch::kLong); // View 2 for both
    torch::Tensor pitch_loc_batch = coords.pitch_location(
            one_batch.index({Slice(), 0}), one_batch.index({Slice(), 1}),
            two_batch.index({Slice(), 0}), two_batch.index({Slice(), 1}),
            view_idx_batch);
    std::cerr << "pitch_location (batch):\n" << pitch_loc_batch << std::endl;
    // For (0,0) and (1,0) crossing (0,0), pitch in view 2 is 0.
    // For (0,1) and (1,0) crossing (1,0), pitch in view 2 is dot((1,0), (0.707,0.707)) = 0.707
    torch::Tensor expected_pitch_loc_batch = torch::tensor({0.0, 0.70710678, 2*0.70710678}, torch::kDouble);
    REQUIRE(are_tensors_close(pitch_loc_batch, expected_pitch_loc_batch));


    // Test pitch_index (scalar input)
    torch::Tensor test_pitch_scalar = torch::tensor(0.75, torch::kDouble);
    torch::Tensor test_view_scalar = torch::tensor(0, torch::kLong);
    torch::Tensor p_idx_scalar = coords.pitch_index(test_pitch_scalar, test_view_scalar);
    std::cerr << "pitch_index (scalar): " << p_idx_scalar << std::endl;
    // pitch_mag[0] is 1.0, so floor(0.75/1.0) = 0
    REQUIRE(p_idx_scalar.item<long>() == 0);

    // Test pitch_index (batched pitch input)
    torch::Tensor test_pitch_batch = torch::tensor({0.2, 1.5, 0.9}, torch::kDouble);
    torch::Tensor test_view_batch = torch::tensor({0, 0, 1}, torch::kLong); // View 0, View 0, View 1
    torch::Tensor p_idx_batch = coords.pitch_index(test_pitch_batch, test_view_batch);
    std::cerr << "pitch_index (batch):\n" << p_idx_batch << std::endl;
    // For pitch 0.2, view 0 (mag 1.0) -> floor(0.2/1.0) = 0
    // For pitch 1.5, view 0 (mag 1.0) -> floor(1.5/1.0) = 1
    // For pitch 0.9, view 1 (mag 1.0) -> floor(0.9/1.0) = 0
    torch::Tensor expected_p_idx_batch = torch::tensor({0, 1, 0}, torch::kLong);
    REQUIRE(are_tensors_close(p_idx_batch.to(torch::kDouble), expected_p_idx_batch.to(torch::kDouble)));


    std::cerr << "All RayGrid tests passed!" << std::endl;

}


TEST_CASE("spng ray grid canned data") 
{
    auto views = get_gcd("views");
    Coordinates coords(views);
    
    auto a = get_gcd("a");
    auto da = a - coords.a;

    std::cerr << "The 'a' tensor from Python:" << a << "\n";
    std::cerr << "The 'a' tensor from C++:" << coords.a << "\n";
    std::cerr << "The difference:" << da << "\n";

    REQUIRE( are_tensors_close(a, coords.a) );

    auto b = get_gcd("b");
    REQUIRE( are_tensors_close(b, coords.b) );

}
