#include "WireCellSpng/RayTest.h"

#include <stdexcept> // For std::invalid_argument
#include <algorithm> // max/min

namespace WireCell::Spng::RayGrid {

    torch::Tensor get_gcd(const std::string& key)
    {
        // abbrev
        const auto& gcd = WireCell::Spng::RayGrid::ray_grid_coordinates_data;
        auto it = gcd.find(key);
        if (it == gcd.end()) {
            return torch::tensor({});
        }
        return it->second;
    }


    Coordinates gcd_coordinates()
    {
        auto views = get_gcd("views");
        return Coordinates(views);
    }

// Helper for comparing tensors
bool are_tensors_close(const torch::Tensor& a, const torch::Tensor& b, double rtol, double atol) {
    if (a.sizes() != b.sizes()) {
        std::cerr << "Tensor sizes mismatch: " << a.sizes() << " vs " << b.sizes() << std::endl;
        return false;
    }
    bool is_close = torch::allclose(a, b, rtol, atol); //.item<bool>();
    if (!is_close) {
        std::cerr << "Tensors are not close:\na=" << a << "\nb=" << b << std::endl;
        return false;
    }
    return true;
}

    // This is dumped from Python like:
    // >>> views = symmetric_views()
    // >>> coords = wirecell.raygrid.coordinates.Coordinates(views)
    // >>> print(coords.as_string("c++"))
    //
    // Or run uv run --with pytest pytest -k 'test_dump_coordinates' wirecell/raygrid/test/test_coordinates.py
    // and insert the generated file test-dump-coordinates.cpp
const std::map<std::string, torch::Tensor> ray_grid_coordinates_data = {
    { "views", torch::tensor({{{ 50.00000000000000000,   0.00000000000000000},
         { 50.00000000000000000, 100.00000000000000000}},

        {{  0.00000000000000000,  50.00000000000000000},
         {100.00000000000000000,  50.00000000000000000}},

        {{  0.75000000000000000,  98.70095825195312500},
         {  2.25000000000000000,  96.10288238525390625}},

        {{  0.75000000000000000,   1.29903805255889893},
         {  2.25000000000000000,   3.89711427688598633}},

        {{  0.00000000000000000,  50.00000000000000000},
         {  3.00000000000000000,  50.00000000000000000}}}, torch::kDouble) },
    { "pitch_mag", torch::tensor({100.00000000000000000, 100.00000000000000000,   2.99999976158142090,
          3.00000000000000000,   3.00000000000000000}, torch::kDouble) },
    { "pitch_dir", torch::tensor({{ 0.00000000000000000,  1.00000000000000000},
        { 1.00000000000000000,  0.00000000000000000},
        { 0.50000005960464478, -0.86602532863616943},
        { 0.50000000000000000,  0.86602544784545898},
        { 1.00000000000000000,  0.00000000000000000}}, torch::kDouble) },
    { "center", torch::tensor({{50.00000000000000000,  0.00000000000000000},
        { 0.00000000000000000, 50.00000000000000000},
        { 0.75000000000000000, 98.70095825195312500},
        { 0.75000000000000000,  1.29903805255889893},
        { 0.00000000000000000, 50.00000000000000000}}, torch::kDouble) },
    { "zero_crossings", torch::tensor({{{  50.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000},
         {-170.20506286621093750,    0.00000000000000000},
         {   3.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000}},

        {{   0.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,   50.00000000000000000},
         {   0.00000000000000000,   98.26794433593750000},
         {   0.00000000000000000,    1.73205184936523438},
         {   0.00000000000000000,    0.00000000000000000}},

        {{-170.20506286621093750,    0.00000000000000000},
         {   0.00000000000000000,   98.26794433593750000},
         {   0.75000000000000000,   98.70095825195312500},
         { -83.60253906250000000,   49.99999618530273438},
         {   0.00000000000000000,   98.26794433593750000}},

        {{   3.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,    1.73205184936523438},
         { -83.60253906250000000,   49.99999618530273438},
         {   0.75000000000000000,    1.29903805255889893},
         {   0.00000000000000000,    1.73205065727233887}},

        {{   0.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,   98.26794433593750000},
         {   0.00000000000000000,    1.73205065727233887},
         {   0.00000000000000000,   50.00000000000000000}}},
       torch::kDouble) },
    { "ray_jump", torch::tensor({{{  -1.00000000000000000,    0.00000000000000000},
         { 100.00000000000000000,    0.00000000000000000},
         {   5.99998474121093750,    0.00000000000000000},
         {   6.00000000000000000,    0.00000000000000000},
         {   3.00000000000000000,    0.00000000000000000}},

        {{   0.00000000000000000,  100.00000000000000000},
         {   0.00000000000000000,    1.00000000000000000},
         {   0.00000000000000000,   -3.46409606933593750},
         {   0.00000000000000000,    3.46410369873046875},
         {   0.00000000000000000,    0.00000000000000000}},

        {{ 173.20506286621093750,  100.00000000000000000},
         { 100.00000000000000000,   57.73503112792968750},
         {   0.86602538824081421,    0.50000005960464478},
         {   3.00000762939453125,    1.73205184936523438},
         {   3.00000000000000000,    1.73205566406250000}},

        {{-173.20507812500000000,  100.00000000000000000},
         { 100.00000000000000000,  -57.73502349853515625},
         {   3.00000762939453125,   -1.73205184936523438},
         {  -0.86602544784545898,    0.50000000000000000},
         {   3.00000000000000000,   -1.73205065727233887}},

        {{   0.00000000000000000,  100.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000},
         {   0.00000000000000000,   -3.46409606933593750},
         {   0.00000000000000000,    3.46410489082336426},
         {   0.00000000000000000,    1.00000000000000000}}},
       torch::kDouble) },
    { "ray_dir", torch::tensor({{-1.00000000000000000,  0.00000000000000000},
        {-0.00000000000000000,  1.00000000000000000},
        { 0.86602532863616943,  0.50000005960464478},
        {-0.86602544784545898,  0.50000000000000000},
        {-0.00000000000000000,  1.00000000000000000}}, torch::kDouble) },
    { "a", torch::tensor({{{   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
            50.00000762939453125,   50.00000000000000000,
           100.00000000000000000},
         {   0.00000000000000000,    5.99998474121093750,
             0.00000000000000000,    2.99999237060546875,
             5.99998474121093750},
         {   0.00000000000000000,    6.00000000000000000,
             3.00000047683715820,    0.00000000000000000,
             6.00000000000000000},
         {   0.00000000000000000,    3.00000000000000000,
             1.50000023841857910,    1.50000000000000000,
             0.00000000000000000}},

        {{   0.00000000000000000,    0.00000000000000000,
           -86.60253143310546875,   86.60254669189453125,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {  -3.46409606933593750,    0.00000000000000000,
             0.00000000000000000,   -2.99999523162841797,
             0.00000000000000000},
         {   3.46410369873046875,    0.00000000000000000,
            -3.00000143051147461,    0.00000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000}},

        {{   0.00000000000000000,  173.20506286621093750,
             0.00000000000000000,  173.20507812500000000,
           173.20506286621093750},
         {  57.73503112792968750,    0.00000000000000000,
             0.00000000000000000,  100.00000762939453125,
           100.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {   1.73205184936523438,    3.00000762939453125,
             0.00000000000000000,    0.00000000000000000,
             3.00000762939453125},
         {   1.73205566406250000,    3.00000000000000000,
             0.00000000000000000,    3.00000429153442383,
             0.00000000000000000}},

        {{   0.00000000000000000, -173.20507812500000000,
          -173.20507812500000000,    0.00000000000000000,
          -173.20507812500000000},
         { -57.73502349853515625,    0.00000000000000000,
           100.00000000000000000,    0.00000000000000000,
           100.00000000000000000},
         {  -1.73205184936523438,    3.00000762939453125,
             0.00000000000000000,    0.00000000000000000,
             3.00000762939453125},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {  -1.73205065727233887,    3.00000000000000000,
             3.00000000000000000,    0.00000000000000000,
             0.00000000000000000}},

        {{   0.00000000000000000,    0.00000000000000000,
           -86.60253143310546875,   86.60254669189453125,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {  -3.46409606933593750,    0.00000000000000000,
             0.00000000000000000,   -2.99999523162841797,
             0.00000000000000000},
         {   3.46410489082336426,    0.00000000000000000,
            -3.00000262260437012,    0.00000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000}}}, torch::kDouble) },
    { "b", torch::tensor({{{   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
            85.10253143310546875,   -1.50000000000000000,
             0.00000000000000000},
         {   0.00000000000000000, -170.20506286621093750,
             0.00000000000000000,  -86.60253143310546875,
          -170.20506286621093750},
         {   0.00000000000000000,    3.00000000000000000,
            86.60253143310546875,    0.00000000000000000,
             3.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
            85.10253143310546875,   -1.50000000000000000,
             0.00000000000000000}},

        {{   0.00000000000000000,    0.00000000000000000,
            85.10253143310546875,   -1.50000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {  98.26794433593750000,    0.00000000000000000,
             0.00000000000000000,   83.60253906250000000,
             0.00000000000000000},
         {   1.73205184936523438,    0.00000000000000000,
            83.60253143310546875,    0.00000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
            85.10253143310546875,   -1.50000000000000000,
             0.00000000000000000}},

        {{   0.00000000000000000, -170.20506286621093750,
             0.00000000000000000,  -86.60253143310546875,
          -170.20506286621093750},
         {  98.26794433593750000,    0.00000000000000000,
             0.00000000000000000,   83.60253906250000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {  49.99999618530273438,  -83.60253906250000000,
             0.00000000000000000,    0.00000000000000000,
           -83.60253906250000000},
         {  98.26794433593750000,    0.00000000000000000,
             0.00000000000000000,   83.60253906250000000,
             0.00000000000000000}},

        {{   0.00000000000000000,    3.00000000000000000,
            86.60253143310546875,    0.00000000000000000,
             3.00000000000000000},
         {   1.73205184936523438,    0.00000000000000000,
            83.60253143310546875,    0.00000000000000000,
             0.00000000000000000},
         {  49.99999618530273438,  -83.60253906250000000,
             0.00000000000000000,    0.00000000000000000,
           -83.60253906250000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000},
         {   1.73205065727233887,    0.00000000000000000,
            83.60253143310546875,    0.00000000000000000,
             0.00000000000000000}},

        {{   0.00000000000000000,    0.00000000000000000,
            85.10253143310546875,   -1.50000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
            85.10253143310546875,   -1.50000000000000000,
             0.00000000000000000},
         {  98.26794433593750000,    0.00000000000000000,
             0.00000000000000000,   83.60253906250000000,
             0.00000000000000000},
         {   1.73205065727233887,    0.00000000000000000,
            83.60253143310546875,    0.00000000000000000,
             0.00000000000000000},
         {   0.00000000000000000,    0.00000000000000000,
             0.00000000000000000,    0.00000000000000000,
             0.00000000000000000}}}, torch::kDouble) }
};

torch::Tensor random_groups(
    int64_t ngroups,
    int64_t n_in_group,
    double sigma,
    std::array<double, 2> ll,
    std::array<double, 2> ur,
    torch::TensorOptions options) {



    // 1. Generate `ngroups` "group points" uniformly random within the bounding box.
    double x_range = ur[0] - ll[0];
    double y_range = ur[1] - ll[1];

    if (x_range < 0 || y_range < 0) {
        throw std::invalid_argument("Upper-right corner must be greater than or equal to lower-left corner.");
    }

    // Generate points in [0, 1) using `options` for dtype and device
    torch::Tensor group_means = torch::rand({ngroups, 2}, options);

    // Scale and shift the group means using broadcasting
    // Create tensors for ranges and lower-left corner, ensuring they match `options`
    torch::Tensor range_tensor = torch::tensor({x_range, y_range}, options);
    torch::Tensor ll_tensor = torch::tensor({ll[0], ll[1]}, options);

    group_means = group_means * range_tensor + ll_tensor;

    // 2. For each group point, generate `n_in_group` points
    // according to a Gaussian distribution centered at the group point.

    // Generate standard normal noise: (ngroups, n_in_group, 2)
    // The options are reused to ensure noise is on the same device and has the same dtype
    torch::Tensor noise = torch::randn({ngroups, n_in_group, 2}, options);

    // Scale the noise by sigma
    torch::Tensor scaled_noise = noise * sigma;

    // Add the group means to the scaled noise.
    // group_means (ngroups, 2) needs to be unsqueezed to (ngroups, 1, 2)
    // for broadcasting with scaled_noise (ngroups, n_in_group, 2).
    torch::Tensor points_in_groups = group_means.unsqueeze(1) + scaled_noise;

    return points_in_groups;
}

std::vector<torch::Tensor> fill_activity(const Coordinates& coords, const torch::Tensor& points) {
    int64_t nviews = coords.nviews();
    std::vector<torch::Tensor> activity_tensors;

    // Flatten points to (num_points, 2) if it's (ngroups, n_in_group, 2)
    torch::Tensor flattened_points = points;
    if (points.dim() == 3 && points.size(2) == 2) {
        flattened_points = points.reshape({-1, 2});
    } else if (points.dim() != 2 || points.size(1) != 2) {
        throw std::invalid_argument("Input 'points' tensor must have shape (N, 2) or (N, M, 2).");
    }

    // Get all pitch indices for all points across all views in one vectorized call
    // Shape (num_points, nviews)
    torch::Tensor all_pitch_indices = coords.point_indices(flattened_points);

    // std::cerr << "indices: " << all_pitch_indices << "\n";

    // Use the torch::indexing namespace for Slice and Ellipsis
    using namespace torch::indexing;

    auto bogus = torch::zeros({}, points.options().dtype(torch::kBool));

    auto ab = coords.active_bounds();

    // std::cerr << "active bounds: " << ab << "\n";
        
    // Iterate through each view to fill its individual activity tensor
    for (int64_t view_idx = 0; view_idx < nviews; ++view_idx) {

        const int64_t b_beg = std::max(0L, ab.index({view_idx, 0}).item<int64_t>());
        const int64_t b_end = std::max(0L, ab.index({view_idx, 1}).item<int64_t>());
        if (b_beg >= b_end) {
            std::cerr << "view " << view_idx << " has no active bounds: [" << b_beg << ", " << b_end << "]\n";
            throw std::invalid_argument("view has no active bounds");
        }

        torch::Tensor current_view_indices = all_pitch_indices.index({Slice(), view_idx});

        // Handle case where no points are present for this view (e.g., if flattened_points was empty)
        if (current_view_indices.numel() == 0) {
            // If no points, the activity tensor should be an empty boolean tensor.
            std::cerr << "no points in view " << view_idx << "\n";
            activity_tensors.push_back(bogus);
            continue; // Move to the next view
        }
        
        {
            // Calculate the minimum and maximum pitch indices encountered for
            // this view as a check, bail if activity is bogus.
            const int64_t min_idx = std::max(0L, current_view_indices.min().item<int64_t>());
            const int64_t max_idx = std::max(0L, current_view_indices.max().item<int64_t>());

            // Possible that both were negative, ie all points are outside bounds
            if (max_idx < 0) {
                activity_tensors.push_back(bogus);
                std::cerr << "fill_activity: all points below view " << view_idx << " minmax=["<<min_idx<<","<<max_idx<<"]\n";
                continue; // Move to the next view
            }
            if (min_idx > max_idx) {
                activity_tensors.push_back(bogus);
                std::cerr << "fill_activity: all points miss in view " << view_idx << " minmax=["<<min_idx<<","<<max_idx<<"]\n";
                continue; // Move to the next view
            }
        }        

        int64_t activity_tensor_size = b_end; // b_end is "hi" side of half-open range

        // Create the activity tensor initialized to False (all zeros for bool)
        torch::Tensor activity_tensor = torch::zeros({activity_tensor_size}, points.options().dtype(torch::kBool));

        torch::Tensor clamped = torch::clamp(current_view_indices, 0, activity_tensor_size - 1);

        // Set the corresponding positions in the activity tensor to True.
        // This is equivalent to Python's `activity_tensor[indices] = True`.
        // `index_put_` requires the value to be a tensor.
        activity_tensor.index_put_({clamped}, torch::tensor(true, points.options().dtype(torch::kBool)));
        activity_tensors.push_back(activity_tensor);
    }
    return activity_tensors;
}

torch::Tensor symmetric_views(
    double width,
    double height,
    double pitch_mag,
    double angle, // Already in radians
    torch::TensorOptions options) {

    // Initialize pitches tensor with the specified options
    torch::Tensor pitches = torch::zeros({5, 2, 2}, options);

    // horizontal ray bounds, pitch is vertical
    // pitches[0] = torch.tensor([[width/2.0, 0], [width/2.0, height]])
    pitches.index_put_({0}, torch::tensor(
        {{width/2.0, 0.0},
         {width/2.0, height}}, options));

    // vertical ray bounds, pitch is horizontal
    // pitches[1] = torch.tensor([[0, height/2.0], [width, height/2.0]])
    pitches.index_put_({1}, torch::tensor(
        {{0.0, height/2.0},
         {width, height/2.0}}, options));

    // corners
    // Ensure these tensors are created with the same options as 'pitches'
    torch::Tensor ll = torch::tensor({0.0, 0.0}, options);
    torch::Tensor ul = torch::tensor({0.0, height}, options);
    torch::Tensor lr = torch::tensor({width, 0.0}, options);

    // /-wires
    // w = torch.tensor([math.sin(angle), math.cos(angle)])
    torch::Tensor w_slash = torch::tensor({std::sin(angle), std::cos(angle)}, options);
    // p = torch.tensor([w[1], -w[0]])
    // Accessing elements of w_slash and creating a new tensor for p_slash
    torch::Tensor p_slash = torch::tensor({w_slash.index({1}).item<double>(), -w_slash.index({0}).item<double>()}, options);

    // pitches[2] = torch.vstack([ul + 0.5*pitch_mag*p, ul + 1.5*pitch_mag*p])
    pitches.index_put_({2}, torch::vstack(
        {ul + 0.5 * pitch_mag * p_slash,
         ul + 1.5 * pitch_mag * p_slash}));

    // the symmetry (angle *= -1)
    angle *= -1.0; // Modify the angle directly for the next calculation

    // \-wires
    // w = torch.tensor([math.sin(angle), math.cos(angle)])
    torch::Tensor w_backslash = torch::tensor({std::sin(angle), std::cos(angle)}, options);
    // p = torch.tensor([w[1], -w[0]])
    // Accessing elements of w_backslash and creating a new tensor for p_backslash
    torch::Tensor p_backslash = torch::tensor({w_backslash.index({1}).item<double>(), -w_backslash.index({0}).item<double>()}, options);

    // pitches[3] = torch.vstack([ll + 0.5*pitch_mag*p, ll + 1.5*pitch_mag*p])
    pitches.index_put_({3}, torch::vstack(
        {ll + 0.5 * pitch_mag * p_backslash,
         ll + 1.5 * pitch_mag * p_backslash}));

    // |-wires
    // pitches[4] = torch.tensor([[0, height/2.0], [pitch_mag, height/2.0]])
    pitches.index_put_({4}, torch::tensor(
        {{0.0, height/2.0},
         {pitch_mag, height/2.0}}, options));

    return pitches;
}



} // namespace
