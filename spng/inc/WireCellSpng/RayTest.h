#ifndef WIRECELL_SPNG_RAYTEST
#define WIRECELL_SPNG_RAYTEST

#include <map>
#include <string>
#include "WireCellSpng/Torch.h"

#include "RayGrid.h"

namespace WireCell::Spng::RayGrid {

    // Example data members dumped from the Python version of Coordinates for a
    // microboone-like detector.
    extern const std::map<std::string, torch::Tensor> ray_grid_coordinates_data;


    /// Helper to get one key from the const ray_grid_coordinates_data.
    torch::Tensor get_gcd(const std::string& key);


    /// Return a coordinates made from the ray_grid_coordinates_data.
    Coordinates gcd_coordinates();

    // Helper for comparing tensors
    bool are_tensors_close(const torch::Tensor& a, const torch::Tensor& b, double rtol = 1e-05, double atol = 1e-08);

    /**
     * @brief Generates "clumpy" 2D points by first creating group centers uniformly
     * and then sampling points around these centers from a Gaussian distribution.
     *
     * @param ngroups The number of distinct groups (clusters) of points.
     * @param n_in_group The number of points to generate within each group.
     * @param sigma The standard deviation (width) of the Gaussian distribution
     * for points within each group.
     * @param ll A std::array (x_min, y_min) representing the lower-left corner
     * of the bounding box for generating group centers.
     * @param ur A std::array (x_max, y_max) representing the upper-right corner
     * of the bounding box for generating group centers.
     * @param options TensorOptions to control the dtype and device of the output tensor.
     * Defaults to float32 on CPU.
     * @return A torch::Tensor of shape (ngroups, n_in_group, 2) containing
     * the generated 2D points. Points outside the bounding box are
     * included as per the requirement.
     */
    torch::Tensor random_groups(
        int64_t ngroups = 10,
        int64_t n_in_group = 10,
        double sigma = 10.0f,
        std::array<double, 2> ll = {0.0f, 0.0f},
        std::array<double, 2> ur = {100.0f, 100.0f},
        torch::TensorOptions options = torch::TensorOptions().dtype(torch::kDouble));

    /**
     * @brief Produces a list of "activity" tensors by projecting a tensor of points.
     *
     * This function processes a tensor of 2D points, projecting them onto
     * different "views" (as defined by the `coords` object) and generating
     * boolean activity tensors for each view.
     *
     * @param coords An object providing view information and point indexing capabilities.
     * @param points A torch::Tensor of 2D points. Expected shape is (N, 2) or (N, M, 2).
     * If (N, M, 2), it will be internally reshaped to (N*M, 2).
     * @return A std::vector of torch::Tensor, where each tensor represents
     * the activity (boolean presence) for a specific view.
     */
    std::vector<torch::Tensor> fill_activity(const Coordinates& coords, const torch::Tensor& points);


    /**
     * @brief Returns a (5,2,2) pitch ray tensor suitable for giving to Coordinates().
     *
     * Angle is from y-axis to wire direction.
     *
     * This function is equivalent to symmetric_raypairs() in WCT's
     * util/src/RayHelpers.cxx except it returns a pitch ray instead of a pair of
     * wire rays and it operates in 2D (x_rg, y_rg) instead of 3D (x_wc, y_wc,
     * z_wc). The coordinates correspondence is:
     *
     * y_rg == y_wc
     * x_rg == z_wc
     * z_rg == -x_wc
     *
     * Note, z_rg is not used here but may be defined by right-hand-rule. Also by
     * RHR, "pitch cross wire = z_rg" (WCT it's opposite: "wire cross pitch = x_wc").
     *
     * Rays generally start on left hand side, pitch increasing toward right
     * (x_rg). This is equivalent to increasing with z_wc. Special case
     * horizontal rays start at bottom and go up.
     *
     * @param width The width of the detector area.
     * @param height The height of the detector area.
     * @param pitch_mag The magnitude of the pitch (e.g., wire pitch).
     * @param angle The angle in radians, from y-axis to wire direction.
     * @param options TensorOptions to control the dtype and device of the output tensor.
     * Defaults to double on CPU.
     * @return A torch::Tensor of shape (5, 2, 2) representing pitch rays.
     */
    torch::Tensor symmetric_views(
        double width = 100.0f,
        double height = 100.0f,
        double pitch_mag = 3.0f,
        double angle = 1.0471975511965976, // math.radians(60.0)
        torch::TensorOptions options = torch::TensorOptions().dtype(torch::kDouble));
    



} // namespace

#endif
