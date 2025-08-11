// RayGrid.h
#ifndef WIRECELL_SPNG_RAYGRID_RAYGRID_H
#define WIRECELL_SPNG_RAYGRID_RAYGRID_H

#include <torch/torch.h>
#include <stdexcept> // For std::runtime_error

namespace WireCell {
namespace Spng {
namespace RayGrid {

// --- Helper functions from funcs.py ---

/**
 * @brief Returns the perpendicular vector giving the separation between two parallel rays.
 *
 * @param r0 A tensor representing the first ray, shape (2, 2).
 * @param r1 A tensor representing the second ray, shape (2, 2).
 * @return A tensor representing the pitch vector.
 */
torch::Tensor pitch(const torch::Tensor& r0, const torch::Tensor& r1);

/**
 * @brief Returns the 2D crossing point of two non-parallel rays.
 *
 * @param r0 A tensor representing the first ray, shape (2, 2).
 * @param r1 A tensor representing the second ray, shape (2, 2).
 * @return A tensor representing the intersection point, shape (2).
 * @throws std::runtime_error if the lines are parallel.
 */
torch::Tensor crossing(const torch::Tensor& r0, const torch::Tensor& r1);

/**
 * @brief Returns the vector along the ray direction.
 *
 * @param ray A tensor representing the ray, shape (2, 2).
 * @return A tensor representing the direction vector, shape (2).
 */
torch::Tensor vector(const torch::Tensor& ray);

/**
 * @brief Returns the unit vector along the ray direction.
 *
 * @param ray A tensor representing the ray, shape (2, 2).
 * @return A tensor representing the unit direction vector, shape (2).
 */
torch::Tensor ray_direction(const torch::Tensor& ray);


// --- Coordinates class ---

class Coordinates {

public:
    // Public tensor attributes.  All are dtype kDouble.

    // (Nview,) magnitude of pitch of each view
    torch::Tensor pitch_mag;
    // (Nview, 2) unit vector along the pitch direction of each view
    torch::Tensor pitch_dir;
    // (Nview, 2) origin vector of each view
    torch::Tensor center;
    // (Nview, Nview, 2) crossing point of a "ray zero" from a pair of views.
    // Undefined for views which are mutually parallel.
    torch::Tensor zero_crossings;
    // (Nview, Nview, 2) displacement vector along ray direction of the first
    // view between two crossing of that ray and two consecutive rays in the
    // second view.
    torch::Tensor ray_jump;
    // The ray grid "A" tensor.  See raygrid.pdf doc.
    torch::Tensor a;
    // The ray grid "B" tensor.  See raygrid.pdf doc.
    torch::Tensor b;

    // (Nview,2), Per-view unit vectors in direction of ray_jump.  Helper, not
    // needed for main ray-grid coordinate calculations
    torch::Tensor ray_dir;

    /**
     * @brief Constructs Ray Grid coordinates specified by views.
     *
     * The views is a 3-D tensor of shape:
     * (N-views, 2 endpoints, 2 coordinates)
     * Each view is a pair of endpoints. The first point marks the origin of
     * the view. The relative vector from first to second point is in the
     * direction of the pitch. The magnitude of the vector is the pitch.
     *
     * @param views A tensor of shape (Nviews, 2, 2) representing the views.
     */
    Coordinates(const torch::Tensor& views);

    /**
     * @brief Move tensors to a device.
     *
     * Tensors are created on the CPU device but after construction can be moved.
     */
    void to(torch::Device);

    /**
     * @brief Returns the number of views.
     * @return The number of views.
     */
    long nviews() const;

    /**
     * @brief Returns the bounding box of the ray grid.
     * @return A tensor of shape (2, 2) holding [ (x0,x1), (y0,y1) ] bounds in Cartesian space.
     */
    torch::Tensor bounding_box() const;

    /**
     * @brief Returns the pitch location measured in each view for a batch of 2D Cartesian points.
     *
     * @param points A tensor of shape (nbatch, 2) providing 2D Cartesian coordinates.
     * @return A tensor of floating point values and shape (nbatch, nview)
     * giving the per-view pitch for each point.
     * @throws std::invalid_argument if 'points' is not of shape (nbatch, 2).
     */
    torch::Tensor point_pitches(const torch::Tensor& points) const;

    /**
     * @brief Returns the integer pitch index in each view for a batch of 2D Cartesian points.
     *
     * @param points A tensor of shape (nbatch, 2) providing 2D Cartesian coordinates.
     * @return A tensor of integer values and shape (nbatch, nviews)
     * giving the pitch index in each view for each point.
     */
    torch::Tensor point_indices(const torch::Tensor& points) const;

    /**
     * @brief Returns the 2D crossing point(s) of ray grid coordinates "one" and "two".
     *
     * Each coordinate is given as a pair (view,ray) of indices. These may be scalar or batched array.
     *
     * @param view1 A scalar tensor or of shape (nbatch,) giving view indices for first ray.
     * @param ray1 A scalar tensor or of shape (nbatch,) giving ray indices for first ray.
     * @param view2 A scalar tensor or of shape (nbatch,) giving view indices for second ray.
     * @param ray2 A scalar tensor or of shape (nbatch,) giving ray indices for second ray.

     * @return A tensor representing the crossing point(s).
     */
    torch::Tensor ray_crossing(torch::Tensor view1, torch::Tensor ray1,
                               torch::Tensor view2, torch::Tensor ray2) const;

    /**
     * @brief Returns the pitch location measured in the given view (an index) of
     * the crossing point of ray grid coordinates one and two.
     *
     * @param view1 A scalar tensor or of shape (nbatch,) giving view indices for first ray.
     * @param ray1 A scalar tensor or of shape (nbatch,) giving ray indices for first ray.
     * @param view2 A scalar tensor or of shape (nbatch,) giving view indices for second ray.
     * @param ray2 A scalar tensor or of shape (nbatch,) giving ray indices for second ray.
     * @param view3 A scalar tensor or of shape (nboatch,) representing the third view index.
     * @return A tensor representing the pitch location.
     */
    torch::Tensor pitch_location(torch::Tensor view1, torch::Tensor ray1,
                                 torch::Tensor view2, torch::Tensor ray2,
                                 torch::Tensor view3) const;

    /**
     * @brief Returns the index of the closest ray at a location in the view that
     * is less than or equal to the given pitch.
     *
     * @param pitch A tensor representing the pitch value(s).
     * @param view A scalar tensor representing the view index.
     * @return A tensor of long integer values representing the pitch index.
     */
    torch::Tensor pitch_index(const torch::Tensor& pitch, const torch::Tensor& view) const;


    /**
     * @brief Returns the half-open indices bounds on each view that contain the
     * crossing points of the first two views.  Returns shape (nviews, 2)
     * giving.  Each row is (min,max).
     */
    torch::Tensor active_bounds() const;

private:
    /**
     * @brief Initializes or reinitializes the coordinate system.
     * @param pitches A tensor of shape (Nviews, 2, 2) representing the views.
     */
    void init(const torch::Tensor& pitches);
};

} // namespace RayGrid
} // namespace Spng
} // namespace WireCell

#endif // WIRECELL_SPNG_RAYGRID_RAYGRID_H

