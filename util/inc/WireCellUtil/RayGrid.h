/** A ray grid is defined by pair-wise layers of parallel, ordered pairs of rays.
 *
 * This file includes the implemenation for a ray grid as well as
 * classes and functions that operate with it.
 *
 * See section in the WCT manual for details.
 */

#ifndef WIRECELL_RAYGRID
#define WIRECELL_RAYGRID

#include "WireCellUtil/Point.h"
#include "WireCellUtil/ObjectArray2d.h"
#include "WireCellUtil/MultiArray.h"
#include <vector>
#include <map>

namespace WireCell {

    namespace RayGrid {

        /// A ray grid layer associates the set of parallel rays
        /// generated from the seeding pair.  All layers in a ray grid
        /// are considered co-planar.
        typedef int layer_index_t;

        /// Within one layer any ray may be identified by its index in
        /// a uniform linear grid.  Grid index 0 is the first ray of
        /// the defining pair, grid index 1 is the second.
        typedef int grid_index_t;

        /// A ray is located in the ray grid through its layer and its
        /// grid indices.
        struct coordinate_t {
            layer_index_t layer;
            grid_index_t grid;
        };

        /// Identify a pair of layers.
        typedef std::pair<layer_index_t, layer_index_t> layer_pair_t;

        /// Identify a triple of layers.
        typedef std::tuple<layer_index_t, layer_index_t, layer_index_t> layer_triple_t;

        /// A tensor type used to connect two ray coordinates and their
        /// crossing point location in a third layer.
        typedef boost::multi_array<double, 3> tensor_t;

        /// A 1D array of Vectors.
        typedef std::vector<Vector> vector_array1d_t;
        /// A 2D array'ish of Vectors.
        typedef ObjectArray2d<Vector> vector_array2d_t;

        /// A half open range between two grid indices.
        typedef std::pair<grid_index_t, grid_index_t> grid_range_t;

        /// A crossing is identified by the addresses of two rays.
        typedef std::pair<coordinate_t, coordinate_t> crossing_t;
        typedef std::vector<crossing_t> crossings_t;

        /// Bundle ray grid coordinate operations.
        class Coordinates {
           public:
            Coordinates();
            ~Coordinates();
            Coordinates(const Coordinates& other);
            Coordinates& operator=(const Coordinates& other);

            /// Create a ray grid by specifying the axis of projection.
            Coordinates(const ray_pair_vector_t& rays, int normal_axis = 0, double normal_location = 0.0);

            /// Return the crossing point of the index=0 rays for two
            /// layers.  As a special case, if one==two, the center
            /// point of ray 0 is returned.
            Vector zero_crossing(layer_index_t one, layer_index_t two) const;

            /// Return the crossing point of two rays.
            Vector ray_crossing(const coordinate_t& one, const coordinate_t& two) const;

            /// Return points of corners sorted by their angle about their center point.
            vector_array1d_t ring_points(const crossings_t& corners) const;

            /// Return the pitch location of the crossing point of two coordinate rays as measured in an other layer.
            double pitch_location(const coordinate_t& one, const coordinate_t& two, layer_index_t other) const;

            // fixme: this is in PR #187
            // Return points of corners sorted by their angle about their center point.
            // vector_array1d_t ring_points(const crossings_t& corners) const;

            // Return relative number of pitch steps from ray 0 to the pitch location.
            double pitch_relative(double pitch, layer_index_t layer) const {
                return pitch / m_pitch_mag[layer];
            }

            // Return the index of the ray bin covering the pitch location.
            int pitch_index(double pitch, layer_index_t layer) const 
            {
                return std::floor(pitch_relative(pitch, layer));
            }

            /// Number of layers
            int nlayers() const { return m_nlayers; }

            /// Return array of pitch magnitudes for the layers.
            const std::vector<double>& pitch_mags() const { return m_pitch_mag; }

            /// Return an array of vector.  The vector at layer index
            /// "l" gives a unit vector in the pitch direction of the
            /// layer.
            const vector_array1d_t& pitch_dirs() const { return m_pitch_dir; }

            /// Return an array of vector.  The vector at layer index
            /// "l" gives the center of ray 0 for that layer.
            const vector_array1d_t& centers() const { return m_center; }

            /// The ray_jumps collection provide a 2D array of
            /// vectors.  The vector at indices [a][b] is along the
            /// direction of rays in layer "a" has a magnitude which
            /// is the distance between two neighboring rays in layer
            /// "b" as measured along the direction or "a".
            ///
            /// As a special case, the elements [a][a] on the diagonal
            /// supply a unit vector in the direction of rays of layer
            /// "a".
            const vector_array2d_t& ray_jumps() const { return m_ray_jump; }

            /// Return the underlying ray grid "a" tensor.
            const tensor_t a() const { return m_a; }

            /// Return the underlying ray grid "b" tensor.
            const tensor_t b() const { return m_b; }

           private:
            int m_nlayers{0};

            // Pitch magnitude for each layer
            std::vector<double> m_pitch_mag;

            // The unit vector in the pitch direction for each layer
            vector_array1d_t m_pitch_dir;

            // A point (center point) on ray 0 of each layer
            vector_array1d_t m_center;

            // Zero-rays crossing points indexed by layer index pairs.
            // Symmetric array, diagonal is invalid.
            vector_array2d_t m_zero_crossing{};

            // Element (l,m) holds a relative vector which jumps along ray
            // direction of layer l between crossing points of neighboring
            // rays of layer m.  Not symmectric, and diagonal is invalid.
            vector_array2d_t m_ray_jump{};

            // Coefficients for fast pitch location calculation.  These
            // are scalar values indexed by three different layer
            // indicies.
            tensor_t m_a, m_b;
        };

        inline std::ostream& operator<<(std::ostream& os, const WireCell::RayGrid::coordinate_t& ra)
        {
            os << "<rayaddr {L" << ra.layer << ",G" << ra.grid << "}>";
            return os;
        }



    }  // namespace RayGrid
}  // namespace WireCell

#endif
