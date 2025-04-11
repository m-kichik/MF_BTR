#include <cmath>

#include "mf/MF.hpp"
#include "interpolation/Interpolation.hpp"

/**
 * @param x - x coordinate in absolute values
 * @param grid - corresponding grid
 * @return x_idx_left, x_idx_right, x_delta - left and right boundaries 
 * and delta between x and left boundary,
 * scaled between 0 and 1 (for interpolation)
 */
std::tuple<uint, uint, double> find_nodes(const double &x, const Grid &grid) {
    if (x <= grid.min) {
        return std::make_tuple(0u, 0u, 0.);
    }

    if (x >= grid.max) {
        return std::make_tuple(grid.size - 1, grid.size - 1, 0.);
    }

    auto idx = (uint)((x - grid.min) / grid.step);
    if (std::abs(x / grid.step - std::round(x / grid.step)) < 1e-12) {
        if (x - grid.grid[idx] > 1e-12) {
            idx += 1;
        }
        if (x - grid.grid[idx] < -1e-12) {
            idx -= 1;
        }
        return std::make_tuple(idx, idx, 0.);
    } else {
        return std::make_tuple(idx, idx + 1, (x - grid.grid[idx]) / grid.step);
    }
}

/**
 * Returns the magnetic field at the given coordinates.
 *
 * @param coords The coordinates of the point.
 * @return The magnetic field at the given coordinates.
 */
std::array<double, 3> MagneticField::get_field(const std::tuple<double, double, double> &coords) {
    auto [x, y, z] = coords;

    auto [x_left, x_right, x_delta] = find_nodes(x, this->x_grid_);
    auto [y_left, y_right, y_delta] = find_nodes(y, this->y_grid_);
    auto [z_left, z_right, z_delta] = find_nodes(z, this->z_grid_);

    if (x_left == x_right) {
        if (y_left == y_right) {
            if (z_left == z_right) {
                return this->field_[x_left][y_left][z_left];
            } else {
                return linear_1d(
                    this->field_[x_left][y_left][z_left], 
                    this->field_[x_left][y_left][z_right], 
                    z_delta
                    );
            }
        } else {
            if (z_left == z_right) {
                return linear_1d(
                    this->field_[x_left][y_left][z_left], 
                    this->field_[x_left][y_right][z_left], 
                    y_delta
                    );
            } else {
                return linear_2d(
                    this->field_[x_left][y_left][z_left],
                    this->field_[x_left][y_left][z_right],
                    this->field_[x_left][y_right][z_left],
                    this->field_[x_left][y_right][z_right],
                    y_delta,
                    z_delta
                );
            }
        }
    } else {
        if (y_left == y_right) {
            if (z_left == z_right) {
                return linear_1d(
                    this->field_[x_left][y_left][z_left], 
                    this->field_[x_right][y_left][z_left], 
                    x_delta
                    );
            } else {
                return linear_2d(
                    this->field_[x_left][y_left][z_left],
                    this->field_[x_left][y_left][z_right],
                    this->field_[x_right][y_left][z_left],
                    this->field_[x_right][y_left][z_right],
                    x_delta,
                    z_delta
                );
            }
        } else {
            if (z_left == z_right) {
                return linear_2d(
                    this->field_[x_left][y_left][z_left],
                    this->field_[x_left][y_right][z_left],
                    this->field_[x_right][y_left][z_left],
                    this->field_[x_right][y_right][z_left],
                    x_delta,
                    y_delta
                );
            } else {
                return linear_3d(
                    this->field_[x_left][y_left][z_left],
                    this->field_[x_left][y_left][z_right],
                    this->field_[x_left][y_right][z_left],
                    this->field_[x_left][y_right][z_right],
                    this->field_[x_right][y_left][z_left],
                    this->field_[x_right][y_left][z_right],
                    this->field_[x_right][y_right][z_left],
                    this->field_[x_right][y_right][z_right],
                    x_delta,
                    y_delta,
                    z_delta
                );
            }

        }
    }
}
