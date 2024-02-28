#include "mf_3d.h"

std::array<double, 3> MagneticField::get_field(const float& x, const float& y, const float& z) {
    auto x_ids = this->lower_id_is_same(this->x_grid, x);
    auto y_ids = this->lower_id_is_same(this->y_grid, y);
    auto z_ids = this->lower_id_is_same(this->z_grid, z);
    
    return compute_field(x, y, z, x_ids, y_ids, z_ids);
};

bool MagneticField::is_in_node(const float& point, const float& step, const float& epsilon = 1e-5) {
    return (point / step) - std::floor(point / step) <= epsilon;
}

std::tuple<int, bool> MagneticField::lower_id_is_same(const Grid& grid, const float& point) {
    if (point <= grid.min_value) {
        return std::tuple<int, bool> {0, true};
    }

    if (point >= grid.max_value) {
        return std::tuple<int, bool> {grid.size - 1, true};
    }

    if (is_in_node(point, grid.step)) {
        auto idx = std::round((point - grid.min_value) / grid.step);
        return std::tuple<int, bool> {idx, true};
    }

    auto idx = (int)(std::floor((point - grid.min_value) / grid.step));
    return std::tuple<int, bool> {idx, false};
}

/*
Which interpolation tree:
Upper: YES, lower: NO:
                              field[x][y][z]
                    z0 ?= z1 /
                             \
                              1d-interpolation with z
          y0 ?= y1 /
                   \
                              1d-interpolation with y
                    z0 ?= z1 /
                             \
                              2d-interpolation with y-z
x0 ?= x1 /
         \
                              1d-interpolation with x
                    z0 ?= z1 /
                             \
                              2d-interpolation with x-z
          y0 ?= y1 /
                   \
                              2d-interpolation with x-y
                    z0 ?= z1 /
                             \
                              3d-interpolation with x-y-z
*/

std::array<double, 3> MagneticField::compute_field(
    const float& x, const float& y, const float& z,
    const std::tuple<int, bool>& x_idxs,    // idxs is a tuple with lower bound an bool
    const std::tuple<int, bool>& y_idxs,    // (true if upper bound is same, false else)
    const std::tuple<int, bool>& z_idxs
) {
    if (!std::get<1>(x_idxs)) {
        if (!std::get<1>(y_idxs)) {
            if (!std::get<1>(z_idxs)) {
                // 3d with x-y-z
                return interp_3d(
                    x, y, z,
                    this->x_grid.grid[std::get<0>(x_idxs)],
                    this->y_grid.grid[std::get<0>(y_idxs)],
                    this->z_grid.grid[std::get<0>(z_idxs)],
                    this->x_grid.step, this->y_grid.step, this->z_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs) + 1],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs) + 1],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs)][std::get<0>(z_idxs) + 1],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs) + 1]
                );
            } else {
                // 2d with x-y
                return interp_2d(
                    x, y,
                    this->x_grid.grid[std::get<0>(x_idxs)], this->y_grid.grid[std::get<0>(y_idxs)],
                    this->x_grid.step, this->y_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs)]
                );
            }
        } else {
            if (!std::get<1>(z_idxs)) {
                // 2d with x-z
                return interp_2d(
                    x, z,
                    this->x_grid.grid[std::get<0>(x_idxs)], this->z_grid.grid[std::get<0>(z_idxs)],
                    this->x_grid.step, this->z_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs) + 1],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs)][std::get<0>(z_idxs) + 1]
                );
            } else {
                // 1d with x
                return interp_1d(
                    x, this->x_grid.grid[std::get<0>(x_idxs)], this->x_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs) + 1][std::get<0>(y_idxs)][std::get<0>(z_idxs)]
                );
            }
        } 
    } else {
        if (!std::get<1>(y_idxs)) {
            if (!std::get<1>(z_idxs)) {
                // 2d with y-z
                return interp_2d(
                    y, z,
                    this->y_grid.grid[std::get<0>(y_idxs)], this->z_grid.grid[std::get<0>(z_idxs)],
                    this->y_grid.step, this->z_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs) + 1],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs) + 1]
                );
            } else {
                // 1d with y
                return interp_1d(
                    y, this->y_grid.grid[std::get<0>(y_idxs)], this->y_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs) + 1][std::get<0>(z_idxs)]
                );
            }
        } else {
            if (!std::get<1>(z_idxs)) {
                // 1d with z
                return interp_1d(
                    z, this->z_grid.grid[std::get<0>(z_idxs)], this->z_grid.step,
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)],
                    this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs) + 1]
                );
            } else {
                return this->field[std::get<0>(x_idxs)][std::get<0>(y_idxs)][std::get<0>(z_idxs)];
            }
        }
    }
}
