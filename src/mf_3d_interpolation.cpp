#include "mf_3d.h"

bool is_in_node(const float& point, const float& step, const float& epsilon = 1e-5) {
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
    return std::tuple<int, int> {idx, false};
}

Grid make_new_grid(const Grid &old_grid, const std::array<int, 3> &decomp_coeffs, const int &idx) {
    uint new_size = (old_grid.size - 1) * decomp_coeffs[idx] + 1;
    float new_step = old_grid.step / decomp_coeffs[idx];
    std::vector<float> new_grid;
    for (auto i = 0; i < new_size; ++i) {
        auto iidx = i / decomp_coeffs[idx];
        auto offset = (i % decomp_coeffs[idx]) * new_step;
        new_grid.push_back(old_grid.grid[iidx] + offset);
    }
    return Grid{
        new_size,
        new_step,
        old_grid.min_value,
        old_grid.max_value,
        new_grid
    };
}

/*
Interpolation tree:
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

/*
Interpolation tree:
Upper: YES, lower: NO:
                                    field[node]
                        z in node? /
                                   \
                                    cubic interpolation with z
            y in node? /
                       \
                                    cubic interpolation with y
                        z in node? /
                                   \
                                    bicubic interpolation with y-z
x in node? /
           \
                                    cubic interpolation with x
                        z in node? /
                                   \
                                    bicubic interpolation with x-z
            y in node? /
                       \
                                    bicubic interpolation with x-y
                        z in node? /
                                   \
                                   my_function interpolation with x-y-z
*/

void MagneticField::compute_interpolated_matrix() {
    std::vector<std::vector<std::vector<std::array<double, 3>>>> extended_field;

    Grid new_x_grid = make_new_grid(this->x_grid, this->decomp_coeffs, 0);
    Grid new_y_grid = make_new_grid(this->y_grid, this->decomp_coeffs, 1);
    Grid new_z_grid = make_new_grid(this->z_grid, this->decomp_coeffs, 2);

    extended_field.resize(
        new_x_grid.size,
        std::vector<std::vector<std::array<double, 3>>> (
            new_y_grid.size,
            std::vector<std::array<double, 3>>(
                new_z_grid.size,
                std::array<double, 3>{1e20, 1e20, 1e20})));

    for (auto i = 0; i < new_x_grid.size; i++) {
        for (auto j = 0; j < new_y_grid.size; j++) {
            for (auto k = 0; k < new_z_grid.size; k++) {
                if (i % decomp_coeffs[0] == 0) {
                    if (j % decomp_coeffs[1] == 0) {
                        if (k % decomp_coeffs[2] == 0) {
                            // field in node
                            extended_field[i][j][k] = this->field[i / decomp_coeffs[0]][j / decomp_coeffs[1]][k / decomp_coeffs[2]];
                        } else {
                            // cubic interpolation with z
                            extended_field[i][j][k] = cubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                2,
                                this->field
                            );
                        }
                    } else {
                        if (k % decomp_coeffs[2] == 0) {
                            // cubic interpolation with y
                            extended_field[i][j][k] = cubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                1,
                                this->field
                            );
                        } else {
                            // bicubic interpolation with y-z
                            extended_field[i][j][k] = bicubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                {0, 1, 1},
                                this->field
                            );
                        }
                    }
                } else {
                    if (j % decomp_coeffs[1] == 0) {
                        if (k % decomp_coeffs[2] == 0) {
                            // cubic interpolation with x
                            extended_field[i][j][k] = cubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                0,
                                this->field
                            );
                        } else {
                            // bicubic interpolation with x-z
                            extended_field[i][j][k] = bicubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                {1, 0, 1},
                                this->field
                            );
                        }
                    } else {
                        if (k % decomp_coeffs[2] == 0) {
                            // bicubic interpolation with x-y
                            extended_field[i][j][k] = bicubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                {1, 1, 0},
                                this->field
                            );
                        } else {
                            // tricubic interpolation with x-y-z
                            extended_field[i][j][k] = tricubic(
                                new_x_grid.grid[i], new_y_grid.grid[j], new_z_grid.grid[k],
                                i, j, k,
                                this->x_grid, this->y_grid, this->z_grid,
                                this->decomp_coeffs,
                                this->field
                            );
                        }
                    }
                }
            }
        }
    }

    this->x_grid = new_x_grid;
    this->y_grid = new_y_grid;
    this->z_grid = new_z_grid;

    this->field = extended_field;
}
