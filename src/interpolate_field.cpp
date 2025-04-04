#include "MF.hpp"
#include "Interpolation.hpp"


/**
* Creates a new grid by scaling the old grid by the given coefficient.
*
* @param old_grid The old grid.
* @param scaling_coeff The scaling coefficient.
* @return The new grid.
*/
// Grid make_new_grid(Grid old_grid, int &scaling_coeff) {
//     auto size = (old_grid.size - 1) * scaling_coeff + 1;
//     auto min = old_grid.min;
//     auto max = old_grid.max;
//     auto step = old_grid.step / scaling_coeff;
//     std::vector<double> grid(size);

//     for (auto i = 0; i < size; i++) {
//         if (i % scaling_coeff == 0) {
//             grid[i] = old_grid.grid[i / scaling_coeff];
//         } else {
//             grid[i] = old_grid.grid[i / scaling_coeff] + step;
//         }
//     }

//     return Grid{size, min, max, step, grid};
// }
Grid make_new_grid(Grid old_grid, int scaling_coeff) {
    auto new_size = (old_grid.size - 1) * scaling_coeff + 1;
    std::vector<double> grid(new_size);

    for (int i = 0; i < new_size; i++) {
        int old_idx = i / scaling_coeff;
        double t = static_cast<double>(i % scaling_coeff) / scaling_coeff;
        grid[i] = old_grid.grid[old_idx] + t * (old_grid.grid[old_idx + 1] - old_grid.grid[old_idx]);
    }

    double new_step = (old_grid.grid[1] - old_grid.grid[0]) / scaling_coeff;
    return Grid{new_size, old_grid.min, old_grid.max, new_step, grid};
}


/**
* Interpolates the field. 
* Creates new grid based on the provided coefficients.
* Creates new resized field. Interpolates field sequentially in each cell of the new grid.
*/
void MagneticField::interpolate_field() {
    auto new_x_grid = make_new_grid(this->x_grid, this->scaling_coeffs[0]);
    auto new_y_grid = make_new_grid(this->y_grid, this->scaling_coeffs[1]);
    auto new_z_grid = make_new_grid(this->z_grid, this->scaling_coeffs[2]);

    MagneticField::field_type new_field;
    new_field.resize(
        new_x_grid.size, 
        std::vector<std::vector<std::array<double, 3>>> (
            new_y_grid.size, 
            std::vector<std::array<double, 3>> (
                new_z_grid.size, std::array<double, 3> {
                    std::numeric_limits<double>::infinity(), 
                    std::numeric_limits<double>::infinity(), 
                    std::numeric_limits<double>::infinity()
                    }
                )
            )
        );

    // std::cout << new_x_grid.size << " " << new_y_grid.size << " " << new_z_grid.size << std::endl;

    for (auto i = 0; i < this->x_grid.size - 1; i ++) { // -1 because we process the last cell on the max_size - 1 step
        for (auto j = 0; j < this->y_grid.size - 1; j ++) {
            for (auto k = 0; k < this->z_grid.size - 1; k ++) {
                // here we interpolate all values in cell (i, j, k) - (i + 1, j + 1, k + 1)
                interpolate_in_cell({i, j, k}, new_field, new_x_grid, new_y_grid, new_z_grid);
            }
        }
    };

    this->x_grid = new_x_grid;
    this->y_grid = new_y_grid;
    this->z_grid = new_z_grid;
    this->field = new_field;
}


/**
* Computes the derivatives of the field at the given point.
*
* @param x_idx The x index of the point.
* @param y_idx The y index of the point.
* @param z_idx The z index of the point.
* @param x_grid The x grid.
* @param y_grid The y grid.
* @param z_grid The z grid.
* @param field The field.
* @return The derivatives of the field at the given point: 
* array {f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz}.
*/
std::array<std::array<double, 3>, 8> compute_derivatives(
    const uint64_t &x_idx, const uint64_t &y_idx, const uint64_t &z_idx,
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid,
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
    ) {
        std::array<std::array<double, 3>, 8> derivatives;

        // auto x_idx_prev = (x_idx == 0) ? 0 : x_idx - 1;
        // auto y_idx_prev = (y_idx == 0) ? 0 : y_idx - 1;
        // auto z_idx_prev = (z_idx == 0) ? 0 : z_idx - 1;

        // auto x_idx_next = (x_idx == x_grid.size - 1) ? x_idx : x_idx + 1;
        // auto y_idx_next = (y_idx == y_grid.size - 1) ? y_idx : y_idx + 1;
        // auto z_idx_next = (z_idx == z_grid.size - 1) ? z_idx : z_idx + 1;

        // Boundary handling
        auto x_idx_prev = (x_idx == 0) ? x_idx + 1 : x_idx - 1;  // Forward difference at left boundary
        auto x_idx_next = (x_idx == x_grid.size - 1) ? x_idx - 1 : x_idx + 1;  // Backward difference at right boundary
        auto y_idx_prev = (y_idx == 0) ? y_idx + 1 : y_idx - 1;
        auto y_idx_next = (y_idx == y_grid.size - 1) ? y_idx - 1 : y_idx + 1;
        auto z_idx_prev = (z_idx == 0) ? z_idx + 1 : z_idx - 1;
        auto z_idx_next = (z_idx == z_grid.size - 1) ? z_idx - 1 : z_idx + 1;

        for (auto i = 0; i < 3; i ++) {
            // calculate f
            derivatives[0][i] = field[x_idx][y_idx][z_idx][i];

            // calculate dfdx, dfdy, dfdz
            derivatives[1][i] = (field[x_idx_next][y_idx][z_idx][i] - field[x_idx_prev][y_idx][z_idx][i]) / (2 * x_grid.step);
            derivatives[2][i] = (field[x_idx][y_idx_next][z_idx][i] - field[x_idx][y_idx_prev][z_idx][i]) / (2 * y_grid.step);
            derivatives[3][i] = (field[x_idx][y_idx][z_idx_next][i] - field[x_idx][y_idx][z_idx_prev][i]) / (2 * z_grid.step);

            // calculate d2fdxdy, d2fdxdz, d2fdydz
            derivatives[4][i] = (
                field[x_idx_next][y_idx_next][z_idx][i] - field[x_idx_next][y_idx_prev][z_idx][i]
                - field[x_idx_prev][y_idx_next][z_idx][i] + field[x_idx_prev][y_idx_prev][z_idx][i]
                ) / (4 * x_grid.step * y_grid.step);

            derivatives[5][i] = (
                field[x_idx_next][y_idx][z_idx_next][i] - field[x_idx_next][y_idx][z_idx][i] 
                - field[x_idx][y_idx][z_idx_next][i] + field[x_idx_prev][y_idx][z_idx_prev][i]
                ) / (2 * x_grid.step * 2 * z_grid.step);

            derivatives[6][i] = (
                field[x_idx][y_idx_next][z_idx_next][i] - field[x_idx][y_idx_next][z_idx][i] 
                - field[x_idx][y_idx][z_idx_next][i] + field[x_idx][y_idx_prev][z_idx_prev][i]
                ) / (2 * y_grid.step * 2 * z_grid.step);
            
            // calculate d3fdxdydz
            // derivatives[7][i] = (
            //     field[x_idx_next][y_idx_next][z_idx_next][i] - field[x_idx_next][y_idx][z_idx_prev][i]
            //     - field[x_idx][y_idx_next][z_idx_prev][i] + field[x_idx_prev][y_idx_prev][z_idx_next][i]
            //     - field[x_idx_next][y_idx_next][z_idx_prev][i] + field[x_idx_next][y_idx][z_idx_prev][i]
            //     + field[x_idx][y_idx_next][z_idx_prev][i] - field[x_idx_prev][y_idx_prev][z_idx_prev][i]
            //     ) / (2 * x_grid.step * 2 * y_grid.step * 2 * z_grid.step);

            derivatives[7][i] = (
                field[x_idx_next][y_idx_next][z_idx_next][i] - field[x_idx_next][y_idx_next][z_idx_prev][i]
                - field[x_idx_next][y_idx_prev][z_idx_next][i] + field[x_idx_next][y_idx_prev][z_idx_prev][i]
                - field[x_idx_prev][y_idx_next][z_idx_next][i] + field[x_idx_prev][y_idx_next][z_idx_prev][i]
                + field[x_idx_prev][y_idx_prev][z_idx_next][i] - field[x_idx_prev][y_idx_prev][z_idx_prev][i]
                ) / (8 * x_grid.step * y_grid.step * z_grid.step);

            // Scale derivatives for unit cube
            derivatives[1][i] *= x_grid.step;
            derivatives[2][i] *= y_grid.step;
            derivatives[3][i] *= z_grid.step;
            derivatives[4][i] *= x_grid.step * y_grid.step;
            derivatives[5][i] *= x_grid.step * z_grid.step;
            derivatives[6][i] *= y_grid.step * z_grid.step;
            derivatives[7][i] *= x_grid.step * y_grid.step * z_grid.step;
        }

        return derivatives;
}

// std::array<std::array<double, 3>, 8> compute_derivatives(
//     const uint64_t &x_idx, const uint64_t &y_idx, const uint64_t &z_idx,
//     const Grid &x_grid, const Grid &y_grid, const Grid &z_grid,
//     const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
//     ) {
//         std::array<std::array<double, 3>, 8> derivatives;

//         // Corrected boundary handling for x-direction
//         uint64_t x_idx_prev, x_idx_next;
//         if (x_idx == 0) {
//             x_idx_prev = x_idx;
//             x_idx_next = x_idx + 1;
//         } else if (x_idx == x_grid.size - 1) {
//             x_idx_prev = x_idx - 1;
//             x_idx_next = x_idx;
//         } else {
//             x_idx_prev = x_idx - 1;
//             x_idx_next = x_idx + 1;
//         }

//         // Corrected boundary handling for y-direction
//         uint64_t y_idx_prev, y_idx_next;
//         if (y_idx == 0) {
//             y_idx_prev = y_idx;
//             y_idx_next = y_idx + 1;
//         } else if (y_idx == y_grid.size - 1) {
//             y_idx_prev = y_idx - 1;
//             y_idx_next = y_idx;
//         } else {
//             y_idx_prev = y_idx - 1;
//             y_idx_next = y_idx + 1;
//         }

//         // Corrected boundary handling for z-direction
//         uint64_t z_idx_prev, z_idx_next;
//         if (z_idx == 0) {
//             z_idx_prev = z_idx;
//             z_idx_next = z_idx + 1;
//         } else if (z_idx == z_grid.size - 1) {
//             z_idx_prev = z_idx - 1;
//             z_idx_next = z_idx;
//         } else {
//             z_idx_prev = z_idx - 1;
//             z_idx_next = z_idx + 1;
//         }

//         for (auto i = 0; i < 3; i ++) {
//             // Calculate f
//             derivatives[0][i] = field[x_idx][y_idx][z_idx][i];

//             // Calculate dfdx using forward/backward/central differences
//             if (x_idx == 0 || x_idx == x_grid.size - 1) {
//                 derivatives[1][i] = (field[x_idx_next][y_idx][z_idx][i] - field[x_idx_prev][y_idx][z_idx][i]) / (x_grid.step);
//             } else {
//                 derivatives[1][i] = (field[x_idx_next][y_idx][z_idx][i] - field[x_idx_prev][y_idx][z_idx][i]) / (2 * x_grid.step);
//             }

//             // Similarly for dfdy and dfdz
//             if (y_idx == 0 || y_idx == y_grid.size - 1) {
//                 derivatives[2][i] = (field[x_idx][y_idx_next][z_idx][i] - field[x_idx][y_idx_prev][z_idx][i]) / (y_grid.step);
//             } else {
//                 derivatives[2][i] = (field[x_idx][y_idx_next][z_idx][i] - field[x_idx][y_idx_prev][z_idx][i]) / (2 * y_grid.step);
//             }

//             if (z_idx == 0 || z_idx == z_grid.size - 1) {
//                 derivatives[3][i] = (field[x_idx][y_idx][z_idx_next][i] - field[x_idx][y_idx][z_idx_prev][i]) / (z_grid.step);
//             } else {
//                 derivatives[3][i] = (field[x_idx][y_idx][z_idx_next][i] - field[x_idx][y_idx][z_idx_prev][i]) / (2 * z_grid.step);
//             }

//             // Calculate d2fdxdy (ensure proper handling)
//             derivatives[4][i] = (
//                 field[x_idx_next][y_idx_next][z_idx][i] - field[x_idx_next][y_idx_prev][z_idx][i]
//                 - field[x_idx_prev][y_idx_next][z_idx][i] + field[x_idx_prev][y_idx_prev][z_idx][i]
//             ) / (4 * x_grid.step * y_grid.step);

//             // Similar corrections for d2fdxdz and d2fdydz
//             derivatives[5][i] = (
//                 field[x_idx_next][y_idx][z_idx_next][i] - field[x_idx_next][y_idx][z_idx_prev][i]
//                 - field[x_idx_prev][y_idx][z_idx_next][i] + field[x_idx_prev][y_idx][z_idx_prev][i]
//             ) / (4 * x_grid.step * z_grid.step);

//             derivatives[6][i] = (
//                 field[x_idx][y_idx_next][z_idx_next][i] - field[x_idx][y_idx_prev][z_idx_next][i]
//                 - field[x_idx][y_idx_next][z_idx_prev][i] + field[x_idx][y_idx_prev][z_idx_prev][i]
//             ) / (4 * y_grid.step * z_grid.step);

//             // Correct d3fdxdydz calculation
//             derivatives[7][i] = (
//                 field[x_idx_next][y_idx_next][z_idx_next][i] - field[x_idx_next][y_idx_next][z_idx_prev][i]
//                 - field[x_idx_next][y_idx_prev][z_idx_next][i] + field[x_idx_next][y_idx_prev][z_idx_prev][i]
//                 - field[x_idx_prev][y_idx_next][z_idx_next][i] + field[x_idx_prev][y_idx_next][z_idx_prev][i]
//                 + field[x_idx_prev][y_idx_prev][z_idx_next][i] - field[x_idx_prev][y_idx_prev][z_idx_prev][i]
//             ) / (8 * x_grid.step * y_grid.step * z_grid.step);

//             // Scale derivatives for unit cube
//             // derivatives[1][i] *= x_grid.step;
//             // derivatives[2][i] *= y_grid.step;
//             // derivatives[3][i] *= z_grid.step;
//             // derivatives[4][i] *= x_grid.step * y_grid.step;
//             // derivatives[5][i] *= x_grid.step * z_grid.step;
//             // derivatives[6][i] *= y_grid.step * z_grid.step;
//             // derivatives[7][i] *= x_grid.step * y_grid.step * z_grid.step;
//         }

//         return derivatives;
// }


void MagneticField::interpolate_in_cell(
    std::tuple<uint64_t, uint64_t, uint64_t> coords,
    field_type &new_field,
    Grid &new_x_grid, Grid &new_y_grid, Grid &new_z_grid
    ) {
    auto [cell_i, cell_j, cell_k] = coords;

    // calculate field derivatives in cube
    //
    //    Z
    //    |   p7---------p8
    //    p5---------p6  |
    //    |   |      |   |
    //    |   |   Y  |   |
    //    |   p3-----|---p4
    //    p1---------p2 ------ X
    //
    // each vector has 8 components: f, df/dx, df/dy, df/dz, d2f/dxdy, d2f/dxdz, d2f/dydz, d3f/dxdydz
    // for full scheme see F. Lekien et. al. "Tricubic interpolation in three dimensions"

    auto derivs1 = compute_derivatives(cell_i, cell_j, cell_k, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs2 = compute_derivatives(cell_i + 1, cell_j, cell_k, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs3 = compute_derivatives(cell_i, cell_j + 1, cell_k, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs4 = compute_derivatives(cell_i + 1, cell_j + 1, cell_k, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs5 = compute_derivatives(cell_i, cell_j, cell_k + 1, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs6 = compute_derivatives(cell_i + 1, cell_j, cell_k + 1, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs7 = compute_derivatives(cell_i, cell_j + 1, cell_k + 1, this->x_grid, this->y_grid, this->z_grid, this->field);
    auto derivs8 = compute_derivatives(cell_i + 1, cell_j + 1, cell_k + 1, this->x_grid, this->y_grid, this->z_grid, this->field);
    
    // Solve tricubic interpolation in 3D cell volume
    // fill vector b
    std::array<std::array<double, 3>, 64> b_3d;
    for (auto i = 0; i < 64; i += 8) {
        b_3d[i] = derivs1[i / 8];
        b_3d[i + 1] = derivs2[i / 8];
        b_3d[i + 2] = derivs3[i / 8];
        b_3d[i + 3] = derivs4[i / 8];
        b_3d[i + 4] = derivs5[i / 8];
        b_3d[i + 5] = derivs6[i / 8];
        b_3d[i + 6] = derivs7[i / 8];
        b_3d[i + 7] = derivs8[i / 8];
    }

    // solve alpha = B^-1 * b from the paper "Tricubic interpolation in three dimensions"
    std::array<std::array<double, 3>, 64> alpha_3d;
    for (auto i = 0; i < 64; i++) {
        alpha_3d[i] = {0, 0, 0};
        for (auto j = 0; j < 64; j++) {
            for (auto fid = 0; fid < 3; fid ++) {
                alpha_3d[i][fid] += tricubic_B[i][j] * b_3d[j][fid];
            } 
        }
    }

    // iterate through new nodes in the cell and interpolate the field
    for (auto node_i = 0; node_i <= this->scaling_coeffs[0]; node_i++) {
        auto x_tail = double(node_i) / this->scaling_coeffs[0]; // works only for regular grid
        for (auto node_j = 0; node_j <= this->scaling_coeffs[1]; node_j++) {
            auto y_tail = double(node_j) / this->scaling_coeffs[1];
            for (auto node_k = 0; node_k <= this->scaling_coeffs[2]; node_k++) {
                auto z_tail = double(node_k) / this->scaling_coeffs[2];
                // solve basic equation from the paper
                std::array<double, 3> f {0., 0., 0.};
                for (auto i = 0; i <= 3; i++) {
                    for (auto j = 0; j <= 3; j++) {
                        for (auto k = 0; k <= 3; k++) {
                            // finally iterate through components of the field
                            for (auto cidx = 0; cidx < 3; cidx++) {
                                f[cidx] += alpha_3d[i + 4 * j + 16 * k][cidx] * (
                                    std::pow(x_tail, i) * std::pow(y_tail, j) * std::pow(z_tail, k)
                                );
                            }
                        }
                    }
                }
                new_field[
                    node_i + cell_i * this->scaling_coeffs[0]
                ][
                    node_j + cell_j * this->scaling_coeffs[1]
                ][
                    node_k + cell_k * this->scaling_coeffs[2]
                ] = f;
                // std::cout << node_i + cell_i * this->scaling_coeffs[0] 
                //         << " " << node_j + cell_j * this->scaling_coeffs[1] 
                //         << " " << node_k + cell_k * this->scaling_coeffs[2] 
                //         << std::endl;
            }
        }
    }

    // fill node corners with original field (not necessary, we didn't see any corruption of the original values)
    for (auto i = cell_i; i <= cell_i + 1; i++) {
        for (auto j = cell_j; j <= cell_j + 1; j++) {
            for (auto k = cell_k; k <= cell_k + 1; k++) {
                new_field[
                    i * this->scaling_coeffs[0]
                ][
                    j * this->scaling_coeffs[1]
                ][
                    k * this->scaling_coeffs[2]
                ] = this->field[i][j][k];
            }
        }
    }
}