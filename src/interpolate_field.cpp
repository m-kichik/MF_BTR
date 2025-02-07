#include "MF.hpp"
#include "Interpolation.hpp"


/**
* Creates a new grid by scaling the old grid by the given coefficient.
*
* @param old_grid The old grid.
* @param scaling_coeff The scaling coefficient.
* @return The new grid.
*/
Grid make_new_grid(Grid old_grid, int &scaling_coeff) {
    auto size = old_grid.size * scaling_coeff - scaling_coeff + 1;
    auto min = old_grid.min;
    auto max = old_grid.max;
    auto step = old_grid.step / scaling_coeff;
    std::vector<double> grid(size);

    for (auto i = 0; i < size; i++) {
        if (i % scaling_coeff == 0) {
            grid[i] = old_grid.grid[i / scaling_coeff];
        } else {
            grid[i] = old_grid.grid[i / scaling_coeff] + step;
        }
    }

    return Grid{size, min, max, step, grid};
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

    for (auto i = 0; i < this->x_grid.size - 1; i ++) { // -1 because we process the last cell on the max_size - 1 step
        for (auto j = 0; j < this->y_grid.size - 1; j ++) {
            for (auto k = 0; k < this->z_grid.size - 1; k ++) {
                interpolate_in_cell({i, j, k}, new_field, new_x_grid, new_y_grid, new_z_grid);
            }
        }

    }
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
* @return The derivatives (Derivatives struct) of the field at the given point.
*/
Derivatives compute_derivatives(
    const uint64_t &x_idx, const uint64_t &y_idx, const uint64_t &z_idx,
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid,
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
    ) {
        Derivatives derivs;

        auto x_idx_prev = (x_idx == 0) ? 0 : x_idx - 1;
        auto y_idx_prev = (y_idx == 0) ? 0 : y_idx - 1;
        auto z_idx_prev = (z_idx == 0) ? 0 : z_idx - 1;

        auto x_idx_next = (x_idx == x_grid.size - 1) ? x_idx : x_idx + 1;
        auto y_idx_next = (y_idx == y_grid.size - 1) ? y_idx : y_idx + 1;
        auto z_idx_next = (z_idx == z_grid.size - 1) ? z_idx : z_idx + 1;

        for (auto i = 0; i < 3; i ++) {
            derivs.f[i] = field[x_idx][y_idx][z_idx][i];

            derivs.dfdx[i] = (field[x_idx_next][y_idx][z_idx][i] - field[x_idx_prev][y_idx][z_idx][i]) / (2 * x_grid.step);
            derivs.dfdy[i] = (field[x_idx][y_idx_next][z_idx][i] - field[x_idx][y_idx_prev][z_idx][i]) / (2 * y_grid.step);
            derivs.dfdz[i] = (field[x_idx][y_idx][z_idx_next][i] - field[x_idx][y_idx][z_idx_prev][i]) / (2 * z_grid.step);

            derivs.d2fdxdy[i] = (
                field[x_idx_next][y_idx_next][z_idx][i] - field[x_idx_next][y_idx][z_idx][i] 
                - field[x_idx][y_idx_next][z_idx][i] + field[x_idx_prev][y_idx_prev][z_idx][i]
                ) / (2 * x_grid.step * 2 * y_grid.step);

            derivs.d2fdxdz[i] = (
                field[x_idx_next][y_idx][z_idx_next][i] - field[x_idx_next][y_idx][z_idx][i] 
                - field[x_idx][y_idx][z_idx_next][i] + field[x_idx_prev][y_idx][z_idx_prev][i]
                ) / (2 * x_grid.step * 2 * z_grid.step);

            derivs.d2fdydz[i] = (
                field[x_idx][y_idx_next][z_idx_next][i] - field[x_idx][y_idx_next][z_idx][i] 
                - field[x_idx][y_idx][z_idx_next][i] + field[x_idx][y_idx_prev][z_idx_prev][i]
                ) / (2 * y_grid.step * 2 * z_grid.step);

            derivs.d3fdxdydz[i] = (
                field[x_idx_next][y_idx_next][z_idx_next][i] - field[x_idx_next][y_idx][z_idx_prev][i]
                - field[x_idx][y_idx_next][z_idx_prev][i] + field[x_idx_prev][y_idx_prev][z_idx_next][i]
                - field[x_idx_next][y_idx_next][z_idx_prev][i] + field[x_idx_next][y_idx][z_idx_prev][i]
                + field[x_idx][y_idx_next][z_idx_prev][i] - field[x_idx_prev][y_idx_prev][z_idx_prev][i]
                ) / (2 * x_grid.step * 2 * y_grid.step * 2 * z_grid.step);
        }

        return derivs;
}


void MagneticField::interpolate_in_cell(
    std::tuple<uint64_t, uint64_t, uint64_t> coords,
    field_type &field,
    Grid &x_grid, Grid &y_grid, Grid &z_grid
    ) {
    auto [i, j, k] = coords;
    auto derivs1 = compute_derivatives(i, j, k, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs2 = compute_derivatives(i + 1, j, k, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs3 = compute_derivatives(i, j + 1, k, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs4 = compute_derivatives(i + 1, j + 1, k, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs5 = compute_derivatives(i, j, k + 1, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs6 = compute_derivatives(i + 1, j, k + 1, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs7 = compute_derivatives(i, j + 1, k + 1, this->x_grid, this->y_grid, this->z_grid, field);
    auto derivs8 = compute_derivatives(i + 1, j + 1, k + 1, this->x_grid, this->y_grid, this->z_grid, field);



    std::cout << i << " " << j << " " << k << std::endl;
}