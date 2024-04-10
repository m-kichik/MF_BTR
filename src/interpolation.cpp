#include "mf_3d.h"

std::array<double, 3> interp_1d(
    const float& x, const float& x0, const float& step,
    const std::array<double, 3>& f0, const std::array<double, 3>& f1
) {
    auto k = (x - x0) / step;

    std::array<double, 3> f;

    for (auto i = 0; i < 3; i++) {
        f[i] = f0[i] + (f1[i] - f0[i]) * k;
    }

    return f;
}

std::array<double, 3> interp_2d(
    const float& x, const float& y, const float& x0, const float& y0,
    const float& step_x, const float& step_y,
    const std::array<double, 3>& f00, const std::array<double, 3>& f01,
    const std::array<double, 3>& f10, const std::array<double, 3>& f11
) {
    auto x1 = x0 + step_x;
    auto y1 = y0 + step_y;
    auto k = 1 / (step_x * step_y);

    std::array<double, 3> f;

    for (auto i = 0; i < 3; i++) {
        f[i] = k * (f00[i] * (x1 - x) * (y1 - y) + 
                    f10[i] * (x - x0) * (y1 - y) +
                    f01[i] * (x1 - x) * (y - y0) +
                    f11[i] * (x - x0) * (y - y0));
    }

    return f;
}

std::array<double, 3> interp_3d(
    const float& x, const float& y, const float& z,
    const float& x0, const float& y0, const float& z0,
    const float& step_x, const float& step_y, const float& step_z,
    const std::array<double, 3>& f000, const std::array<double, 3>& f001,
    const std::array<double, 3>& f010, const std::array<double, 3>& f011,
    const std::array<double, 3>& f100, const std::array<double, 3>& f101,
    const std::array<double, 3>& f110, const std::array<double, 3>& f111
) {
    auto x1 = x0 + step_x;
    auto y1 = y0 + step_y;
    auto z1 = z0 + step_z;
    auto k = 1 / (step_x * step_y * step_z);

    std::array<double, 3> f;

    for (auto i = 0; i < 3; i++) {
        f[i] = k * (f000[i] * (x1 - x) *  (y1 - y) * (z1 - z) +
                    f001[i] * (x1 - x) *  (y1 - y) * (z - z0) +
                    f010[i] * (x1 - x) *  (y - y0) * (z1 - z) +
                    f011[i] * (x1 - x) *  (y - y0) * (z - z0) +
                    f100[i] * (x - x0) *  (y1 - y) * (z1 - z) +
                    f101[i] * (x - x0) *  (y1 - y) * (z - z0) +
                    f110[i] * (x - x0) *  (y - y0) * (z1 - z) +
                    f111[i] * (x - x0) *  (y - y0) * (z - z0));
    }

    return f;
}

Derivatives compute_derivatives(
    const uint &x_idx, const uint &y_idx, const uint &z_idx,
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
            derivs.dfdx[i] = (field[x_idx][y_idx][z_idx_next][i] - field[x_idx][y_idx][z_idx_prev][i]) / (2 * z_grid.step);

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

            derivs.d3fdxdyxz[i] = (
                field[x_idx_next][y_idx_next][z_idx_next][i] - field[x_idx_next][y_idx][z_idx_prev][i]
                - field[x_idx][y_idx_next][z_idx_prev][i] + field[x_idx_prev][y_idx_prev][z_idx_next][i]
                - field[x_idx_next][y_idx_next][z_idx_prev][i] + field[x_idx_next][y_idx][z_idx_prev][i]
                + field[x_idx][y_idx_next][z_idx_prev][i] - field[x_idx_prev][y_idx_prev][z_idx_prev][i]
                ) / (2 * x_grid.step * 2 * y_grid.step * 2 * z_grid.step);
        }

        return derivs;
}

std::array<double, 3> cubic(
    const double &x, const double &y, const double &z, // coords in new grid
    const uint &x_idx, const uint &y_idx, const uint &z_idx, // ids in new grid
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid, // old grids
    const std::array<int, 3> &dec_coeffs,
    const uint &mask1d, // 0 for x cubic, 1 for y, 2 for z
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
) {
    std::array<std::array<double, 3>, 4> b;
    std::array<std::array<double, 3>, 4> alpha;
    std::array<double, 3> f {0., 0., 0.};

    for (auto i = 0; i < 2; i++) {
        auto px = x_idx / dec_coeffs[0];
        auto py = y_idx / dec_coeffs[1];
        auto pz = z_idx / dec_coeffs[2];

        if (mask1d == 0) { px += i; }
        if (mask1d == 1) { py += i; }
        if (mask1d == 2) { pz += i; }

        auto derivs = compute_derivatives(px, py, pz, x_grid, y_grid, z_grid, field);
        b[2 * i] = derivs.f;
        if (mask1d == 0) { b[2 * i + 1] = derivs.dfdx; }
        if (mask1d == 1) { b[2 * i + 1] = derivs.dfdy; }
        if (mask1d == 2) { b[2 * i + 1] = derivs.dfdz; }
    }

    // calculate values in vector alpha

    for (auto i = 0; i < 4; i++) {
        alpha[i] = {0, 0, 0};
        for (auto j = 0; j < 4; j++) {
            for (auto fid = 0; fid < 3; fid ++) {
                alpha[i][fid] += cubic_B[i][j] * b[j][fid];
            } 
        }
    }

    // finally fill f values
    auto x_tail = std::fmod(x, x_grid.step);
    auto y_tail = std::fmod(y, y_grid.step);
    auto z_tail = std::fmod(z, z_grid.step);

    for (auto i = 0; i <= 3; i++) {
        for (auto fid = 0; fid < 3; fid++) {
            auto i_ = 0;
            auto j_ = 0;
            auto k_ = 0;
            if (mask1d == 0) { i_ = i; j_ = 0; k_ = 0; }
            if (mask1d == 1) { i_ = 0; j_ = i; k_ = 0; }
            if (mask1d == 2) { i_ = 0; j_ = 0; k_ = i; }

            f[fid] += alpha[i][fid]
                    * (std::pow(x_tail, i_) * std::pow(y_tail, j_) * std::pow(z_tail, k_));
        }
    }

    return f;
}

std::array<double, 3> bicubic(
    const double &x, const double &y, const double &z, // coords in new grid
    const uint &x_idx, const uint &y_idx, const uint &z_idx, // ids in new grid
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid, // old grids
    const std::array<int, 3> &dec_coeffs,
    const std::array<int, 3> &mask2d, // {1, 1, 0} for x-y bicubic, {1, 0, 1} for x-z, {0, 1, 1} for y-z
    //                   y-z bicubic
    //                  / yes
    // mask2d[0] = 0 ?                 x-z bicubic
    //                  \ no          /yes
    //                   mask2d[1] = 0?
    //                                \ no
    //                                 x-y bicubic
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
) {
    std::array<std::array<double, 3>, 16> b;
    std::array<std::array<double, 3>, 16> alpha;
    std::array<double, 3> f {0., 0., 0.};

    // calculate field derivatives in square + write them to vector b
    //
    //    Y
    //    |
    //    p3---------p4
    //    |          |
    //    |          |
    //    |          |
    //    p1---------p2 ------ X

    for (auto i = 0; i < 4; i++) {
        auto px = x_idx / dec_coeffs[0];
        auto py = y_idx / dec_coeffs[1];
        auto pz = z_idx / dec_coeffs[2];

        if (mask2d[0] == 0) {
            py += (i % 2 == 1) ? 1 : 0;
            pz += (i > 1) ? 1 : 0;
        } else {
            if (mask2d[1] == 0) {
                px += (i % 2 == 1) ? 1 : 0;
                pz += (i > 1) ? 1 : 0;
            } else {
                px += (i % 2 == 1) ? 1 : 0;
                py += (i > 1) ? 1 : 0;
            }
        }

        auto derivs = compute_derivatives(px, py, pz, x_grid, y_grid, z_grid, field);

        b[i] = derivs.f;
        if (mask2d[0] == 0) {
            b[i + 4] = derivs.dfdy;
            b[i + 8] = derivs.dfdz;
            b[i + 12] = derivs.d2fdydz;
        } else {
            if (mask2d[1] == 0) {
                b[i + 4] = derivs.dfdx;
                b[i + 8] = derivs.dfdz;
                b[i + 12] = derivs.d2fdxdz;
            } else {
                b[i + 4] = derivs.dfdx;
                b[i + 8] = derivs.dfdy;
                b[i + 12] = derivs.d2fdxdy;
            }
        }
    }

    // calculate values in vector alpha

    for (auto i = 0; i < 16; i++) {
        alpha[i] = {0, 0, 0};
        for (auto j = 0; j < 16; j++) {
            for (auto fid = 0; fid < 3; fid ++) {
                alpha[i][fid] += bicubic_B[i][j] * b[j][fid];
            } 
        }
    }

    // finally fill f values

    auto x_tail = std::fmod(x, x_grid.step);
    auto y_tail = std::fmod(y, y_grid.step);
    auto z_tail = std::fmod(z, z_grid.step);

    for (auto i = 0; i <= 3; i++) {
        for (auto j = 0; j <= 3; j++) {
            for (auto fid = 0; fid < 3; fid++) {
                auto i_ = 0;
                auto j_ = 0;
                auto k_ = 0;
                if (mask2d[0] == 0) {
                    i_ = 0; j_ = i; k_ = j;
                } else {
                    if (mask2d[1] == 0) {
                        i_ = i; j_ = 0; k_ = j;
                    } else {
                        i_ = i; j_ = j; k_ = 0;
                    }
                }
                f[fid] += alpha[i + 4 * j][fid]
                        * (std::pow(x_tail, i_) * std::pow(y_tail, j_) * std::pow(z_tail, k_));
            }
        }
    }

    return f;
}

std::array<double, 3> tricubic(
    const double &x, const double &y, const double &z, // coords in new grid
    const uint &x_idx, const uint &y_idx, const uint &z_idx, // ids in new grid
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid, // old grids
    const std::array<int, 3> &dec_coeffs,
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
) {
    // for full scheme see F. Lekien et. al. "Tricubic interpolation in three dimensions"

    std::array<std::array<double, 3>, 64> b;
    std::array<std::array<double, 3>, 64> alpha;
    std::array<double, 3> f {0., 0., 0.};

    // calculate field derivatives in cube + write them to vector b
    //
    //    Z
    //    |   p7---------p8
    //    p5---------p6  |
    //    |   |      |   |
    //    |   |   Y  |   |
    //    |   p3-----|---p4
    //    p1---------p2 ------ X

    for (auto i = 0; i < 8; i++) {
        auto px = x_idx / dec_coeffs[0];
        auto py = y_idx / dec_coeffs[1];
        auto pz = z_idx / dec_coeffs[2];

        px += (i % 2 == 1) ? 1 : 0;
        py += ((i % 4 == 0) || (i % 4 == 1)) ? 1 : 0;
        pz += (i > 3) ? 1 : 0;

        auto derivs = compute_derivatives(px, py, pz, x_grid, y_grid, z_grid, field);

        b[i] = derivs.f;
        b[i + 8] = derivs.dfdx;
        b[i + 16] = derivs.dfdy;
        b[i + 24] = derivs.dfdz;
        b[i + 32] = derivs.d2fdxdy;
        b[i + 40] = derivs.d2fdxdz;
        b[i + 48] = derivs.d2fdydz;
        b[i + 56] = derivs.d3fdxdyxz;
    }

    // calculate values in vector alpha

    for (auto i = 0; i < 64; i++) {
        alpha[i] = {0, 0, 0};
        for (auto j = 0; j < 64; j++) {
            for (auto fid = 0; fid < 3; fid ++) {
                alpha[i][fid] += tricubic_B[i][j] * b[j][fid];
            } 
        }
    }

    // finally fill f values

    auto x_tail = std::fmod(x, x_grid.step);
    auto y_tail = std::fmod(y, y_grid.step);
    auto z_tail = std::fmod(z, z_grid.step);

    for (auto i = 0; i <= 3; i++) {
        for (auto j = 0; j <= 3; j++) {
            for (auto k = 0; k <= 3; k++) {
                for (auto fid = 0; fid < 3; fid++) {
                    f[fid] += alpha[i + 4 * j + 16 * k][fid]
                            * (std::pow(x_tail, i) * std::pow(y_tail, j) * std::pow(z_tail, k));
                }
            }
        }
    }

    return f;
}
