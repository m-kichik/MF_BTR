#include "mf_3d.h"

std::array<double, 3> interp_1d_linear(
    const float& x, const float& x0, const float& step,
    const std::array<double, 3>& f0, const std::array<double, 3>& f1);

std::array<double, 3> interp_1d(
    const float& x, const float& x0, const float& step,
    const std::array<double, 3>& f0, const std::array<double, 3>& f1
) {
    return interp_1d_linear(x, x0, step, f0, f1);
}

std::array<double, 3> interp_2d_linear(
    const float& x, const float& y, const float& x0, const float& y0,
    const float& step_x, const float& step_y,
    const std::array<double, 3>& f00, const std::array<double, 3>& f01,
    const std::array<double, 3>& f10, const std::array<double, 3>& f11);

std::array<double, 3> interp_2d(
    const float& x, const float& y, const float& x0, const float& y0,
    const float& step_x, const float& step_y,
    const std::array<double, 3>& f00, const std::array<double, 3>& f01,
    const std::array<double, 3>& f10, const std::array<double, 3>& f11
) {
    return interp_2d_linear(x, y, x0, y0, step_x, step_y, f00, f01, f10, f11);
}

std::array<double, 3> interp_3d_linear(
    const float& x, const float& y, const float& z,
    const float& x0, const float& y0, const float& z0,
    const float& step_x, const float& step_y, const float& step_z,
    const std::array<double, 3>& f000, const std::array<double, 3>& f001,
    const std::array<double, 3>& f010, const std::array<double, 3>& f011,
    const std::array<double, 3>& f100, const std::array<double, 3>& f101,
    const std::array<double, 3>& f110, const std::array<double, 3>& f111);

std::array<double, 3> interp_3d(
    const float& x, const float& y, const float& z,
    const float& x0, const float& y0, const float& z0,
    const float& step_x, const float& step_y, const float& step_z,
    const std::array<double, 3>& f000, const std::array<double, 3>& f001,
    const std::array<double, 3>& f010, const std::array<double, 3>& f011,
    const std::array<double, 3>& f100, const std::array<double, 3>& f101,
    const std::array<double, 3>& f110, const std::array<double, 3>& f111
) {
    return interp_3d_linear(
        x, y, z, x0, y0, z0, step_x, step_y, step_z,
        f000, f001, f010, f011, f100, f101, f110, f111
    );
}

std::array<double, 3> interp_1d_linear(
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

std::array<double, 3> interp_2d_linear(
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

std::array<double, 3> interp_3d_linear(
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
