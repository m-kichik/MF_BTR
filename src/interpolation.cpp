#include "Interpolation.hpp"


/**
* Performs linear interpolation in 1D.
*
* @param f0 The value at the left endpoint.
* @param f1 The value at the right endpoint.
* @param delta The distance from the left endpoint.
* @return The interpolated value.
*/
std::array<double, 3> linear_1d(
    std::array<double, 3> f0,
    std::array<double, 3> f1,
    double delta
) {
    std::array<double, 3> f;
    for (int i = 0; i < 3; i++) {
        f[i] = (f0[i] + delta * (f1[i] - f0[i]));
    }
    return f;
}


/**
* Performs linear interpolation in 2D.
*
* @param f00 The value at the bottom-left corner.
* @param f01 The value at the bottom-right corner.
* @param f10 The value at the top-left corner.
* @param f11 The value at the top-right corner.
* @param delta0 The distance from the bottom-left corner.
* @param delta1 The distance from the top-left corner.
* @return The interpolated value.
*/
std::array<double, 3> linear_2d(
    std::array<double, 3> f00,
    std::array<double, 3> f01,
    std::array<double, 3> f10,
    std::array<double, 3> f11,
    double delta0, double delta1
) {
        std::array<double, 3> f;
        for (int i = 0; i < 3; i++) {
            f[i] = (1 - delta1) * ((1 - delta0) * f00[i] + delta0 * f10[i]) + 
                    delta1 * ((1 - delta0) * f01[i] + delta0 * f11[i]);
        }
        return f;
};


/**
* Performs linear interpolation in 3D.
*
* @param f000 The value at the bottom-left-front corner.
* @param f001 The value at the bottom-left-back corner.
* @param f010 The value at the bottom-right-front corner.
* @param f011 The value at the bottom-right-back corner.
* @param f100 The value at the top-left-front corner.
* @param f101 The value at the top-left-back corner.
* @param f110 The value at the top-right-front corner.
* @param f111 The value at the top-right-back corner.
* @param delta0 The distance from the bottom-left-front corner.
* @param delta1 The distance from the bottom-left-back corner.
* @param delta2 The distance from the top-left-front corner.
* @return The interpolated value.
*/
std::array<double, 3> linear_3d(
    std::array<double, 3> f000,
    std::array<double, 3> f001,
    std::array<double, 3> f010,
    std::array<double, 3> f011,
    std::array<double, 3> f100,
    std::array<double, 3> f101,
    std::array<double, 3> f110,
    std::array<double, 3> f111,
    double delta0,
    double delta1,
    double delta2
) {
    std::array<double, 3> f00;
    std::array<double, 3> f01;
    std::array<double, 3> f10;
    std::array<double, 3> f11;

    for (auto i = 0u; i < 3; i++) {
        f00[i] = (1 - delta0) * f000[i] + delta0 * f100[i];
        f01[i] = (1 - delta0) * f001[i] + delta0 * f101[i];
        f10[i] = (1 - delta0) * f010[i] + delta0 * f110[i];
        f11[i] = (1 - delta0) * f011[i] + delta0 * f111[i];
    }

    std::array<double, 3> f0;
    std::array<double, 3> f1;

    for (auto i = 0u; i < 3; i++) {
        f0[i] = (1 - delta1) * f00[i] + delta1 * f001[i];
        f1[i] = (1 - delta1) * f01[i] + delta1 * f11[i];
    }

    std::array<double, 3> f;

    for (auto i = 0u; i < 3; i++) {
        f[i] = (1 - delta2) * f0[i] + delta2 * f1[i];
    }

    return f;
}
