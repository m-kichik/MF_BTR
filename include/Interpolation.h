#include <vector>

std::array<double, 3> linear_1d(
    std::array<double, 3> f0, 
    std::array<double, 3> f1, 
    double delta
);

std::array<double, 3> linear_2d(
    std::array<double, 3> f00,
    std::array<double, 3> f01,
    std::array<double, 3> f10,
    std::array<double, 3> f11,
    double delta0, double delta1
);

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
);
