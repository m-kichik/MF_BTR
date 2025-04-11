#include <cmath>
#include <ctime>
#include <iostream>

#include "mf/MF.hpp"

int main(int argc, char *argv[]) {
    std::time_t start = std::time(nullptr);
    auto mf = MagneticField(argv[1], true, {10, 10, 10});
    
    auto b1 = mf.get_field({0., 0., 0.0});
    for (auto &v : b1) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    auto b2 = mf.get_field({0.5, 0.5, 0.5});
    for (auto &v : b2) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    auto b3 = mf.get_field({1, 1, 1});
    for (auto &v : b3) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}