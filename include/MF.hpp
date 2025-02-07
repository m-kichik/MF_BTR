// #include <algorithm>
#include <array>
// #include <chrono>
#include <cmath>
// #include <format>
#include <fstream>
// #include <iterator>
#include <iostream>
#include <limits>
// #include <numeric>
#include <regex>
#include <set>
#include <sstream>
// #include <stdexcept>
// #include <string>
// #include <thread>
// #include <typeinfo>
#include <tuple>
#include <vector>


struct Grid {
    uint64_t size;
    double min, max, step;
    std::vector<double> grid;
};

struct Derivatives {
    std::array<double, 3> f;

    std::array<double, 3> dfdx;
    std::array<double, 3> dfdy;
    std::array<double, 3> dfdz;
    
    std::array<double, 3> d2fdxdy;
    std::array<double, 3> d2fdxdz;
    std::array<double, 3> d2fdydz;

    std::array<double, 3> d3fdxdydz;
};


class MagneticField {
    public:
        explicit MagneticField(
            const std::string &file_path,
            const bool &interpolate=false,
            const std::array<int, 3> scaling_coeffs={1, 1, 1}
            ) : interpolate(interpolate), scaling_coeffs(scaling_coeffs) {
            this->read_mf(file_path);

            if (interpolate) {
                this->interpolate_field();
            }
        };

    private:
        using raw_mf_and_corrupted_lines = std::tuple<
        std::vector<std::vector<double>>, std::vector<uint>>;
        using field_type = std::vector<
        std::vector<std::vector<std::array<double, 3>>>>;

        void read_mf(const std::string &file_path);
        raw_mf_and_corrupted_lines read_raw(const std::string &file_path);
        Grid x_grid, y_grid, z_grid;
        std::tuple<bool, bool, bool> define_grid(
            std::vector<std::vector<double>> raw_mf
            );
        std::vector<std::vector<double>> fill_field(std::vector<std::vector<double>> raw_mf);

        field_type field;

        bool interpolate;
        std::array<int, 3> scaling_coeffs;
        void interpolate_field();
        void interpolate_in_cell(
            std::tuple<uint64_t, uint64_t, uint64_t> coords,
            field_type &field,
            Grid &x_grid, Grid &y_grid, Grid &z_grid
            );

    public:
        std::array<double, 3> get_field(const std::tuple<double, double, double> &coords);

    public:
        ~MagneticField() {};
};