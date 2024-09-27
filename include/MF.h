// #include <algorithm>
#include <array>
// #include <chrono>
#include <cmath>
// #include <format>
#include <fstream>
// #include <iterator>
#include <iostream>
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
    uint size;
    double min, max, step;
    std::vector<double> grid;
};

class MagneticField {
    public:
        explicit MagneticField(
            const std::string &file_path,
            const bool &interpolate=false,
            const std::array<int, 3> scaling_coeffs={1, 1, 1}
            ) : interpolate(interpolate), scaling_coeffs(scaling_coeffs) {
            this->read_mf(file_path);
        };

    private:
        using raw_mf_and_corrupted_lines = std::tuple<std::vector<std::vector<double>>, std::vector<uint>>;

        void read_mf(const std::string &file_path);
        raw_mf_and_corrupted_lines read_raw(const std::string &file_path);
        Grid x_grid, y_grid, z_grid;
        std::tuple<bool, bool, bool> define_grid(
            std::vector<std::vector<double>> raw_mf
            );
        std::vector<std::vector<double>> fill_field(std::vector<std::vector<double>> raw_mf);

        std::vector<std::vector<std::vector<std::array<double, 3>>>> field;

        bool interpolate;
        std::array<int, 3> scaling_coeffs;

    public:
        std::array<double, 3> get_field(const std::tuple<double, double, double> &coords);

    public:
        ~MagneticField() {};
};