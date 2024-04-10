#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iterator>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <typeinfo>
#include <tuple>
#include <vector>

struct Grid {
    uint size;
    float step;
    float min_value;
    float max_value;
    std::vector<float> grid;
};

class MagneticField {
    // Constructor
    public:
        explicit MagneticField(
            const std::string &file_path,
            const bool &interpolate=false,
            const std::array<int, 3> decomp_coeffs={0, 0, 0}
            ) : interpolate(interpolate), decomp_coeffs(decomp_coeffs) {
            this->read_mf(file_path);

            if (this->interpolate) {
                this->compute_interpolated_matrix();
            }
        };

    private:
        void read_mf(const std::string &file_path);
        std::vector<std::string> read_raw(const std::string &file_path);
        Grid define_grid(const int& ax, const std::vector<std::string>& raw_mf);
        void fill_field(const std::vector<std::string>& raw_mf);
        void compute_interpolated_matrix();

        Grid x_grid;
        Grid y_grid;
        Grid z_grid;

        bool interpolate;
        std::array<int, 3> decomp_coeffs;

        // X * Y * Z * (Bx, By, Bz)
        std::vector<std::vector<std::vector<std::array<double, 3>>>> field;

    // Get method
    public:
        std::array<double, 3> get_field(const float& x, const float& y, const float& z);

    private:
        std::tuple<int, bool> lower_id_is_same(const Grid& grid, const float& point);
        // Compute linearly interpolated field
        std::array<double, 3> compute_field(
            const float& x, const float& y, const float& z,
            const std::tuple<int, bool>& x_idxs,
            const std::tuple<int, bool>& y_idxs,
            const std::tuple<int, bool>& z_idxs
            );
        
    // Destructor
    public:
        ~MagneticField() {};
};

struct Derivatives {
    std::array<double, 3> f;

    std::array<double, 3> dfdx;
    std::array<double, 3> dfdy;
    std::array<double, 3> dfdz;
    
    std::array<double, 3> d2fdxdy;
    std::array<double, 3> d2fdxdz;
    std::array<double, 3> d2fdydz;

    std::array<double, 3> d3fdxdyxz;
};

extern std::array<std::array<int, 4>, 4> cubic_B;
extern std::array<std::array<int, 16>, 16> bicubic_B;
extern std::array<std::array<int, 64>, 64> tricubic_B;

std::array<double, 3> interp_1d(
    const float& x, const float& x0, const float& step,
    const std::array<double, 3>& f0, const std::array<double, 3>& f1);

std::array<double, 3> interp_2d(
    const float& x, const float& y, const float& x0, const float& y0,
    const float& step_x, const float& step_y,
    const std::array<double, 3>& f00, const std::array<double, 3>& f01,
    const std::array<double, 3>& f10, const std::array<double, 3>& f11);

std::array<double, 3> interp_3d(
    const float& x, const float& y, const float& z,
    const float& x0, const float& y0, const float& z0,
    const float& step_x, const float& step_y, const float& step_z,
    const std::array<double, 3>& f000, const std::array<double, 3>& f001,
    const std::array<double, 3>& f010, const std::array<double, 3>& f011,
    const std::array<double, 3>& f100, const std::array<double, 3>& f101,
    const std::array<double, 3>& f110, const std::array<double, 3>& f111);

std::array<double, 3> cubic(
    const double &x, const double &y, const double &z,
    const uint &x_idx, const uint &y_idx, const uint &z_idx,
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid,
    const std::array<int, 3> &dec_coeffs,
    const uint &mask1d,
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
);

std::array<double, 3> bicubic(
    const double &x, const double &y, const double &z,
    const uint &x_idx, const uint &y_idx, const uint &z_idx,
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid,
    const std::array<int, 3> &dec_coeffs,
    const std::array<int, 3> &mask2d,
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
);

std::array<double, 3> tricubic(
    const double &x, const double &y, const double &z,
    const uint &x_idx, const uint &y_idx, const uint &z_idx,
    const Grid &x_grid, const Grid &y_grid, const Grid &z_grid,
    const std::array<int, 3> &dec_coeffs,
    const std::vector<std::vector<std::vector<std::array<double, 3>>>> &field
);
