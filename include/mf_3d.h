#include <algorithm>
#include <cmath>
// #include <format>
#include <fstream>
#include <iterator>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
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
        explicit MagneticField(const std::string &file_path){
            this->read_mf(file_path);
        };

    private:
        void read_mf(const std::string &file_path);
        std::vector<std::string> read_raw(const std::string &file_path);
        Grid define_grid(const int& ax, const std::vector<std::string>& raw_mf);
        void fill_field(const std::vector<std::string>& raw_mf);

        Grid x_grid;
        Grid y_grid;
        Grid z_grid;

        // X * Y * Z * (Bx, By, Bz)
        std::vector<std::vector<std::vector<std::array<double, 3>>>> field;

    // Get method
    public:
        std::array<double, 3> get_field(const float& x, const float& y, const float& z);

    private:
        // Returns lower idx + true if point in grid node, false else
        std::tuple<int, bool> lower_id_is_same(const Grid& grid, const float& point);
        bool is_in_node(const float& point, const float& step, const float& epsilon);
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
