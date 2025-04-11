#include <array>
#include <sstream>
#include <tuple>
#include <vector>


struct Grid {
    Grid(
        uint64_t size = 0, 
        double min = 0, 
        double max = 0, 
        double step = 0, 
        std::vector<double> grid={}
    ) : size(size), 
        min(min), 
        max(max), 
        step(step), 
        grid(grid) 
    {}
        
    uint64_t size;
    double min, max, step;
    std::vector<double> grid;
};


struct RawMFData {
    std::vector<std::vector<double>> field_data;
    std::vector<uint> corrupted_lines;
};


using FieldCell = std::array<double, 3>;
using FieldLayer = std::vector<std::vector<FieldCell>>;
using Field3D = std::vector<FieldLayer>;


class MagneticField {
    public:
        explicit MagneticField(
            std::string file_path,
            bool interpolate = false,
            std::array<int, 3> scaling_coeffs={1, 1, 1}
        ) : interpolate_(interpolate), 
            scaling_coeffs_(std::move(scaling_coeffs)) 
        {
            read_mf(file_path);

            if (interpolate) {
                this->interpolate_field();
            }
        };

    private:
        void read_mf(const std::string &file_path);
        RawMFData read_raw(const std::string &file_path);
        Grid x_grid_, y_grid_, z_grid_;
        std::tuple<bool, bool, bool> define_grid(
            std::vector<std::vector<double>> raw_mf
            );
        std::vector<std::vector<double>> fill_field(std::vector<std::vector<double>> raw_mf);

        Field3D field_;

        bool interpolate_;
        std::array<int, 3> scaling_coeffs_;
        void interpolate_field();
        void interpolate_in_cell(
            std::tuple<uint64_t, uint64_t, uint64_t> coords,
            Field3D &field,
            Grid &x_grid, Grid &y_grid, Grid &z_grid
            );

    public:
        std::array<double, 3> get_field(const std::tuple<double, double, double> &coords);

    public:
        ~MagneticField() = default;
        MagneticField(const MagneticField&) = delete;
        MagneticField& operator=(const MagneticField&) = delete;
        MagneticField(MagneticField&&) = default;
        MagneticField& operator=(MagneticField&&) = default;
};
