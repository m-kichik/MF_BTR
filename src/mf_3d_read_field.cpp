#include "mf_3d.h"

// TODO: checks for missed values (e.g. one missed value in string)

void MagneticField::read_mf(const std::string &file_path) {
    std::vector<std::string> mf_raw = this->read_raw(file_path);
    if (!mf_raw.size()) {
        throw std::runtime_error("Empty MF");
    }

    this->x_grid = this->define_grid(0, mf_raw);
    this->y_grid = this->define_grid(1, mf_raw);
    this->z_grid = this->define_grid(2, mf_raw);

    this->fill_field(mf_raw);
};

std::vector<std::string> MagneticField::read_raw(const std::string &file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::invalid_argument("Can't open MF file");
    }

    std::vector<std::string> mf_raw;

    std::string num_pattern = R"([+-]?\d+(\.\d+)?[Ee]?[+\-]?(\d+)?)";
    std::string mf_pattern = R"(\s*)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s*)";
    std::regex mf_regex_pattern(mf_pattern);

    std::regex empty_str_regex_pattern(R"(\s*)");

    auto reading_mf = false;
    auto line_idx = 1;
    std::string line;
    while (std::getline(file, line)) {
        auto match = std::regex_match(line, mf_regex_pattern);

        if (!reading_mf && match) {
            reading_mf = true;
        }

        if (reading_mf && match) {
            std::istringstream values(line);
            std::string value;
            while (values >> value) {
                mf_raw.push_back(value);
            }
        }

        if (reading_mf && !match) {
            if (!std::regex_match(line, empty_str_regex_pattern)) {
                auto exception_msg = "Error in MF file content: line " +
                                     std::to_string(line_idx) +
                                     " doesn't satisfy MF 3D format pattern";
                throw std::invalid_argument(exception_msg);
            }
        }

        line_idx++;
    }

    file.close();

    return mf_raw;
}

Grid MagneticField::define_grid(const int &ax, const std::vector<std::string> &raw_mf) {
    std::set<float> points;
    std::vector<float> grid;
    for (auto i = ax; i < raw_mf.size(); i += 6) {
        points.insert(std::stof(raw_mf[i]));
    }

    std::copy(points.begin(), points.end(), std::back_inserter(grid));

    auto n_nodes = points.size();

    auto min_value = *points.begin();
    auto max_value = *points.rbegin();

    auto step = (max_value - min_value) / (n_nodes - 1);

    // TODO: check that step remains the same

    return Grid{(uint)(n_nodes), step, min_value, max_value, grid};
};

void MagneticField::fill_field(const std::vector<std::string> &raw_mf) {
    this->field.resize(
        this->x_grid.size,
        std::vector<std::vector<std::array<double, 3>>>(
            this->y_grid.size,
            std::vector<std::array<double, 3>>(
                this->z_grid.size,
                std::array<double, 3>{1e20, 1e20, 1e20})));

    for (auto l = 0; l < raw_mf.size(); l += 6) {
        auto i = std::distance(
            this->x_grid.grid.begin(),
            std::find(
                this->x_grid.grid.begin(),
                this->x_grid.grid.end(),
                std::stof(raw_mf[l])));

        auto j = std::distance(
            this->y_grid.grid.begin(),
            std::find(
                this->y_grid.grid.begin(),
                this->y_grid.grid.end(),
                std::stof(raw_mf[l + 1])));

        auto k = std::distance(
            this->z_grid.grid.begin(),
            std::find(
                this->z_grid.grid.begin(),
                this->z_grid.grid.end(),
                std::stof(raw_mf[l + 2])));

        this->field[i][j][k][0] = std::stod(raw_mf[l + 3]);
        this->field[i][j][k][1] = std::stod(raw_mf[l + 4]);
        this->field[i][j][k][2] = std::stod(raw_mf[l + 5]);
    }

    // Final check that every value was set, this doesn't save from missed values that cause abortion.
    for (auto i = 0; i < this->x_grid.size; i++) {
        for (auto j = 0; j < this->y_grid.size; j++) {
            for (auto k = 0; k < this->z_grid.size; k++) {
                auto s = std::accumulate(this->field[i][j][k].begin(), this->field[i][j][k].end(), 0.);
                if (s > 1e20) {
                    std::string coord_str = std::to_string(this->x_grid.grid[i]) + ", " +
                                            std::to_string(this->y_grid.grid[j]) + ", " +
                                            std::to_string(this->z_grid.grid[k]);
                    throw std::invalid_argument(
                        "Some value(s) was/were not set in the coordinate " + coord_str + ".");
                }
            }
        }
    }
}
