#include "MF.hpp"

void MagneticField::read_mf(const std::string &file_path) {
    auto [raw_mf, corrupted_lines] = this->read_raw(file_path);

    if (corrupted_lines.size()) {
        std::cerr << "Warning: following lines may be corrupted and were not parsed: ";
        for (auto i = 0; i < corrupted_lines.size(); i++) {
            std::cerr <<  std::to_string(corrupted_lines[i]);
            if (i!= corrupted_lines.size() - 1) {
                std::cerr << ", ";
            }
        }
        std::cerr << std::endl;
    }

    if (!raw_mf.size()) {
        throw std::runtime_error("Empty MF file");
    }

    auto [x_grid_ok, y_grid_ok, z_grid_ok] = this->define_grid(raw_mf);
    if (! x_grid_ok) {
        std::cerr << 
        "Warning: X grid has non-constant steps. A part lying between " << 
        this->x_grid.min << " and " << this->x_grid.max << " was parsed." << std::endl;
    }
    if (! y_grid_ok) {
        std::cerr << 
        "Warning: Y grid has non-constant steps. A part lying between " << 
        this->y_grid.min << " and " << this->y_grid.max << " was parsed." << std::endl;
    }
    if (! z_grid_ok) {
        std::cerr << 
        "Warning: Z grid has non-constant steps. A part lying between " << 
        this->z_grid.min << " and " << this->z_grid.max << " was parsed." << std::endl;
    }

    auto missed_points = this->fill_field(raw_mf);
    if (missed_points.size()) {
        std::cerr << "Warning: following points are missed and magnetic field in these points is set to zero: ";
        for (auto i = 0; i < missed_points.size(); i++) {
            auto point = missed_points[i];
            for (auto &v: point) {
                if (std::abs(v) < 1e-8) {v = 0.;}
            }
            std::cout << "(" << point[0] << " " << point[1] << " " << point[2] << ")";
            if (i != missed_points.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
}

MagneticField::raw_mf_and_corrupted_lines MagneticField::read_raw(const std::string &file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Can't open MF file");
    }

    std::vector<std::vector<double>> raw_mf;
    std::string line;

    std::string num_pattern = R"([+-]?\d+(\.\d+)?[Ee]?[+\-]?(\d+)?)";
    std::string mf_pattern = R"(\s*)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s+)" +
                          num_pattern + R"(\s*)";
    std::regex mf_regex_pattern(mf_pattern);
    std::regex empty_string_pattern("");

    auto line_counter = 0u;
    auto read_mf = false;
    std::vector<uint> corrupted_lines;

    while (std::getline(file, line)) {
        auto match = std::regex_match(line, mf_regex_pattern);
        if (!read_mf && match) {
            read_mf = true;
        }
        if (match) {
            std::vector<double> point;
            std::istringstream values(line);
            std::string value;
            while (values >> value) {
                point.push_back(std::stod(value));
            }
            if (point.size() != 6) {
                corrupted_lines.push_back(line_counter);
            } else {
                raw_mf.push_back(point);
            }
        } else {
            auto empty_match = std::regex_match(line, empty_string_pattern);
            if (read_mf && !empty_match) {
                corrupted_lines.push_back(line_counter);
            }
        }

        line_counter += 1;
    }

    return std::make_tuple(raw_mf, corrupted_lines);
}

std::tuple<bool, std::vector<double>> check_constant_step(std::vector<double> grid) {
    if (grid.size() <= 1) {
        return std::make_tuple(true, grid);
    }

    auto step = grid[1] - grid[0];
    std::vector<double> verified_grid {grid.front()};

    for (auto i = 0; i < grid.size() - 1; i++) {
        if (std::abs((grid[i + 1] - grid[i]) - step) < 1e-10) {
            verified_grid.push_back(grid[i + 1]);
        } else {
            return std::make_tuple(false, verified_grid);
        }
    }
    
    return std::make_tuple(true, verified_grid);
}

std::tuple<bool, bool, bool> MagneticField::define_grid(std::vector<std::vector<double>> raw_mf) {
    std::set<double> x_grid_set, z_grid_set, y_grid_set;

    for (auto &point: raw_mf) {
        for (auto &value: point) {
            x_grid_set.insert(point[0]);
            y_grid_set.insert(point[1]);
            z_grid_set.insert(point[2]);
        }
    }

    auto x_grid_values = std::vector<double>(x_grid_set.begin(), x_grid_set.end());
    auto y_grid_values = std::vector<double>(y_grid_set.begin(), y_grid_set.end());
    auto z_grid_values = std::vector<double>(z_grid_set.begin(), z_grid_set.end());

    auto [x_grid_ok, verified_x_grid_values] = check_constant_step(x_grid_values);
    auto [y_grid_ok, verified_y_grid_values] = check_constant_step(y_grid_values);
    auto [z_grid_ok, verified_z_grid_values] = check_constant_step(z_grid_values);

    this->x_grid.min = verified_x_grid_values.front();
    this->x_grid.max = verified_x_grid_values.back();
    if (verified_x_grid_values.size() > 1) {
        this->x_grid.step = verified_x_grid_values[1] - verified_x_grid_values[0];
    } else {
        this->x_grid.step = 0;
    }
    this->x_grid.grid = verified_x_grid_values;
    this->x_grid.size = verified_x_grid_values.size();

    this->y_grid.min = verified_y_grid_values.front();
    this->y_grid.max = verified_y_grid_values.back();
    if (verified_y_grid_values.size() > 1) {
        this->y_grid.step = verified_y_grid_values[1] - verified_y_grid_values[0];
    } else {
        this->y_grid.step = 0;
    }
    this->y_grid.grid = verified_y_grid_values;
    this->y_grid.size = verified_y_grid_values.size();

    this->z_grid.min = verified_z_grid_values.front();
    this->z_grid.max = verified_z_grid_values.back();
    if (verified_z_grid_values.size() > 1) {
        this->z_grid.step = verified_z_grid_values[1] - verified_z_grid_values[0];
    } else {
        this->z_grid.step = 0;
    }
    this->z_grid.grid = verified_z_grid_values;
    this->z_grid.size = verified_z_grid_values.size();

    return std::make_tuple(x_grid_ok, y_grid_ok, z_grid_ok);
}

std::vector<std::vector<double>> MagneticField::fill_field(std::vector<std::vector<double>> raw_mf) {
    auto raw_mf_max = [](std::vector<double> point_1, std::vector<double> point_2) {
        auto max_1 = std::max({std::abs(point_1[3]), std::abs(point_1[4]), std::abs(point_1[5])});
        auto max_2 = std::max({std::abs(point_2[3]), std::abs(point_2[4]), std::abs(point_2[5])});

        return max_1 < max_2;
    };

    // Define big value to fill field matrix and thus check for missings after filling values
    auto max_point = *std::max_element(raw_mf.begin(), raw_mf.end(), raw_mf_max);
    auto big_value = (
        std::abs(max_point[3]) + std::abs(max_point[4]) + std::abs(max_point[5])
        ) * 1000 + 1000;

    // Resize field matrix
    this->field.resize(this->x_grid.size);
    for (auto &row: this->field) {
        row.resize(this->y_grid.size);
        for (auto &col: row) {
            col.resize(this->z_grid.size);
            for (auto &point: col) {
                point = std::array<double, 3>{big_value, big_value, big_value};
            }
        }
    }

    // Fill field matrix
    for (auto &point: raw_mf) {
        auto x_coord = point[0];
        auto y_coord = point[1];
        auto z_coord = point[2];

        if (
            ((this->x_grid.min <= x_coord) && (x_coord <= this->x_grid.max)) 
            && ((this->y_grid.min <= y_coord) && (y_coord <= this->y_grid.max)) 
            && ((this->z_grid.min <= z_coord) && (z_coord <= this->z_grid.max))
            ) {
            auto x_index = std::round((x_coord - this->x_grid.min) / this->x_grid.step);
            auto y_index = std::round((y_coord - this->y_grid.min) / this->y_grid.step);
            auto z_index = std::round((z_coord - this->z_grid.min) / this->z_grid.step);

            this->field[x_index][y_index][z_index] = {point[3], point[4], point[5]};
        }
    }

    // Check for missings in field matrix and fill them with zeros
    std::vector<std::vector<double>> missed_points;
    for (auto i = 0; i < this->x_grid.size; i++) {
        for (auto j = 0; j < this->y_grid.size; j++) {
            for (auto k = 0; k < this->z_grid.size; k++) {
                if (this->field[i][j][k][0] + this->field[i][j][k][1] + this->field[i][j][k][2] >= big_value) {
                    auto x_coord = (i * this->x_grid.step) + this->x_grid.min;
                    auto y_coord = (j * this->y_grid.step) + this->y_grid.min;
                    auto z_coord = (k * this->z_grid.step) + this->z_grid.min;

                    missed_points.push_back({x_coord, y_coord, z_coord});
                    this->field[i][j][k] = {0., 0., 0.};
                }
            }
        }
    }

    return missed_points;
}
