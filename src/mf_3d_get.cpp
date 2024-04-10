#include "mf_3d.h"

std::array<double, 3> MagneticField::get_field(const float& x, const float& y, const float& z) {
    auto x_ids = this->lower_id_is_same(this->x_grid, x);
    auto y_ids = this->lower_id_is_same(this->y_grid, y);
    auto z_ids = this->lower_id_is_same(this->z_grid, z);
    
    return compute_field(x, y, z, x_ids, y_ids, z_ids);
};
