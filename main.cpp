#include "mf_3d.h"

int main() {
    auto mf = MagneticField("/home/rito4ka/dev/magnetic_field_BTR/Basic_PMS.txt");
    auto field = mf.get_field(0.25, -0.15, 0.42);
    std::cout << field[0] << " " << field[1] << " " << field[2] << std::endl;

    return EXIT_SUCCESS;
}