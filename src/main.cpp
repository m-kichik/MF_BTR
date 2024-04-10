#include <cmath>
#include "mf_3d.h"

int main(int argc, char *argv[]) {
    auto mf = MagneticField(argv[1], true, {2, 2, 2});
    // auto mf = MagneticField(argv[1]);
    auto field = mf.get_field(0., 0., 0.35);
    std::cout << field[0] << " " << field[1] << " " << field[2] << std::endl;

    return EXIT_SUCCESS;
}