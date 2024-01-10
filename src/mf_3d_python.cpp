#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "mf_3d.h"

PYBIND11_MODULE(mf_3d, m) {
    pybind11::class_<MagneticField>(m, "MagneticField")
        .def(pybind11::init<const std::string &>())
        .def("get_field", &MagneticField::get_field);
}