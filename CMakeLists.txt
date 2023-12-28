cmake_minimum_required(VERSION 3.5)

project(BTR_Field_3D VERSION 0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# find_package(pybind11 REQUIRED)

# pybind11_add_module(mf_3d_python mf_3d_python.cpp)

set(PROJECT_SOURCES
    main.cpp

    mf_3d.h
    mf_3d_read_field.cpp
    mf_3d_get.cpp
    mf_3d_interpolation.cpp
    mf_3d_python.cpp
)

add_executable(BTR_Field_3D ${PROJECT_SOURCES})
