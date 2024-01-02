cmake_minimum_required(VERSION 3.5)

project(BTR_Field_3D VERSION 0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(extern/pybind11)

set(PROJECT_SOURCES
    main.cpp

    mf_3d.h
    mf_3d_read_field.cpp
    mf_3d_get.cpp
    mf_3d_interpolation.cpp
)

pybind11_add_module(
    mf_3d

    mf_3d_python.cpp
    
    mf_3d.h
    mf_3d_read_field.cpp
    mf_3d_get.cpp
    mf_3d_interpolation.cpp
    )

add_executable(BTR_Field_3D ${PROJECT_SOURCES})
