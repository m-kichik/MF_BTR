cmake_minimum_required(VERSION 3.5)

project(BTR_Field_3D VERSION 0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(extern/pybind11)

# set(LIBRARY_NAME
#     MF
# )

include_directories(
    include
)

set(PROJECT_SOURCES
    src/main.cpp

    include/MF.h
    include/Interpolation.h

    src/read_field.cpp
    src/get_field.cpp
    src/interpolation.cpp
)

pybind11_add_module(
    mf_3d

    src/mf_python.cpp
    
    include/MF.h
    include/Interpolation.h

    src/read_field.cpp
    src/get_field.cpp
    src/interpolation.cpp
    )

add_executable(BTR_Field_3D ${PROJECT_SOURCES})
